# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## 项目概览

PEC 导体球单站 RCS 计算 — 矩量法 (MoM) + RWG 基函数 + BiCGSTAB + FMM 加速。Makefile 构建，C++11，无外部依赖。

教学演示项目，参数硬编码在 `main.cpp` 中，非精确工程计算。

## 构建与运行

```bash
make                    # 编译 sphere_rcs.exe
make run                # 编译并运行
make clean              # 清理 .o 和 .exe
```

## 计算流水线

1. **网格生成** — 正二十面体 → 递归细分 → 归一化投影到球面（Icosphere），参数 `refinement=2` 产生约 320 个三角形
2. **RWG 基函数** — 遍历三角形边，找相邻三角对构成公共边，每条边对应一个 RWG 基函数（`N ≈ 1.5 × 三角形数`）
3. **EFIE 填充** — `Z_{ij} = ∫∫ [f_i·f_j - (1/k²)(∇·f_i)(∇·f_j)] G(R) dS' dS`，用 7 点 Gauss 双重积分
4. **求解 Z·I = b** — BiCGSTAB（复非对称 EFIE 阻抗阵）或 CG（Hermitian 正定阵）
5. **远场/RCS** — 由电流系数 I 通过辐射矢量 N(θ,φ) 计算散射远场和单站 RCS

## 关键类型

| 类型 | 定义位置 | 用途 |
|---|---|---|
| `Vertex` (≡ `Vec3`) | `mesh.h` | 三维点/向量 `{x,y,z}` |
| `Triangle` | `mesh.h` | 三个顶点索引 + 编号，`.area()` `.center()` `.normal()` |
| `Mesh` | `mesh.h` | `vertices` + `triangles`，`.generate_sphere(radius, refinement)` |
| `Edge` | `rwg.h` | RWG 基函数的几何定义：两个相邻三角形、边长、自由顶点 |
| `RWGData` | `rwg.h` | 所有公共边 + 预计算系数 `l/(2A±)` |
| `EFIE` | `efie.h` | 阻抗矩阵填充 + 激励向量填充 |
| `Complex` | `constants.h` | `std::complex<double>` |

## 物理常量 (`constants.h`)

`PI`, `EPS0`, `MU0`, `Z0 = sqrt(MU0/EPS0) ≈ 376.73 Ω`

## EFIE 奇异积分的处理方式

同一三角形上的自作用 (`same_tri == true`) 将 Green 函数拆分为平滑部分 `(e^{-jkR}-1)/(4πR)` 和静态奇异部分 `1/(4πR)`：
- **平滑部分** — 7 点 Gauss 双重积分（`G_smooth` 在 R→0 处解析，Gauss 积分安全）
- **静态奇异部分** — 利用恒等式 `1/R = ∇'²_s R`（平坦三角形）和表面散度定理，将内层面积分 `∫_T 1/R dS'` 转化为三条边上的解析线积分（公式 `Σ P_i^0·f_i`），外层 r 积分用 7 点 Gauss 求 `I_double = ∫∫ 1/R dS' dS`。核函数 `f_i·f_j - ∇·f_i∇·f_j/k²` 在三角形中心点求值，乘以 `jkη·I_double/(4π)`。

关键函数（`efie.cpp`）：
- `integrate_1_over_R(T, verts, r)` — 解析计算 `∫_T 1/|r-r'| dS'`，每条边贡献 `P_i^0·ln((R_A+R_B+L)/(R_A+R_B-L))`
- 静态修正 = `jkη/(4π) * K(center) * [A * Σ_m w_m * integrate_1_over_R(r_m)]`

相邻三角形（共边/共顶点）不做特殊处理，直接 Gauss 积分。

## FMM 加速 (`fmm.h/cpp`, `fmm_translator.h/cpp`)

单层快速多极子 (Single-Level FMM)，使用对角形式（平面波展开）加速矩阵-向量积。

### 核心结构

| 类型 | 定义位置 | 用途 |
|---|---|---|
| `FMM` | `fmm.h` | FMM 主类：网格构建、方向图预计算、近场填充、matvec |
| `FMMBox` | `fmm.h` | 3D 网格盒子：中心、基函数列表、近/远场邻居 |
| `Translator` | `fmm_translator.h` | 对角转移算子 T_L + 球面 Gauss-Legendre 积分 |
| `PWSDirection` | `fmm_translator.h` | 角谱采样方向 `{k_hat, weight}` |
| `Vec3` 运算符 | `vec3_ops.h` | `Vertex` 的 `+`, `-`, `*`, `dot`, `cross`, `norm`, `normalize` |

### FMM 算法流程

1. **空间分组** — 3D 规则网格，盒子尺寸 0.3λ，基函数按边中点分配到盒子
2. **近场预计算** — ±1 邻域盒子间用精确 EFIE 双重积分 (`EFIE::compute_interaction`)，CSR 稀疏存储
3. **方向图预计算** — 4 通道 (矢量 x,y,z + 散度) × K 方向 (Gauss-Legendre 球面积分)
4. **每次迭代 matvec**:
   - **近场**: CSR 稀疏 matvec (O(N_near))
   - **聚合**: 每个盒子对 K 方向求和 4 通道源贡献 (O(N·K))
   - **转移**: 远场盒子对间施加预计算 T_L 算子 (O(P·K), P=远场盒子对数)
   - **解聚**: 测试函数处叠加各方向入射场 (O(N·K))

### EFIE 核的 4 通道分解

Z_ij = jkη [f_i·f_j - k⁻² div_i div_j] G(R)

对角形式将 G(R) 的远场部分分解为 4 个独立标量通道，每个基函数由 (f_x, f_y, f_z, div/k0) 表征。方向图用三角形重心处的 RWG 中心点近似。

### 关键参数

- `box_size`: 盒子尺寸 (~0.3λ)，影响近场/远场划分
- `L`: 截断数，由多余带宽公式 `L ≈ kD/2 + 1.8(kD/2)^{⅓}` 确定
- `K ≈ 2(L+1)²`: 单位球面采样方向数
- 近场覆盖率: ~34% (盒子尺寸 0.3λ 时)

### 精度

- matvec 相对误差 < 2% vs 稠密
- 电流解相对误差 < 2%
- RCS 差异 < 0.1 dB

## 已知限制

- 所有参数硬编码在 `main.cpp`，修改需重新编译
- FMM 远场用中心点近似（三角形重心处求 RWG 值），非完整双面积分
- 单层 FMM，非多层 MLFMA，大规模问题 (N>5000) 效率下降
- 相邻三角形的弱奇异性未单独处理
- 仅支持封闭导体球面网格（`t2 == -1` 的边界边会被跳过）
