# 0511_mm_1 核心模块总结

PEC 导体球单站 RCS 计算 — 矩量法 (MoM) + RWG 基函数 + BiCGSTAB 求解。

---

## 一、EFIE 模块 (`efie.h` / `efie.cpp`)

### 1.1 模块职责

构建电场积分方程 (EFIE) 的复稠密线性系统：**Z·I = b**

- **阻抗矩阵 Z**：N×N 复矩阵，`Z[i][j]` 表示第 j 个 RWG 基函数对第 i 个基函数的电磁耦合
- **激励向量 b**：N 维复向量，`b[i]` 表示入射平面波在第 i 个基函数上的投影

### 1.2 阻抗矩阵填充 — `fill_matrix()`

双重循环遍历所有边对 `(i, j)`，调用 `compute_interaction(i, j)` 逐元素填充。

**核心积分公式**（`integrateTrianglePair`）：

$$
Z_{ij}^{pq} = jk\eta \iint_{T_i^p} \iint_{T_j^q} \left[ \mathbf{f}_i(\mathbf{r}) \cdot \mathbf{f}_j(\mathbf{r}') - \frac{1}{k^2} (\nabla \cdot \mathbf{f}_i)(\nabla \cdot \mathbf{f}_j) \right] G(R) \, dS' \, dS
$$

其中：
- $\mathbf{f}_i(\mathbf{r}) = c_i^{\pm} \cdot (\mathbf{r} - \mathbf{r}_{\text{apex}})$ — RWG 矢量基函数
- $\nabla \cdot \mathbf{f}_i = \pm l / A^{\pm}$ — 散度为常量
- $G(R) = e^{-jkR} / (4\pi R)$ — 自由空间 Green 函数
- 积分用 **7 点 Gauss 双重积分**（49 个采样点/三角形对），5 阶代数精度 (Dunavant 1985)

每个 RWG 边有两个三角形（T⁺ 和 T⁻），因此 `compute_interaction` 分解为 4 个子积分：
- `Z_ij = Z(T_i⁺,T_j⁺) + Z(T_i⁺,T_j⁻) + Z(T_i⁻,T_j⁺) + Z(T_i⁻,T_j⁻)`

### 1.3 奇异积分处理（`same_tri == true`）

当源三角形与场三角形为同一个三角形时，Green 函数 $G(R)$ 在 $R \to 0$ 处奇异。

**处理策略**：将 Green 函数拆分为平滑部分 + 静态奇异部分：

$$
G(k,R) = \frac{e^{-jkR} - 1}{4\pi R} + \frac{1}{4\pi R}
$$

| 部分 | 方法 | 说明 |
|------|------|------|
| **平滑部分** `G_smooth` | 7 点 Gauss 双重积分 | $\lim_{R\to0} G_{\text{smooth}} = -jk/(4\pi)$，解析无奇点 |
| **静态奇异部分** | 解析线积分 + 中心点近似 | 利用恒等式 $\nabla'^2_s R = 1/R$（平坦三角形）+ 表面散度定理 |

**静态奇异部分的具体计算**：

1. **内层积分**：$\int_T \frac{1}{|\mathbf{r}-\mathbf{r}'|} dS'$ 用解析线积分：
   $$
   \oint_{\partial T} \mathbf{P}_i^0 \cdot \mathbf{f}_i \, dl, \quad \mathbf{f}_i = \ln\frac{R_A + R_B + L}{R_A + R_B - L}
   $$
   其中 $\mathbf{P}_i^0 = \mathbf{n}_{\text{out}} \cdot (\mathbf{A} - \mathbf{r})$

2. **外层积分**：用 7 点 Gauss 对 `integrate_1_over_R` 的结果再积分，得到 $I_{\text{double}} = \iint \frac{1}{R} dS' dS$

3. **静态修正** = $jk\eta \cdot \frac{K(\mathbf{r}_c)}{4\pi} \cdot I_{\text{double}}$，其中 $K = \mathbf{f}_i \cdot \mathbf{f}_j - \frac{1}{k^2}(\nabla\cdot\mathbf{f}_i)(\nabla\cdot\mathbf{f}_j)$ 在三角形中心 $\mathbf{r}_c$ 取值

**关键函数**：
| 函数 | 位置 | 作用 |
|------|------|------|
| `green_kernel()` | `efie.cpp:41` | 自由空间 Green 函数 $e^{-jkR}/(4\pi R)$ |
| `green_smooth()` | `efie.cpp:49` | 平滑化 Green 函数 $(e^{-jkR}-1)/(4\pi R)$ |
| `integrate_1_over_R()` | `efie.cpp:61` | 解析计算 $\int_T 1/R \, dS'$ |
| `integrateTrianglePair()` | `efie.cpp:139` | 三角形对上的 RWG-EFIE 双重积分 |
| `compute_interaction()` | `efie.cpp:198` | 两条边的完整相互作用 $Z_{ij}$ |

### 1.4 激励向量填充 — `fill_vector()`

$$
b_i = \int_{T_i^+ \cup T_i^-} \mathbf{f}_i(\mathbf{r}) \cdot \mathbf{E}^{\text{inc}}(\mathbf{r}) \, dS
$$

入射场为平面波：$\mathbf{E}^{\text{inc}}(\mathbf{r}) = \mathbf{E}_0 \cdot e^{-j\mathbf{k}^{\text{inc}} \cdot \mathbf{r}}$

每个 RWG 边的两个三角形（T⁺, T⁻）分别用 7 点 Gauss 积分，内联计算点积 `E0_pol·f(r)` 和相位因子 `exp(-j k_inc·r)`。

### 1.5 无用函数

`plane_wave_e()` (`efie.cpp:313-320`) — 定义了但从未被调用，且返回类型不合理（矢量电场被压缩成标量 Complex），属于遗留死代码。

---

## 二、CG 求解器模块 (`cg_solver.h` / `cg_solver.cpp`)

### 2.1 模块职责

求解复线性系统 **Z·I = b**，提供两种迭代算法。

### 2.2 CG（共轭梯度法）— `solve()`

**适用范围**：仅 Hermitian 正定矩阵。EFIE 阻抗阵是复非对称的，因此**本项目中不适用**，仅为教学对比保留。

**算法流程**：

```
x₀ = 0, r₀ = b, p₀ = r₀
for iter = 0, 1, ...:
    α = (rᵀr) / (pᵀZp)
    x ← x + α·p
    r ← r - α·Zp
    if ‖r‖² < tol²·‖b‖²: 收敛退出
    β = ‖r_new‖² / ‖r_old‖²
    p ← r + β·p
```

**收敛条件**：`‖r‖² < tol² × ‖b‖²`（默认 tol=1e-6, max_iter=1000）

### 2.3 BiCGSTAB（双共轭梯度稳定法）— `solve_bicgstab()`

**适用范围**：复非对称矩阵（EFIE 阻抗阵的标准选择）。基于 van der Vorst (1992)。

**算法流程**：

```
x₀ = 0, r₀ = b, r̂₀ = r₀
ρ_old = α = ω = 1
v₀ = p₀ = 0

for iter = 0, 1, ...:
    ρ = ⟨r̂, r⟩
    β = (ρ/ρ_old) × (α/ω)    （iter>0 时）
    p = r + β(p - ω·v)
    v = Z·p
    α = ρ / ⟨r̂, v⟩
    s = r - α·v
    t = Z·s
    ω = ⟨t, s⟩ / ⟨t, t⟩
    x ← x + α·p + ω·s
    r ← s - ω·t
    if ‖r‖ < tol × max(‖b‖, 1e-30): 收敛退出
```

**与 CG 的关键区别**：
- CG 每次迭代 1 次矩阵-向量乘；BiCGSTAB 需要 2 次
- BiCGSTAB 不需要矩阵对称/Hermitian
- BiCGSTAB 有更多的 breakdown 检测（ρ≈0, r̂·v≈0, t·t≈0）

**收敛条件**：`‖r‖ < tol × max(‖b‖, 1e-30)`（默认 tol=1e-6, max_iter=1000）

---

## 三、远场 / RCS 模块 (`farfield.h` / `farfield.cpp`)

### 3.1 模块职责

由 MoM 求得的电流系数 **I** 计算远区散射电场和单站 RCS。

### 3.2 辐射矢量 — `compute_N()`

$$
\mathbf{N}(\theta, \phi) = \sum_{n=1}^{N} I_n \int_{S_n} \mathbf{f}_n(\mathbf{r}') \, e^{j k \, \hat{\mathbf{k}}_s \cdot \mathbf{r}'} \, dS'
$$

- $\hat{\mathbf{k}}_s = (\sin\theta\cos\phi, \sin\theta\sin\phi, \cos\theta)$ — 散射方向
- 对每个 RWG 边的两个支撑三角形分别用 7 点 Gauss 积分
- 被积函数为三维矢量：$\mathbf{f}_n$ 的三个分量 × $e^{j\phi}$

### 3.3 远区散射电场

由辐射矢量投影到球坐标基矢量：

$$
E_\theta = -\frac{jk\eta}{4\pi} \, N_\theta, \quad N_\theta = N_x \cos\theta\cos\phi + N_y \cos\theta\sin\phi - N_z \sin\theta
$$

$$
E_\phi = -\frac{jk\eta}{4\pi} \, N_\phi, \quad N_\phi = -N_x \sin\phi + N_y \cos\phi
$$

远场因子 $e^{-jkr}/r$ 已省略（单位远场），$\eta = Z_0 \approx 376.73\,\Omega$。

### 3.4 RCS 计算

`compute_rcs()` 计算给定观测方向上的 RCS，公式本身对单站和双站通用：

$$
\sigma(\theta, \phi) = 4\pi \frac{|E_\theta(\theta, \phi)|^2 + |E_\phi(\theta, \phi)|^2}{|E_{\text{inc}}|^2}
$$

其中 $|E_{\text{inc}}| = 1$ V/m（单位入射场），RCS 以 dBsm 表示：

$$
\text{RCS}_{\text{dBsm}} = 10 \log_{10}(\sigma)
$$

> **注意**：当前 `main.cpp` 中入射方向固定为 (θ=45°, φ=0°)，而观测方向 θ 从 0° 变到 180°，入射与观测方向不同——因此实际计算的是**双站 RCS**（bistatic）。函数名为 `compute_rcs_monostatic` 是命名上的误导。要获得真正的单站 RCS，需在每个观测角度将入射方向设为同方向后重新求解。

**关键函数**：
| 函数 | 作用 |
|------|------|
| `compute_N()` | 计算辐射矢量 $\mathbf{N}$（内部辅助） |
| `compute_Etheta()` | 远场 $E_\theta$ 分量 |
| `compute_Ephi()` | 远场 $E_\phi$ 分量 |
| `compute_rcs()` | RCS (dBsm)，入射方向由求解时的 I 决定 |

---

## 四、Main 函数 (`main.cpp`)

### 4.1 整体计算流水线

```
参数设置 → 网格生成 → RWG构建 → EFIE填充 → 线性求解 → RCS计算
```

### 4.2 逐步流程

| 步骤 | 操作 | 关键参数/输出 |
|------|------|---------------|
| **1. 参数** | 硬编码物理参数 | `radius=0.5m`, `freq=300MHz`, `ka≈3.14`, `refinement=2` (≈320 三角) |
| **2. 网格** | `mesh.generate_sphere()` | 正二十面体 → 递归细分 → 投影到球面 |
| **3. RWG** | `rwg.build(mesh)` | 约 480 个未知量 (≈1.5×三角形数) |
| **4. EFIE** | `efie.fill_matrix(Z)` | 480×480 复稠密矩阵，O(N²) 填充，7-pt Gauss |
| **5. 激励** | `efie.fill_vector(k_inc, E0, b)` | 入射方向 θ=45°, φ=0°, φ-极化 (E₀ = ŷ) |
| **6. 求解** | `CGSolver::solve_bicgstab(Z,b,I)` | BiCGSTAB, tol=1e-6, max_iter=1000 |
| **7. RCS** | `farfield.compute_rcs()` | θ = 0°~180°, 固定入射方向 → **实际为双站 RCS** |

### 4.3 硬编码参数一览

| 参数 | 值 | 说明 |
|------|-----|------|
| `radius` | 0.5 m | 导体球半径 |
| `freq` | 300 MHz | 入射波频率 |
| `lambda` | 1.0 m | 波长 (c₀/f) |
| `k0` | 2π | 波数 |
| `refinement` | 2 | 网格细分次数 |
| `inc_theta` | 45° | 入射俯仰角 |
| `inc_phi` | 0° | 入射方位角 |
| `E0_pol` | (0, 1, 0) | φ 极化（垂直极化） |
| `n_theta` | 5 | RCS θ 采样点数 |
| `solver tol` | 1e-6 | BiCGSTAB 收敛容差 |
| `solver max_iter` | 1000 | BiCGSTAB 最大迭代次数 |

### 4.4 输出

- 频率、波长、ka 值
- 阻抗矩阵填充进度（每 50 行）
- BiCGSTAB 收敛信息（迭代次数，相对残差）
- RCS 表格：θ（度）vs RCS（dBsm）

---

## 五、关键设计决策与限制

1. **双重 Gauss 积分**：选择 7 点规则（5 阶精度）在精度与计算量之间折中
2. **奇异积分**：仅处理同三角形自作用（`same_tri`），相邻三角形（共边/共顶点）不做特殊处理
3. **矩阵存储**：`vector<vector<Complex>>` 满存储，无压缩/稀疏格式，O(N²) 内存
4. **EFIE 直接填充**：无 FMM/AIM 等加速算法，O(N²) 矩阵填充 + O(N²) 每迭代
5. **极化假设**：入射场为 $\hat{y}$ 方向线极化（φ 极化），对所有散射角适用
6. **教学导向**：所有参数硬编码，侧重于算法清晰性而非数值精度或计算效率
