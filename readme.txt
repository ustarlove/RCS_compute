改动总结

  1. rwg.h / rwg.cpp — RWG 基函数完整实现

  - Edge 结构添加了 apex1、apex2（存储三角形中与边相对的顶点坐标），这是 RWG 矢量基函数 f(r) = ±l/(2A) * (r - apex) 的必要几何信息
  - get_rwg() 从空壳改为完整的矢量计算
  - 新增 get_divergence() 计算散度 ∇·f = ±l/A（EFIE 标量势项需要）

  2. efie.cpp — EFIE 阻抗矩阵完整实现（核心改动）

  - 7 点 Gauss 积分（Dunavant 1985，5 阶精度）替代原来的单中心点近似
  - 双重面积分：对每对三角形 (T_i^p, T_j^q) 使用 7×7=49 点 Gauss 积分，正确计算 Galerkin 测试
  - 完整 EFIE 核：包含矢量势项 jkη f_i·f_j 和标量势项 (η/jk)(∇·f_i)(∇·f_j)
  - 奇异积分处理：提取静态核 1/(4πR) 单独处理，平滑部分用 Gauss 积分
  - 激励向量：对每个 RWG 边做 7 点 Gauss 积分计算 b_i = ∫ f_i·E_inc dS

  3. cg_solver.h / cg_solver.cpp — 新增 BiCGSTAB 求解器

  - 原始 CG 仅适用于 Hermitian 正定矩阵，EFIE 阻抗阵为复非对称
  - 实现标准 BiCGSTAB 算法，正确处理非对称复矩阵

  4. farfield.cpp — 远场积分完善

  - 辐射积分 ∫ f_n(r) exp(jk k̂·r) dS 使用 7 点 Gauss 积分替代中心点近似
  - 正确转换为球坐标分量 E_θ、E_φ

  5. main.cpp — 切换到 BiCGSTAB

  - CGSolver::solve → CGSolver::solve_bicgstab，最大迭代从 500 增加到 1000