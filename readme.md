sphere\_rcs/
│
├── main.cpp                 # 主程序入口
├── mesh.h                   # 网格数据结构声明
├── mesh.cpp                 # 球面三角剖分实现
├── rwg.h                    # RWG基函数声明
├── rwg.cpp                  # RWG基函数实现
├── efie.h                   # EFIE矩阵和向量声明
├── efie.cpp                 # EFIE矩阵填充实现
├── cg\_solver.h              # 共轭梯度法求解器声明
├── cg\_solver.cpp            # 共轭梯度法实现
├── farfield.h               # 远场计算声明
├── farfield.cpp             # 远场计算实现
├── constants.h              # 常量定义（π, ε₀, μ₀等）
├── Makefile                 # 编译脚本
└── README.md                # 说明文件



项目中几乎所有涉及电磁计算的量都是复数的——Green 函数、相位因子、阻抗矩阵、电流系数、辐射矢量、散射电场等，因为频域 Maxwell

&#x20; 方程组天然产生复数运算。只有几何量（坐标、边长、面积）和物理常数保持实数。

