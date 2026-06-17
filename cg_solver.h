#ifndef CG_SOLVER_H
#define CG_SOLVER_H

#include <vector>
#include <functional>
#include "constants.h"

class CGSolver {
public:
    // 求解 Z x = b 用共轭梯度法（Z 须为 Hermitian 正定）
    static bool solve(const std::vector<std::vector<Complex>>& Z,
                      const std::vector<Complex>& b,
                      std::vector<Complex>& x,
                      double tol = 1e-6,
                      int max_iter = 1000);

    // 求解 Z x = b 用 BiCGSTAB（适用于复非对称矩阵，如 EFIE 阻抗阵）
    static bool solve_bicgstab(const std::vector<std::vector<Complex>>& Z,
                               const std::vector<Complex>& b,
                               std::vector<Complex>& x,
                               double tol = 1e-6,
                               int max_iter = 1000);

    // FMM 加速版 BiCGSTAB: 用 matvec 回调替代稠密 Z 矩阵
    static bool solve_bicgstab_functor(
        const std::function<void(const std::vector<Complex>& x,
                                  std::vector<Complex>& y)>& matvec,
        const std::vector<Complex>& b,
        std::vector<Complex>& x,
        double tol = 1e-6,
        int max_iter = 1000);

    // FMM 加速版 GMRES(m): restarted GMRES，残差单调下降
    static bool solve_gmres_functor(
        const std::function<void(const std::vector<Complex>& x,
                                  std::vector<Complex>& y)>& matvec,
        const std::vector<Complex>& b,
        std::vector<Complex>& x,
        double tol = 1e-6,
        int max_iter = 1000,
        int restart = 30);
};

#endif // CG_SOLVER_H