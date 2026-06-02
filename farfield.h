#ifndef FARFIELD_H
#define FARFIELD_H

#include "mesh.h"
#include "rwg.h"
#include "constants.h"

class FarField {
public:
    FarField(const Mesh& mesh, const RWGData& rwg, double k0);
    
    // 计算远区散射电场（单位远场，省略 e^{-jkr}/r 因子）
    Complex compute_Etheta(double theta, double phi,
                           const std::vector<Complex>& I) const;
    Complex compute_Ephi(double theta, double phi,
                         const std::vector<Complex>& I) const;
    
    // 计算单站 RCS（单位 dBsm）
    double compute_rcs(double theta, double phi,
                                  const std::vector<Complex>& I) const;
    
private:
    const Mesh& mesh_;
    const RWGData& rwg_;
    double k0_;
};

#endif // FARFIELD_H