#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <complex>
#include <cmath>

const double PI = 3.14159265358979323846;
const double EPS0 = 8.854187817e-12;   // 真空介电常数 F/m
const double MU0 = 4.0 * PI * 1e-7;    // 真空磁导率 H/m
const double Z0 = sqrt(MU0 / EPS0);    // 真空波阻抗 ≈ 376.73 Ω

using Complex = std::complex<double>;

//度数转为弧度
inline double deg2rad(double deg) {
    return deg * PI / 180.0;
}

#endif // CONSTANTS_H