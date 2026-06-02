#ifndef VEC3_OPS_H
#define VEC3_OPS_H

#include "mesh.h"
#include <cmath>

inline Vertex operator+(const Vertex& a, const Vertex& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

inline Vertex operator-(const Vertex& a, const Vertex& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

inline Vertex operator*(double s, const Vertex& v) {
    return {s * v.x, s * v.y, s * v.z};
}

inline Vertex operator*(const Vertex& v, double s) {
    return {s * v.x, s * v.y, s * v.z};
}

inline Vertex operator/(const Vertex& v, double s) {
    return {v.x / s, v.y / s, v.z / s};
}

inline double dot(const Vertex& a, const Vertex& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline Vertex cross(const Vertex& a, const Vertex& b) {
    return {a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x};
}

inline double norm(const Vertex& v) {
    return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

inline Vertex normalize(const Vertex& v) {
    double n = norm(v);
    if (n < 1e-15) return {0, 0, 0};
    return {v.x / n, v.y / n, v.z / n};
}

#endif
