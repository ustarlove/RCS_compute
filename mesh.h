#ifndef MESH_H
#define MESH_H

#include <vector>
#include "constants.h"

struct Vertex {
    double x, y, z;
    Vertex() : x(0), y(0), z(0) {}
    Vertex(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
};

struct Triangle {
    int v0, v1, v2;   // 顶点索引
    int idx;           // 三角形编号
    Triangle(int v0_, int v1_, int v2_, int idx_) 
        : v0(v0_), v1(v1_), v2(v2_), idx(idx_) {}
    
    // 计算三角形面积
    double area(const std::vector<Vertex>& verts) const;
    
    // 计算三角形中心
    Vertex center(const std::vector<Vertex>& verts) const;
    
    // 计算三角形法向（单位）
    Vertex normal(const std::vector<Vertex>& verts) const;
};

struct Mesh {
    std::vector<Vertex> vertices;
    std::vector<Triangle> triangles;
    int n_vertices;
    int n_triangles;
    
    // 生成单位球面网格（通过递归细分正二十面体）
    void generate_sphere(double radius, int refinement_level);
    
private:
    void subdivide(int level, int max_level);
    //void add_triangle(int v0, int v1, int v2);
    void normalize_vertices(double radius);
};

#endif // MESH_H