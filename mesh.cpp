#include "mesh.h"
#include <cmath>
#include <unordered_map>
#include <iostream>

//面积
double Triangle::area(const std::vector<Vertex>& verts) const {
    const Vertex& a = verts[v0];
    const Vertex& b = verts[v1];
    const Vertex& c = verts[v2];
    //三角形面积 = 两条边向量的叉积长度的一半
    double ax = b.x - a.x, ay = b.y - a.y, az = b.z - a.z;
    double bx = c.x - a.x, by = c.y - a.y, bz = c.z - a.z;
    double cx = ay * bz - az * by;
    double cy = az * bx - ax * bz;
    double cz = ax * by - ay * bx;
    return 0.5 * sqrt(cx*cx + cy*cy + cz*cz);
}

//中心（重心）三个顶点坐标的平均值，将三角形分成面积相等的三个小三角形
Vertex Triangle::center(const std::vector<Vertex>& verts) const {
    const Vertex& a = verts[v0];
    const Vertex& b = verts[v1];
    const Vertex& c = verts[v2];
    return Vertex((a.x + b.x + c.x)/3.0,
                  (a.y + b.y + c.y)/3.0,
                  (a.z + b.z + c.z)/3.0);
}

//三角形平面的法向量
Vertex Triangle::normal(const std::vector<Vertex>& verts) const {
    const Vertex& a = verts[v0];
    const Vertex& b = verts[v1];
    const Vertex& c = verts[v2];
    //两条边的叉积(A->B->C AB×AC方向垂直纸面向外)
    double ax = b.x - a.x, ay = b.y - a.y, az = b.z - a.z;
    double bx = c.x - a.x, by = c.y - a.y, bz = c.z - a.z;
    double nx = ay * bz - az * by;
    double ny = az * bx - ax * bz;
    double nz = ax * by - ay * bx;
    //模
    double len = sqrt(nx*nx + ny*ny + nz*nz);
    //当三角形的几乎为0，避免除于0
    if (len < 1e-12) return Vertex(0,0,0);
    return Vertex(nx/len, ny/len, nz/len);
}

// 生成正二十面体的12个顶点
static void generate_icosahedron_vertices(std::vector<Vertex>& verts) {
    double phi = (1.0 + sqrt(5.0)) / 2.0; // 黄金比例
    double a = 1.0;
    double b = 1.0 / phi;
    verts.clear();
    verts.push_back(Vertex( a,  b,  0));  // 0
    verts.push_back(Vertex( a, -b,  0));  // 1
    verts.push_back(Vertex(-a,  b,  0));  // 2
    verts.push_back(Vertex(-a, -b,  0));  // 3
    verts.push_back(Vertex( b,  0,  a));  // 4
    verts.push_back(Vertex( b,  0, -a));  // 5
    verts.push_back(Vertex(-b,  0,  a));  // 6
    verts.push_back(Vertex(-b,  0, -a));  // 7
    verts.push_back(Vertex( 0,  a,  b));  // 8
    verts.push_back(Vertex( 0,  a, -b));  // 9
    verts.push_back(Vertex( 0, -a,  b));  // 10
    verts.push_back(Vertex( 0, -a, -b));  // 11
}

// 生成正二十面体的20个三角形面
static void generate_icosahedron_triangles(std::vector<Triangle>& tris) {
    // 预定义的20个三角形（顶点索引）
    int idx[20][3] = {
        {0,1,4}, {0,4,8}, {0,8,9}, {0,9,5}, {0,5,1},
        {1,5,11}, {1,11,10}, {1,10,4}, {2,3,6}, {2,6,8},
        {2,8,9}, {2,9,7}, {2,7,3}, {3,7,11}, {3,11,10},
        {3,10,6}, {4,10,6}, {4,6,8}, {5,9,7}, {5,7,11}
    };
    for (int i = 0; i < 20; i++) {
        tris.push_back(Triangle(idx[i][0], idx[i][1], idx[i][2], i));
    }
}

void Mesh::generate_sphere(double radius, int refinement_level) {
    // 1. 生成正二十面体
    generate_icosahedron_vertices(vertices);
    generate_icosahedron_triangles(triangles);
    n_vertices = vertices.size();
    n_triangles = triangles.size();
    
    // 2. 递归细分
    for (int level = 0; level < refinement_level; level++) {
        subdivide(level, refinement_level);
    }
    
    // 3. 归一化到单位球面并缩放
    normalize_vertices(radius);
    
    std::cout << "Grid generated: " << vertices.size() 
              << " vertices, " << triangles.size() << " triangles\n";
}

void Mesh::subdivide(int level, int max_level) {
    // 对每个三角形，在边的中点插入新顶点，将一个三角形分成4个
    std::vector<Triangle> new_triangles;
    //记录每条边对应的中点索引
    std::unordered_map<long long, int> edge_midpoint;
    
    //为一条边生成唯一的ID
    auto get_key = [](int v0, int v1) -> long long {
        if (v0 > v1) std::swap(v0, v1);
        //高32位存v0，低32位存v1
        return (long long)v0 << 32 | v1;
    };
    
    for (const auto& tri : triangles) {
        int mid12 = -1, mid23 = -1, mid31 = -1;
        
        auto get_mid = [&](int v0, int v1) -> int {
            // 检查这条边是否已经创建过中点
            long long key = get_key(v0, v1);
            auto it = edge_midpoint.find(key);
            if (it != edge_midpoint.end()) return it->second;
            Vertex p0 = vertices[v0];
            Vertex p1 = vertices[v1];
            Vertex mid((p0.x+p1.x)/2, (p0.y+p1.y)/2, (p0.z+p1.z)/2);
            int idx = vertices.size();
            vertices.push_back(mid);
            //这条边的中点在新顶点列表中的索引
            edge_midpoint[key] = idx;
            return idx;
        };
        
        mid12 = get_mid(tri.v0, tri.v1);
        mid23 = get_mid(tri.v1, tri.v2);
        mid31 = get_mid(tri.v2, tri.v0);
        
        new_triangles.push_back(Triangle(tri.v0, mid12, mid31, 0));
        new_triangles.push_back(Triangle(tri.v1, mid23, mid12, 0));
        new_triangles.push_back(Triangle(tri.v2, mid31, mid23, 0));
        new_triangles.push_back(Triangle(mid12, mid23, mid31, 0));
    }
    
    // 更新三角形索引
    for (size_t i = 0; i < new_triangles.size(); i++) {
        new_triangles[i].idx = i;
    }
    triangles = std::move(new_triangles);
    n_triangles = triangles.size();
}

void Mesh::normalize_vertices(double radius) {
    for (auto& v : vertices) {
        double len = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
        if (len > 1e-12) {
            v.x = v.x / len * radius;
            v.y = v.y / len * radius;
            v.z = v.z / len * radius;
        }
    }
}