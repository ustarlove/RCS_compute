#include "rwg.h"
#include <unordered_map>
#include <iostream>

void RWGData::build(const Mesh& mesh) {
    // 建立边到三角形的映射
    std::unordered_map<long long, std::pair<int, int>> edge_tri_map;

    auto make_key = [](int v0, int v1) -> long long {
        if (v0 > v1) std::swap(v0, v1);
        return (long long)v0 << 32 | v1;
    };

    for (const auto& tri : mesh.triangles) {
        int verts[3] = {tri.v0, tri.v1, tri.v2};
        for (int i = 0; i < 3; i++) {
            int v0 = verts[i];
            int v1 = verts[(i+1)%3];
            long long key = make_key(v0, v1);
            auto it = edge_tri_map.find(key);
            if (it == edge_tri_map.end()) {
                edge_tri_map[key] = std::make_pair(tri.idx, -1);
            } else {
                it->second.second = tri.idx;
            }
        }
    }

    // 收集公共边（两个相邻三角形都存在）
    edges.clear();
    edge_to_t1.clear();
    edge_to_t2.clear();
    rwg_coeff1.clear();
    rwg_coeff2.clear();

    for (const auto& entry : edge_tri_map) {
        int t1 = entry.second.first;
        int t2 = entry.second.second;
        if (t2 == -1) continue;  // 跳过边界边（封闭球体没有边界）

        long long key = entry.first;
        int v0 = (int)(key >> 32);
        int v1 = (int)(key & 0xFFFFFFFF);

        const Vertex& p0 = mesh.vertices[v0];
        const Vertex& p1 = mesh.vertices[v1];
        double len = sqrt(pow(p1.x-p0.x,2)+pow(p1.y-p0.y,2)+pow(p1.z-p0.z,2));

        double A1 = mesh.triangles[t1].area(mesh.vertices);
        double A2 = mesh.triangles[t2].area(mesh.vertices);

        // 找到每个三角形中与边相对的顶点 (apex)
        auto find_apex = [&](int tri_idx, int ev0, int ev1) -> Vertex {
            const Triangle& tri = mesh.triangles[tri_idx];
            for (int vi : {tri.v0, tri.v1, tri.v2}) {
                if (vi != ev0 && vi != ev1) {
                    return mesh.vertices[vi];
                }
            }
            return Vertex(0,0,0);
        };

        int edge_idx = edges.size();
        Edge e;
        e.i = edge_idx;
        e.t1 = t1;
        e.t2 = t2;
        e.v1 = v0;
        e.v2 = v1;
        e.length = len;
        e.apex1 = find_apex(t1, v0, v1);
        e.apex2 = find_apex(t2, v0, v1);
        edges.push_back(e);
        edge_to_t1.push_back(t1);
        edge_to_t2.push_back(t2);
        rwg_coeff1.push_back(len / (2.0 * A1));   // l / (2A⁺)
        rwg_coeff2.push_back(len / (2.0 * A2));   // l / (2A⁻)
    }

    std::cout << "RWG edges: " << edges.size() << std::endl;
}

void RWGData::get_rwg(int edge_idx, int tri_idx,
                      const Vertex& r, Complex* value) const {
    const Edge& e = edges[edge_idx];

    double coeff = 0.0;
    Vertex apex;
    if (tri_idx == e.t1) {
        coeff = rwg_coeff1[edge_idx];          // +l/(2A⁺)
        apex = e.apex1;
    } else if (tri_idx == e.t2) {
        coeff = -rwg_coeff2[edge_idx];         // -l/(2A⁻)
        apex = e.apex2;
    } else {
        value[0] = Complex(0, 0);
        value[1] = Complex(0, 0);
        value[2] = Complex(0, 0);
        return;
    }

    // f(r) = coeff * (r - apex)
    value[0] = Complex(coeff * (r.x - apex.x), 0);
    value[1] = Complex(coeff * (r.y - apex.y), 0);
    value[2] = Complex(coeff * (r.z - apex.z), 0);
}

double RWGData::get_divergence(int edge_idx, int tri_idx) const {
    // ∇·f⁺ = l/A⁺ = 2 * rwg_coeff1  
    // ∇·f⁻ = -l/A⁻ = -2 * rwg_coeff2
    if (tri_idx == edges[edge_idx].t1) {
        return 2.0 * rwg_coeff1[edge_idx];
    } else if (tri_idx == edges[edge_idx].t2) {
        return -2.0 * rwg_coeff2[edge_idx];
    }
    return 0.0;
}
