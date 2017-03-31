//
// Created by 汪子琦 on 3/23/17.
//

#ifndef SUPPORTER_SLICE_H
#define SUPPORTER_SLICE_H
#include <Eigen/Dense>
#include <igl/triangle/triangulate.h>
#include <vector>
#include <list>
#include <cmath>
#include "clipper.hpp"
#include "settings.h"
#include "hash_edge.h"
class EdgeHash;
class Slice
{
public:
    Slice()
    {
        V.setZero();
        F.setZero();
    }

    ~Slice()
    {

    }

public:

    void set_mesh(Eigen::MatrixXd &in_V, Eigen::MatrixXi &in_F);
    void set_layer_height();
    bool check_polygon(ClipperLib::Path poly);
    inline int number_layer() { return (max_Y() - min_Y()) / settings.mm2int_Even(settings.layer_height) + 1; }
    inline double max_Y() { assert(!V.isZero());return V.colwise().maxCoeff()[1]; }
    inline double min_Y() { assert(!V.isZero());return V.colwise().minCoeff()[1]; }

protected:
    void build_triangle_list(std::vector<std::vector<int>> &L);
    void compute_intersection(int fd, int Y, Eigen::MatrixXi &mat);
    void incremental_slicing(std::vector< std::vector<Eigen::MatrixXi>> &S);

public:
    void contour_construction();
    void get_intersecting_surface(int layer, Eigen::MatrixXd &V2, Eigen::MatrixXi &F2);

protected:
    Eigen::MatrixXi V;
    Eigen::MatrixXi F;
    std::vector<int> P;
    std::vector< ClipperLib::Paths> layer_slices;
    Settings settings;
};

void Slice::set_mesh(Eigen::MatrixXd &in_V, Eigen::MatrixXi &in_F)
{
    F = in_F;
    V.resize(in_V.rows(), in_V.cols());
    for(size_t id = 0; id < V.rows(); id++)
        for(size_t jd = 0; jd < V.cols(); jd++)
            V(id, jd) = settings.mm2int_Even(in_V(id, jd));

    return;
}

void Slice::set_layer_height()
{
    assert(!V.isZero());

    int dy = settings.mm2int_Even(settings.layer_height);
    int sizeL = (max_Y() - min_Y()) / dy + 1;

    P.clear();
    P.resize(sizeL);
    P[0] = min_Y() + 1;
    for(size_t id = 1;id < sizeL; id++)
        P[id] = P[id - 1] + dy;

    return;
}

void Slice::build_triangle_list(std::vector<std::vector<int>> &L)
{
    set_layer_height();
    int sizeL = number_layer();
    int dy = settings.mm2int_Even(settings.layer_height);

    L.clear();
    L.resize(sizeL + 1);

    for(int id = 0; id < F.rows(); id++)
    {
        Eigen::RowVector3i zvalue(V(F(id, 0),1), V(F(id, 1),1), V(F(id, 2),1));
        int zmin = zvalue.minCoeff();
        if(zmin < P[0])
            L[0].push_back(id);
        else if(zmin > P[sizeL - 1])
            L[sizeL].push_back(id);
        else L[(zmin - P[0]) / dy + 1].push_back(id);
    }

    return;
}

bool Slice::check_polygon(ClipperLib::Path poly)
{
    int size = poly.size();
    for(size_t id = 0; id < poly.size(); id++)
    {
        if(poly[id] == poly[(id + 1) % poly.size()])
            size--;
    }

    if(size >= 3) return true;
    else return false;
}

void Slice::compute_intersection(int fd, int Y, Eigen::MatrixXi &mat)
{

    Eigen::RowVector3d p[3], n, dy, dq,q[2];
    for(int id = 0; id < 3; id++)
        p[id] = Eigen::RowVector3d(V.row(F(fd, id))[0],
                                   V.row(F(fd, id))[1],
                                   V.row(F(fd, id))[2]);

    int jd = 0;
    for(int id = 0; id < 3; id++)
    {
        //id and id + 1
        Eigen::RowVector3d u(p[id]),
                v(p[(id + 1) % 3]);
        if((u.y() - Y) * (v.y() - Y) < 0)
        {
            q[jd ++] = (double)(u.y() - Y) / (u.y() - v.y()) * (v - u) + u;
            //std::cout << q[jd - 1] << std::endl;
        }
    }

    //std::cout << q[0](0) << ", " << q[0](1) << ", " << q[0](2) << std::endl
    //          << q[1](0) << ", " << q[1](1) << ", " << q[1](2) << std::endl;
    n = (p[1] - p[0]).cross(p[2] - p[0]);
    n /= n.norm();
    dy = Eigen::RowVector3d(0, 1, 0);
    dq = n.cross(dy);
    mat.resize(2, 2);
    if((q[1] - q[0]).dot(dq) > 0)
    {
        mat <<  q[0][0], q[0][2],
                q[1][0], q[1][2];
    }
    else
    {
        mat <<  q[1][0], q[1][2],
                q[0][0], q[0][2];
    }

    //std::cout << mat << std::endl;
    return;
}

void Slice::incremental_slicing(std::vector< std::vector<Eigen::MatrixXi>> &S)
{
    assert(!V.isZero() && layer_slices.empty());

    std::vector<std::vector<int>> L;
    build_triangle_list(L);
    int sizeL = number_layer();
    S.resize(sizeL);
//        for(int id = 0; id < P.size(); id++)
//        std::cout << P[id] << std::endl;

    std::list<int> A;
    for(int id = 0; id < sizeL; id++)
    {
        for(int jd = 0; jd < L[id].size(); jd++)
            A.push_back(L[id][jd]);
        std::list<int>::iterator it = A.begin();
        while(it != A.end())
        {
            Eigen::RowVector3i zvalue(V(F(*it, 0),1), V(F(*it, 1),1), V(F(*it, 2),1));
            int zmax = zvalue.maxCoeff();
            if(zmax < P[id])
                it = A.erase(it);
            else
            {
                Eigen::MatrixXi mat;
                compute_intersection(*it, P[id], mat);
                if(mat.row(0) != mat.row(1))
                    S[id].push_back(mat);
                //std::cout << mat << std::endl;
                it++;
            }
        }
    }
}

void Slice::contour_construction()
{
    assert(!V.isZero() && layer_slices.empty());
    std::vector< std::vector<Eigen::MatrixXi>> S;
    std::vector<bool> visited;

    incremental_slicing(S);
    layer_slices.resize(S.size());
    for(int id = 0; id < S.size(); id++)
    {
        std::cout << "slice :" << id << "/" << S.size() << " layer" << std::endl;
        ClipperLib::Paths paths;
        visited.clear();
        visited.resize(S[id].size(), false);
        HashEdge hash;
        for(int jd = 0; jd < S[id].size(); jd++)
        {
            Eigen::RowVector2i u = S[id][jd].row(0);
            Eigen::RowVector2i v = S[id][jd].row(1);
            hash.build_connection(u, v, jd);
        }

        for (int jd = 0; jd < S[id].size(); jd++)
        {
            if(!visited[jd])
            {
                Eigen::RowVector2i last = S[id][jd].row(0);
                Eigen::RowVector2i u = S[id][jd].row(0), v;
                ClipperLib::Path path;
                bool flag = false;
                do
                {
                    int idx;
                    flag = hash.find(u, v, idx);
                    if(!flag) break;
                    //std::cout << u << "," << v << std::endl;
                    path.push_back(ClipperLib::IntPoint(v(0), v(1)));
                    u = v;
                    visited[idx] = true;
                }while(u != last);
                if(flag) paths.push_back(path);
            }
        }
        layer_slices[id] = paths;
    }

}

void Slice::get_intersecting_surface(int layer, Eigen::MatrixXd &V2, Eigen::MatrixXi &F2) {

    assert(!V.isZero());
    if(layer_slices.empty())
        contour_construction();
    assert(0 <= layer && layer < layer_slices.size());

    ClipperLib::Paths solution;
    ClipperLib::Clipper clipper;
    clipper.AddPaths(layer_slices[layer], ClipperLib::ptClip, true);
    clipper.Execute(ClipperLib::ctUnion, solution, ClipperLib::pftPositive, ClipperLib::pftPositive);

    std::vector<Eigen::RowVector2i> holes;
    std::vector<Eigen::RowVector2i> vertices;
    std::vector<Eigen::RowVector2i> edges;

    int size = 0;

    for (size_t id = 0; id < solution.size(); id++) {
        Eigen::RowVector2d center(0, 0), hole(0, 0);
        for (size_t jd = 0; jd < solution[id].size(); jd++) {
            Eigen::RowVector2i pt(solution[id][jd].X, solution[id][jd].Y);
            Eigen::RowVector2i edge(size + jd, size + (jd + 1) % solution[id].size());
            vertices.push_back(pt);
            edges.push_back(edge);
            hole[0] = pt[0];
            hole[1] = pt[1];
            center += hole;
        }

        if (!ClipperLib::Orientation(solution[id])) {
            center /= solution[id].size();
            hole = (hole - center) * 0.999f + center;
            holes.push_back(Eigen::RowVector2i(hole[0], hole[1]));
        }
        size += solution[id].size();
    }

    Eigen::MatrixXd V1(vertices.size(), 2);
    Eigen::MatrixXi E(edges.size(), 2);
    Eigen::MatrixXi H(holes.size(), 2);
    for (size_t id = 0; id < vertices.size(); id++)
    {
        V1(id, 0) = vertices[id][0];
        V1(id, 1) = vertices[id][1];
    }
    for (size_t id = 0; id < edges.size(); id++)
        E.row(id) = edges[id];
    for(size_t id = 0; id < holes.size(); id++)
        H.row(id) = holes[id];

    Eigen::MatrixXi in_V2;
    if(E.rows() >= 3 && V1.rows() >= 3)
    {
        igl::triangle::triangulate(V1, E, H, "Q", in_V2, F2);
        V2.resize(in_V2.rows(), 3);
        for (size_t id = 0; id < in_V2.rows(); id++)
        {
            V2(id, 0) = settings.int2mm(in_V2(id, 0));
            V2(id, 1) = settings.int2mm(P[layer]);
            V2(id, 2) = settings.int2mm(in_V2(id, 1));
        }
    }
    else
    {
        V2 = Eigen::MatrixXd();
        F2 = Eigen::MatrixXi();
    }
    return;
}

#endif //SUPPORTER_SLICE_H