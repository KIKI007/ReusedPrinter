//
// Created by 汪子琦 on 3/1/17.
//

#ifndef SUPPORTER_RENDERING_TREE_SUPPORT_H
#define SUPPORTER_RENDERING_TREE_SUPPORT_H

#include <vector>
#include <Eigen/Core>
#include <Eigen/StdVector>
#include <assert.h>

class RenderingTreeSupport {
public:

    RenderingTreeSupport()
    {

    }

public:

    unsigned long add_root(Eigen::Vector3d root_xyz);

    unsigned long add_node(Eigen::Vector3d node_xyz);

    void add_edge(unsigned long u0, unsigned u1);

    unsigned long size()
    {
        return node_coordinates_.size();
    }

public:

    void draw(Eigen::MatrixXd &V, Eigen::MatrixXi &F);

    double calculate_tree_length();

private:

    void generate_node_radius();

    void generate_trunk(unsigned long u, unsigned long v, std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i>> &vF);

    void generate_tip(unsigned long u, unsigned long v, std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i>> &vF);

    void generate_base(unsigned long u, unsigned long v, std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i>> &vF);

private:

    bool is_node_tip(unsigned long u)
    {
        assert(u < graph_.size());
        return graph_[u].empty();
    }

    bool is_node_root(unsigned long u)
    {
        assert(u < is_roots_.size());
        return is_roots_[u];
    }
private:

    std::vector< std::vector<unsigned long>> graph_;

    std::vector<bool> is_roots_;

private:

    std::vector< Eigen::Vector3d, Eigen::aligned_allocator< Eigen::Vector3d>> node_coordinates_;

    std::vector<double> node_radius_;

    std::vector< std::vector<unsigned long>> node_index_;

private:
    Settings settings_;
};

unsigned long RenderingTreeSupport::add_root(Eigen::Vector3d root_xyz)
{
    is_roots_.push_back(true);
    return add_node(root_xyz);
}

unsigned long RenderingTreeSupport::add_node(Eigen::Vector3d node_xyz)
{
    assert(root_r >= 0);

    is_roots_.push_back(false);
    node_coordinates_.push_back(node_xyz);

    graph_.push_back(std::vector<unsigned long>());

    return node_coordinates_.size() - 1;
}


void RenderingTreeSupport::add_edge(unsigned long u0, unsigned u1)
{
    assert(u0 < size() && u1 < size());
    graph_[u0].push_back(u1);
    return;
}

void RenderingTreeSupport::generate_node_radius()
{
    size_t u, v;
    double l,a,k1 = 0.1, k2 = 0.1;
    for (size_t kd = 0; kd < size(); kd++) {
        u = kd;
        for (size_t ed = 0; ed < graph_[u].size(); ed++) {
            v = graph_[u][ed];
            l = (node_coordinates_[v] - node_coordinates_[u]).norm();
            a = std::acos((node_coordinates_[v] - node_coordinates_[u]).y() / l);
            node_radius_[u] = l * k1;
        }
    }
}

void RenderingTreeSupport::draw(Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> vV;

    node_index_.resize(size());
    node_radius_.clear();
    node_radius_.resize(size(), 0);
    generate_node_radius();

    for (size_t kd = 0; kd < size(); kd++) {
        std::cout << node_radius_[kd] << " ";
    }
    std::cout << std::endl;


    //std::cout << size() << std::endl;
    for (size_t kd = 0; kd < size(); kd++) {
        if (is_node_tip(kd)) {
            vV.push_back(node_coordinates_[kd]);
            node_index_[kd].push_back(vV.size() - 1);

        } else {
            Eigen::Vector3d p = node_coordinates_[kd];
            Eigen::Vector3d dp(0, 0, 0);
            for (size_t id = 0; id < settings_.polygon_vertices_num; id++) {
                dp = Eigen::Vector3d(cos(2.0 * settings_.PI * id / settings_.polygon_vertices_num) * node_radius_[kd],
                                     0,
                                     sin(2.0 * settings_.PI * id / settings_.polygon_vertices_num) * node_radius_[kd]);

                vV.push_back(p + dp);
                node_index_[kd].push_back(vV.size() - 1);
            }
        }
    }

    std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i>> vF;
    unsigned long u, v;
    for (size_t kd = 0; kd < size(); kd++) {
        u = kd;
        for (size_t ed = 0; ed < graph_[u].size(); ed++) {
            v = graph_[u][ed];
            Eigen::Vector3d pu = node_coordinates_[u];
            Eigen::Vector3d pv = node_coordinates_[v];
//            std::cout << ( acos((pv - pu)(1) / (pv - pu).norm())  / 3.1415926f * 180 > 31 ) << std::endl;
            if (is_node_root(u))
                generate_base(u, v, vF);

            if (is_node_tip(v))
                generate_tip(u, v, vF);
            else
                generate_trunk(u, v, vF);
        }
    }

    V = Eigen::MatrixXd::Zero(vV.size(), 3);
    F = Eigen::MatrixXi::Zero(vF.size(), 3);

    for(size_t id = 0; id < vV.size(); id++)
        V.row(id) = vV[id];
    for(size_t id = 0; id < vF.size(); id++)
        F.row(id) = vF[id];
    return;
}

void RenderingTreeSupport::generate_base(unsigned long u, unsigned long v, std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i>> &vF)
{
    //bottom
    for(size_t id = 1; id < settings_.polygon_vertices_num - 1; id++)
    {
        vF.push_back(Eigen::Vector3i(node_index_[u][0], node_index_[u][id], node_index_[u][id + 1]));
    }
}

void RenderingTreeSupport::generate_trunk(unsigned long u, unsigned long v, std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i>> &vF)
{
    //lateral
    for(size_t id = 0;id < settings_.polygon_vertices_num; id++)
    {
        unsigned  long u0 = node_index_[u][id];
        unsigned  long u1 = node_index_[u][(id + 1) % settings_.polygon_vertices_num];
        unsigned  long v0 = node_index_[v][id];
        unsigned  long v1 = node_index_[v][(id + 1) % settings_.polygon_vertices_num];
        vF.push_back(Eigen::Vector3i(u0, v0, v1));
        vF.push_back(Eigen::Vector3i(v1, u1, u0));
    }
}

void RenderingTreeSupport::generate_tip(unsigned long u, unsigned long v, std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i>> &vF)
{
    //cone
    for(size_t id = 0;id < settings_.polygon_vertices_num; id++)
    {
        unsigned  long u0 = node_index_[u][id];
        unsigned  long u1 = node_index_[u][(id + 1) % settings_.polygon_vertices_num];
        unsigned  long v0 = node_index_[v][0];
        vF.push_back(Eigen::Vector3i(u0, v0, u1));
    }
}

double RenderingTreeSupport::calculate_tree_length()
{
    double dist = 0;
    size_t u,v;
    for (size_t kd = 0; kd < size(); kd++) {
        u = kd;
        for (size_t ed = 0; ed < graph_[u].size(); ed++) {
            v = graph_[u][ed];
            dist += (node_coordinates_[v] - node_coordinates_[u]).norm();
        }
    }
    return dist;
}


#endif //SUPPORTER_RENDERING_TREE_SUPPORT_H
