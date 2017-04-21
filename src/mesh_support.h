//
// Created by 汪子琦 on 4/21/17.
//

#ifndef EXAMPLE_MESH_SUPPORT_H
#define EXAMPLE_MESH_SUPPORT_H

#include "mesh_layout.h"
#include "fermat_spirals.h"
#include <algorithm>

typedef struct tagConvexHullPoint
{
public:
    ClipperLib::IntPoint pi;
    Eigen::Vector2d pd;
    double yh;

public:
    Eigen::RowVector3d PT(){return Eigen::RowVector3d(pd.x(), yh, pd.y());}
}ConvexHullPoint;

bool MSTGraphNodeComparator(MSTGraphNode &n0, MSTGraphNode &n1)
{
    if(n0.e1 < n1.e1)
        return true;
    if(n0.e1 == n1.e1 && n0.t1 < n1.t1)
        return true;
    return false;
}

bool ConvexHullPoint_Yi_Increase(ConvexHullPoint& p1, ConvexHullPoint& p2)
{
    if(p1.pd.y() < p2.pd.y())
        return true;
    if(p1.pd.y() == p2.pd.y() && p1.pd.x() < p2.pd.x())
        return true;
    return false;
}

bool ConvexHullPoint_TurnLeft(ConvexHullPoint& p1, ConvexHullPoint& p2)
{
    double cross = p1.pd.x() * p2.pd.y() - p1.pd.y() * p2.pd.x();
    if(cross > 0)
        return true;
    if(cross == 0 && p1.pd.norm() < p2.pd.norm())
        return true;
    return false;
}


class MeshSupport
{
public:
    MeshSupport()
    {

    }

    MeshSupport(Settings &s)
    {
        settings = s;
    }

public:

    int number_layer(){return layer_slices.size();}

    double layer_height(int id){return id * settings.layer_height;}

    void level_set(std::vector<Fermat_Level_Set> &fermat);

    void sp_pin_construction(MeshSlicer &slicer, Eigen::MatrixXd &H);

    void convex_hull_construction(std::vector< std::vector<ConvexHullPoint>> &convex);

public:

    void fermat_spiral( std::vector<std::list<FermatEdge>> &path_layer, MeshSlicer &slicer, Eigen::MatrixXd &H);

    void minimum_spanning_tree(std::vector< Fermat_Level_Set> &fermat,
                               std::vector< std::vector<MSTGraphNode>> &T);

    void fermat_spiral(std::vector< Fermat_Level_Set> &fermat,
                       std::vector< std::vector<MSTGraphNode>> &T,
                       int u,
                       int e,
                       double t,
                       std::list<FermatEdge> &path);

    void find_tree_roots(std::vector<bool> &pin_visited,
                         std::vector< std::vector<MSTGraphNode>> &T,
                         std::vector<int> &roots);

    void convert_path_clipper(ClipperLib::Paths &clipper ,std::list<FermatEdge> &fermat);

    void convert_clipper_path(ClipperLib::Paths &clipper ,std::list<FermatEdge> &fermat, int layer);

    void clipper_path(ClipperLib::Paths &support, int layer);

    void connect(std::list<FermatEdge> &path1, std::list<FermatEdge> &path2);

private:

    Settings settings;

    std::vector<ClipperLib::Paths> sp_pin;

    Eigen::MatrixXd platform;

    std::vector<ClipperLib::Paths> layer_slices;
};

void MeshSupport::level_set(std::vector<Fermat_Level_Set> &fermat)
{
    assert(!sp_pin.empty());

    std::vector<std::vector<ConvexHullPoint>> convex;
    convex_hull_construction(convex);
    fermat.resize(settings.pillar_column * settings.pillar_row);
    for(int pin_id = 0; pin_id < convex.size(); pin_id++)
    {
        if(!convex[pin_id].empty())
        {
            ClipperLib::Path outside;
            for(int id = 0; id < convex[pin_id].size(); id++)
            {
                outside.push_back(convex[pin_id][id].pi);
            }

            fermat[pin_id] = Fermat_Level_Set(outside, platform(pin_id / settings.pillar_column, pin_id % settings.pillar_column));

        }
    }

    return;
}


void MeshSupport::sp_pin_construction(MeshSlicer &slicer,  Eigen::MatrixXd &H)
{
    //platform
    if(H.isZero())
    {
        Mesh_Layout layout;
        double dx, dz;
        layout.xy_layout(slicer, dx, dz, H);
        slicer.move_XY(dx, dz);
    }

    platform = H;
    slicer.get_slices(layer_slices);

    //sp_pin
    std::vector<ClipperLib::Paths> bottom_half;
    ClipperLib::Paths downward;
    slicer.get_bottom_half(bottom_half);
    downward = bottom_half[0];
    ClipperLib::Clipper clipper;
    for(int layer = 1; layer < bottom_half.size(); layer++)
    {
        if(!bottom_half[layer].empty())
        {
            clipper.AddPaths(bottom_half[layer], ClipperLib::ptSubject, true);
        }
    }
    clipper.Execute(ClipperLib::ctUnion, downward, ClipperLib::pftPositive, ClipperLib::pftPositive);

    sp_pin.resize(settings.pillar_column * settings.pillar_row);
    for(int ir = 0; ir < settings.pillar_row; ir++)
    {
        for(int ic = 0; ic < settings.pillar_column; ic++)
        {
            ClipperLib::Path square;
            int L =  settings.mm2int(settings.pad_size * ic);
            int R =  settings.mm2int(settings.pad_size * (ic + 1));
            int T =  settings.mm2int(settings.pad_size * ir);
            int B =  settings.mm2int(settings.pad_size * (ir + 1));
            square.push_back(ClipperLib::IntPoint(L, T));
            square.push_back(ClipperLib::IntPoint(R, T));
            square.push_back(ClipperLib::IntPoint(R, B));
            square.push_back(ClipperLib::IntPoint(L, B));

            ClipperLib::Paths solution;
            clipper.Clear();
            clipper.AddPath(square, ClipperLib::ptSubject, true);
            clipper.AddPaths(downward, ClipperLib::ptClip, true);
            clipper.Execute(ClipperLib::ctIntersection, solution, ClipperLib::pftPositive, ClipperLib::pftPositive);
            ClipperLib::SimplifyPolygons(solution);
            if(!solution.empty())
                sp_pin[ir * settings.pillar_column + ic] = solution;
        }
    }

}

void MeshSupport::convex_hull_construction(std::vector<std::vector<ConvexHullPoint>> &convex)
{
    assert(!sp_pin.empty());

    convex.clear();
    convex.resize(settings.pillar_row * settings.pillar_column);

    for(int row = 0; row < settings.pillar_row; row++)
    {
        for(int col = 0; col < settings.pillar_column; col++)
        {

            int pin_id = row * settings.pillar_column + col;
            std::vector<ConvexHullPoint> pointslist;
            for(int id = 0; id < sp_pin[pin_id].size(); id++)
            {
                for(int jd = 0; jd < sp_pin[pin_id][id].size(); jd++)
                {
                    ConvexHullPoint cpoint;
                    cpoint.pd =  Eigen::Vector2d(settings.int2mm(sp_pin[pin_id][id][jd].X),
                                                 settings.int2mm(sp_pin[pin_id][id][jd].Y));
                    cpoint.pi = sp_pin[pin_id][id][jd];
                    cpoint.yh = platform(row, col);
                    pointslist.push_back(cpoint);
                }
            }

            if(pointslist.size() <= 2) continue;

            std::sort(pointslist.begin(), pointslist.end(), ConvexHullPoint_Yi_Increase);

            for(int id = 1; id < pointslist.size(); id++)
                pointslist[id].pd -= pointslist[0].pd;

            std::sort(pointslist.begin() + 1, pointslist.end(), ConvexHullPoint_TurnLeft);

            for(int id = 1; id < pointslist.size(); id++)
                pointslist[id].pd += pointslist[0].pd;

            std::vector<int> stack;
            stack.push_back(0);

            for(int id = 1; id < pointslist.size(); id++)
            {
                Eigen::Vector2d p0, p1, p2 = pointslist[id].pd;
                int p1_id = 0;
                double cross = 0;
                do
                {
                    p1 = pointslist[stack.back()].pd;
                    p1_id = stack.back();
                    stack.pop_back();
                    if(!stack.empty())
                    {
                        p0 = pointslist[stack.back()].pd;
                        cross = (p1 - p0).x() * (p2 - p1).y() - (p1 - p0).y() * (p2 - p1).x();
                    } else
                    {
                        break;
                    }
                }while(cross < 0);
                stack.push_back(p1_id);
                stack.push_back(id);
            }

            for(int id = 0; id < stack.size(); id++)
                convex[pin_id].push_back(pointslist[stack[id]]);
        }
    }

    return;
}

void MeshSupport::find_tree_roots(std::vector<bool> &pin_visited,
                                             std::vector< std::vector<MSTGraphNode>> &T,
                                             std::vector<int> &roots)
{
    std::vector<bool> is_roots;
    is_roots = pin_visited;

    for(int id = 0; id < T.size(); id++)
    {
        for(int jd = 0; jd < T[id].size(); jd++)
        {
            is_roots[T[id][jd].x] = false;
        }
    }

    for(int id = 0; id < T.size(); id++)
    {
        if(is_roots[id]) roots.push_back(id);
    }

    return;
}


void MeshSupport::convert_path_clipper(ClipperLib::Paths &clipper ,std::list<FermatEdge> &fermat)
{
    clipper.clear();
    ClipperLib::Path support;
    std::list<FermatEdge>::iterator it;
    for(it = fermat.begin(); it != fermat.end(); it++)
        support.push_back(ClipperLib::IntPoint(settings.mm2int(it->p0.x()), settings.mm2int(it->p0.y())));
    support.push_back(ClipperLib::IntPoint(settings.mm2int(fermat.back().p1.x()), settings.mm2int(fermat.back().p1.y())));
    clipper.push_back(support);
    return;
}

void MeshSupport::convert_clipper_path(ClipperLib::Paths &clipper ,std::list<FermatEdge> &fermat, int layer)
{
    fermat.clear();
    for(int id = 0; id < clipper.size(); id++) {
        for (int jd = 1; jd < clipper[id].size(); jd++) {
            FermatEdge fe;
            fe.p0 = Eigen::Vector2d(settings.int2mm(clipper[id][jd - 1].X),
                                    settings.int2mm(clipper[id][jd - 1].Y));
            fe.p1 = Eigen::Vector2d(settings.int2mm(clipper[id][jd].X),
                                    settings.int2mm(clipper[id][jd].Y));
            fe.yh = layer_height(layer);
            fermat.push_back(fe);
        }
    }
}

void MeshSupport::clipper_path(ClipperLib::Paths &support, int layer)
{

    ClipperLib::Paths mesh;
    mesh = layer_slices[layer];

    ClipperLib::Clipper clipper;
    ClipperLib::ClipperOffset co;
    co.AddPaths(mesh, ClipperLib::jtRound, ClipperLib::etClosedPolygon);
    co.Execute(mesh, settings.mm2int(settings.overhang_offset));
    ClipperLib::SimplifyPolygons(mesh);

    clipper.AddPaths(support, ClipperLib::ptSubject, false);
    clipper.AddPaths(mesh, ClipperLib::ptClip, true);

    ClipperLib::PolyTree polytree;
    clipper.Execute(ClipperLib::ctDifference, polytree, ClipperLib::pftPositive, ClipperLib::pftPositive);

    ClipperLib::Paths solution;
    ClipperLib::OpenPathsFromPolyTree(polytree, solution);

    support = solution;
    return;
}

void MeshSupport::fermat_spiral(std::vector< Fermat_Level_Set> &fermat,
                                           std::vector< std::vector<MSTGraphNode>> &T,
                                           int u,
                                           int e,
                                           double t,
                                           std::list<FermatEdge> &path)
{
    if(T[u].empty())
    {
        std::list< FermatEdge> present_path;
        fermat[u].fermat_spiral(e, t, present_path);
        connect(path, present_path);
        return;
    }
    else
    {
        bool is_root = t < 0;

        std::vector<MSTGraphNode> child;
        for(int id = 0; id < T[u].size(); id++)
            child.push_back(T[u][id]);

        if(is_root)
            fermat[u].anticlockwise(0, child[0].e1, child[0].t1, settings.fermat_cut_width * 2, e, t);

        for(int id = 0; id < child.size(); id++) {
            child[id].e1 = (child[id].e1 - e + fermat[u].num_level_id(0)) % fermat[u].num_level_id(0);
        }

        std::sort(child.begin(), child.end(), MSTGraphNodeComparator);

        for(int id = 0; id < child.size(); id++)
            child[id].e1  = (child[id].e1 + e) % fermat[u].num_level_id(0);

        std::vector< std::list< FermatEdge>> path_cuts;
        fermat[u].fermat_spiral(e, t, path_cuts, child);

        for(int id = 0; id < child.size(); id++)
        {
            connect(path, path_cuts[id]);
            std::list<FermatEdge> dfs_fermat;
            fermat_spiral(fermat, T, child[id].x, child[id].e2, child[id].t2, dfs_fermat);
            connect(path, dfs_fermat);
        }
        connect(path, path_cuts.back());

    }
}

void MeshSupport::fermat_spiral(std::vector<std::list<FermatEdge>> &path_layer, MeshSlicer &slicer, Eigen::MatrixXd &H)
{
    sp_pin_construction(slicer, H);

    std::vector<Fermat_Level_Set> fermat;
    level_set(fermat);

    //testing code
//    Fermat_Level_Set squre0, squre1, squre2, squre3;
//
//    ClipperLib::Path outer0;
//    outer0.push_back(ClipperLib::IntPoint(0, 0));
//    outer0.push_back(ClipperLib::IntPoint(0, 10000));
//    outer0.push_back(ClipperLib::IntPoint(10000, 10000));
//    outer0.push_back(ClipperLib::IntPoint(10000, 0));
//    squre0 = Fermat_Level_Set(outer0, 0);
//
//    ClipperLib::Path outer1;
//    outer1.push_back(ClipperLib::IntPoint(15000, 0));
//    outer1.push_back(ClipperLib::IntPoint(15000, 10000));
//    outer1.push_back(ClipperLib::IntPoint(25000, 10000));
//    outer1.push_back(ClipperLib::IntPoint(25000, 0));
//    squre1 = Fermat_Level_Set(outer1, 0);
//
//    ClipperLib::Path outer2;
//    outer2.push_back(ClipperLib::IntPoint(30000, 0));
//    outer2.push_back(ClipperLib::IntPoint(30000, 10000));
//    outer2.push_back(ClipperLib::IntPoint(40000, 10000));
//    outer2.push_back(ClipperLib::IntPoint(40000, 0));
//    squre2 = Fermat_Level_Set(outer2, 0);
//
//    ClipperLib::Path outer3;
//    outer3.push_back(ClipperLib::IntPoint(15000, 15000));
//    outer3.push_back(ClipperLib::IntPoint(15000, 25000));
//    outer3.push_back(ClipperLib::IntPoint(25000, 25000));
//    outer3.push_back(ClipperLib::IntPoint(25000, 15000));
//    squre3 = Fermat_Level_Set(outer3, 0);
//
//    fermat.push_back(squre0);
//    fermat.push_back(squre1);
//    fermat.push_back(squre2);
//    fermat.push_back(squre3);
    //testing code

    std::vector< std::vector<MSTGraphNode>> T;
    minimum_spanning_tree(fermat, T);

    //testing code
//    T.resize(4);
//    MSTGraphNode node;
//    roots.push_back(1);
//
//    squre1.closest(squre0, node.e1, node.t1, node.e2, node.t2);
//    node.x = 0;
//    T[1].push_back(node);
//
//    squre1.closest(squre2, node.e1, node.t1, node.e2, node.t2);
//    node.x = 2;
//    T[1].push_back(node);
//
//    squre1.closest(squre3, node.e1, node.t1, node.e2, node.t2);
//    node.x = 3;
//    T[1].push_back(node);
    //testing code


    std::vector<ClipperLib::Paths> clipper_layer;
    std::list<FermatEdge>::iterator it;
    path_layer.resize(number_layer());

    settings.print_N();
    settings.print_TsN("WALL SUPPORT");

    for(int layer_id = 1; layer_id < number_layer(); layer_id++)
    {
        memset(settings.tmp_str, 0, sizeof(settings.tmp_str));
        sprintf(settings.tmp_str, "\t layer %d, total %.3f %%...", layer_id,  100.0f * (double)(layer_id) / (layer_slices.size() - 2));
        settings.print_Ts(settings.tmp_str);

        //a vector for all the pin's height below P[layer_id]
        std::vector<bool> pin_below_layer;
        pin_below_layer.resize(settings.pillar_column * settings.pillar_row, false);

        //if there is a new pin which height is below P[layer_id], the pattern will be changed.
        bool fermat_pattern_change = false;

        for(int pin_id = 0; pin_id < settings.pillar_column * settings.pillar_row; pin_id++)
        {
            int row = pin_id / settings.pillar_column;
            int column = pin_id % settings.pillar_column;
            if(platform(row, column) < layer_height(layer_id) && platform(row, column) > 0)
            {
                pin_below_layer[pin_id] = true;
                if(platform(row, column) >= layer_height(layer_id - 1))
                    fermat_pattern_change = true;
            }

        }

        if(fermat_pattern_change)
        {

            std::vector<int> roots;
            std::vector< std::vector<MSTGraphNode>> T_layer;

            T_layer.resize(T.size());
            for(int id = 0; id < T.size(); id++) {
                if(!pin_below_layer[id]) continue;
                for (int jd = 0; jd < T[id].size(); jd++) {
                    if(pin_below_layer[T[id][jd].x])
                    {
                        T_layer[id].push_back(T[id][jd]);
                    }
                }
            }

            clipper_layer.clear();
            find_tree_roots(pin_below_layer, T_layer, roots);
            for(int id = 0; id < roots.size(); id++)
            {
                std::list<FermatEdge> subtree_path;
                fermat_spiral(fermat, T_layer, roots[id], -1, -1, subtree_path);

                ClipperLib::Paths subtree_clipper;
                convert_path_clipper(subtree_clipper, subtree_path);
                clipper_path(subtree_clipper, layer_id);
                clipper_layer.push_back(subtree_clipper);
            }
        }
        else
        {
            //if the fermat pattern didn't change,
            // the previous layer's fermal spiral has to subtract the new layer's slice polygons
            for(int id = 0; id < clipper_layer.size(); id++)
                if(clipper_layer[id].size() > 0)
                    clipper_path(clipper_layer[id], layer_id);
        }

        for(int id = 0; id < clipper_layer.size(); id++)
            if(clipper_layer[id].size() > 0)
            {
                std::list<FermatEdge> subtree_path;
                convert_clipper_path(clipper_layer[id], subtree_path, layer_id);
                path_layer[layer_id].insert(path_layer[layer_id].end(), subtree_path.begin(), subtree_path.end());
            }

        settings.print_TsN("done");
    }
    return;
}

void MeshSupport::minimum_spanning_tree(std::vector< Fermat_Level_Set> &fermat, std::vector< std::vector<MSTGraphNode>> &T)
{
    std::vector< std::vector<MSTGraphNode>> G;
    T.resize(fermat.size());
    G.resize(fermat.size());

    int num_fermats = 0;

    for(int pin_id = 0; pin_id < fermat.size(); pin_id++)
    {
        if(!fermat[pin_id].empty())
        {
            std::vector<int> neighbor;
            if(pin_id % settings.pillar_column != 0 && !fermat[pin_id - 1].empty())
                neighbor.push_back(pin_id - 1);
            if((pin_id + 1) % settings.pillar_column != 0 && !fermat[pin_id + 1].empty())
                neighbor.push_back(pin_id + 1);
            if(pin_id >= settings.pillar_column && !fermat[pin_id - settings.pillar_column].empty())
                neighbor.push_back(pin_id - settings.pillar_column);
            if(pin_id < settings.pillar_column * (settings.pillar_row - 1) && !fermat[pin_id + settings.pillar_column].empty())
                neighbor.push_back(pin_id + settings.pillar_column);

            for(int kd = 0; kd < neighbor.size(); kd++)
            {
                MSTGraphNode node;
                node.key = fermat[pin_id].closest(fermat[neighbor[kd]], node.e1, node.t1, node.e2, node.t2);
                node.x = neighbor[kd];
                G[pin_id].push_back(node);
            }

            num_fermats++;
        }
    }


    std::vector<double> dist;
    std::vector<int> pin;
    std::vector<bool> visited;

    dist.resize(fermat.size());
    pin.resize(fermat.size());
    visited.resize(fermat.size());

    for(int id = 0; id < dist.size(); id++) {
        dist[id] = settings.MAX_DOUBLE;
        pin[id] = -1;
        visited[id] = false;
    }

    do
    {
        bool non_valid = true;
        double minimum = settings.MAX_DOUBLE;
        int min_idx;
        for(int id = 0;id < dist.size(); id++)
        {
            if(!visited[id] && (minimum > dist[id])) {
                non_valid = false;
                min_idx = id;
                minimum = dist[id];
            }
        }

        if(non_valid)
        {
            for(int id = 0; id < dist.size(); id++)
            {
                if(!visited[id] && !fermat[id].empty()) {
                    dist[id] = 0;
                    pin[id] = id;
                    break;
                }
            }

        }
        else
        {
            visited[min_idx] = true;
            MSTGraphNode node;
            for(int id = 0; id < G[pin[min_idx]].size(); id++)
            {
                node = G[pin[min_idx]][id];
                if(node.x == min_idx)
                {
                    T[pin[min_idx]].push_back(node);
                    break;
                }
            }

            pin[min_idx] = min_idx;
            dist[min_idx] = 0;

            for(int id = 0; id < G[min_idx].size(); id++)
            {
                node = G[min_idx][id];
                if(dist[node.x] > node.key)
                {
                    dist[node.x] = node.key;
                    pin[node.x] = min_idx;
                }
            }

            num_fermats --;
        }
    }while(num_fermats);

    settings.print_N();
    settings.print_TsN("Minimal Spanning Tree: ");
    for(int id = 0; id < T.size(); id++)
    {
        for(int jd = 0; jd < T[id].size(); jd++)
        {
            memset(settings.tmp_str, 0, sizeof(settings.tmp_str));
            sprintf(settings.tmp_str, "\t Tree %d -> Tree %d = %f",id,  T[id][jd].x, T[id][jd].key);
            settings.print_TsN(settings.tmp_str);
        }
    }
    settings.print_N();
    return;
}

void MeshSupport::connect(std::list<FermatEdge> &path1, std::list<FermatEdge> &path2)
{
    if(!path1.empty() && !path2.empty())
    {
        FermatEdge fermatEdge = path1.back();
        fermatEdge.p0 = fermatEdge.p1;
        fermatEdge.p1 = path2.front().p0;
        path1.push_back(fermatEdge);
    }

    path1.insert(path1.end(), path2.begin(), path2.end());

    return;
}

#endif //EXAMPLE_MESH_SUPPORT_H
