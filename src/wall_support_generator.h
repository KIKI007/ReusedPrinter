//
// Created by 汪子琦 on 4/2/17.
//

#ifndef SUPPORTER_WALL_SUPPORT_GENERATOR_H
#define SUPPORTER_WALL_SUPPORT_GENERATOR_H

#include "slice_overhang_detector.h"
#include "fermat_spirals.h"
#include <algorithm>

typedef struct tagSprtPoint
{
    Eigen::Vector2d pt;
    int layer;
    int idx;
}SprtPoint;

bool MSTGraphNodeComparator(MSTGraphNode &n0, MSTGraphNode &n1)
{
    if(n0.e1 < n1.e1)
        return true;
    if(n0.e1 == n1.e1 && n0.t1 < n1.t1)
        return true;
    return false;
}

bool SprtXiComparatorIncrease(SprtPoint& p1, SprtPoint& p2)
{
    return (p1.pt[0] < p2.pt[0]);
}

bool SprtXiComparatorDecrease(SprtPoint& p1, SprtPoint& p2)
{
    return (p1.pt[0] > p2.pt[0]);
}

bool SprtYiComparatorIncrease(SprtPoint& p1, SprtPoint& p2)
{
    if(p1.pt.y() < p2.pt.y())
        return true;
    if(p1.pt.y() == p2.pt.y() && p1.pt.x() < p2.pt.x())
        return true;
    return false;
}

bool SprtCompratorTurnLeft(SprtPoint& p1, SprtPoint& p2)
{
    double cross = p1.pt.x() * p2.pt.y() - p1.pt.y() * p2.pt.x();
    if(cross > 0)
        return true;
    if(cross == 0 && p1.pt.norm() < p2.pt.norm())
        return true;
    return false;
}

class Wall_Support_Generator : public Slice_Overhang_Detector
{
public:

    Wall_Support_Generator()
    :Slice_Overhang_Detector()
    {
        platform = Eigen::MatrixXd::Zero(settings.pillar_row, settings.pillar_column);
    }

public:

    int get_num_vertices_before(int idx);
    void point_in_rectangle(Eigen::Vector2d point, ClipperLib::Path rectangle);
    bool point_in_polygons(Eigen::Vector2d p0, ClipperLib::Paths poly);

    void draw_platform(Eigen::MatrixXd &H);
    void draw_sp(Eigen::MatrixXd &SP);
    void draw_sp_lines(Eigen::MatrixXi &E);
    void draw_support(Eigen::MatrixXd &V, Eigen::MatrixXi &F);

    void convex_hull(std::vector< std::vector< int>> &convex_polys);
    void level_set(std::vector< Fermat_Level_Set> &fermat);

    void fermat_spiral( std::list<FermatEdge> &path);
    void minimum_spanning_tree(std::vector< Fermat_Level_Set> &fermat,
                               std::vector< std::vector<MSTGraphNode>> &T,
                               std::vector<int> &roots);

    void fermat_spiral(std::vector< Fermat_Level_Set> &fermat,
                       std::vector< std::vector<MSTGraphNode>> &T,
                       int u,
                       int e,
                       double t,
                       std::list<FermatEdge> &path);
private:

    void reverse(std::list< FermatEdge> &path);
    void connect(std::list<FermatEdge> &path1, std::list<FermatEdge> &path2);

    void projecting_sp_into_pin();
    void setup_platform();
    bool turnLeft(Eigen::Vector2d pA, Eigen::Vector2d pB, Eigen::Vector2d pt);
    void connecting_points(std::list<SprtPoint> &sprt);

private:
    std::vector<int> num_vertices_before;
    Eigen::MatrixXd platform;
    int valid_sp_size;
    std::vector <std::vector<Eigen::Vector2d>> sp_layer; // sample points for each layers
    std::vector <std::vector<Eigen::Vector2d>> sp_pin;   // sample points for each pin
    std::vector <std::vector<int>> sp_height_pin;        // height of sample points for each pin
};

void Wall_Support_Generator::connect(std::list<FermatEdge> &path1, std::list<FermatEdge> &path2)
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

int Wall_Support_Generator::get_num_vertices_before(int idx)
{
    assert(!sp_pin.empty());

    if(num_vertices_before.empty())
    {
        num_vertices_before.resize(sp_pin.size());
        num_vertices_before[0] = 0;
        for(int id = 1; id < sp_pin.size();id++)
            num_vertices_before[id] += num_vertices_before[id - 1] + sp_pin[id - 1].size();
    }
    assert(idx >= 0 && idx < sp_pin.size());

    return num_vertices_before[idx];
}

bool Wall_Support_Generator::point_in_polygons(Eigen::Vector2d p0, ClipperLib::Paths poly)
{
    if(point_on_polygons_edges(p0.x(), p0.y(), poly))
        return false;

    double wing_angle = 0;
    for(int id = 0; id < poly.size(); id++) {
        for (int jd = 0; jd < poly[id].size(); jd++) {
            Eigen::Vector2d p1(settings.int2mm(poly[id][jd].X),
                               settings.int2mm(poly[id][jd].Y));
            Eigen::Vector2d p2(settings.int2mm(poly[id][(jd + 1) % poly[id].size()].X),
                               settings.int2mm(poly[id][(jd + 1) % poly[id].size()].Y));
            wing_angle += anti_clockwise_angle(p1 - p0, p2 - p0);
        }
    }

    //std::cout << std::round(wing_angle / (2 * settings.PI)) << std::endl;
    if(std::round(wing_angle / (2 * settings.PI)) == 1)
        return true;

    return false;
}

//void Wall_Support_Generator::fermat_spiral(std::vector< Fermat_Level_Set> &fermat,
//                   std::vector< std::vector<MSTGraphNode>> &T,
//                   int u,
//                   int e,
//                   int t,
//                   std::list<FermatEdge> &path);

void Wall_Support_Generator::fermat_spiral(std::vector< Fermat_Level_Set> &fermat,
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

       std::sort(child.begin(), child.end(), MSTGraphNodeComparator);

       std::vector< std::list< FermatEdge>> path_cuts;

       if(is_root)
           fermat[u].anticlockwise(0, child[0].e1, child[0].t1, settings.fermat_cut_width * 2, e, t);

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

void Wall_Support_Generator::fermat_spiral(std::list<FermatEdge> &path)
{
    if(sp_pin.empty())
        setup_platform();

    std::vector<Fermat_Level_Set> fermat;
    level_set(fermat);

//    Fermat_Level_Set squre0, squre1, squre2;
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
//    fermat.push_back(squre0);
//    fermat.push_back(squre1);
//    fermat.push_back(squre2);


//    for(int id = 0; id < fermat.size(); id++)
//    {
//        if(fermat[id].num_level() != 0)
//        {
//            std::cout << "fermat :" << id << std::endl;
//            std::list<FermatEdge> present_path;
//            fermat[id].fermat_spiral(0, 0.5, present_path);
//            path.insert(path.end(), present_path.begin(), present_path.end());
//        }
//    }

    std::vector< std::vector<MSTGraphNode>> T;
    std::vector<int> roots;
    minimum_spanning_tree(fermat, T, roots);

//    T.resize(3);
//
//    MSTGraphNode node;
//    squre0.closest(squre1, node.e1, node.t1, node.e2, node.t2);
//
//    node.x = 1;
//    T[0].push_back(node);
//    roots.push_back(0);
//
//    squre1.closest(squre2, node.e1, node.t1, node.e2, node.t2);
//    node.x = 2;
//    T[1].push_back(node);

    for(int id = 0; id < roots.size(); id++)
    {
        std::list<FermatEdge> present_path;
        fermat_spiral(fermat, T, roots[id], -1, -1, present_path);
        path.insert(path.end(), present_path.begin(), present_path.end());
    }

    return;
}

void Wall_Support_Generator::convex_hull(std::vector<std::vector<int>> &convex_polys) {

    if(sp_pin.empty())
        setup_platform();

    convex_polys.clear();
    convex_polys.resize(settings.pillar_row * settings.pillar_column);

    for(int row = 0; row < settings.pillar_row; row++)
    {
        for(int col = 0; col < settings.pillar_column; col++)
        {

            int pin_id = row * settings.pillar_column + col;
            if(sp_pin[pin_id].size() <= 2) continue;

            std::vector<SprtPoint> sprt_list;
            for(int sp_id = 0; sp_id < sp_pin[pin_id].size(); sp_id++)
            {
                SprtPoint sprt;
                sprt.idx = sp_id;
                sprt.layer = sp_height_pin[pin_id][sp_id];
                sprt.pt = sp_pin[pin_id][sp_id];
                sprt_list.push_back(sprt);
            }

//            for(int id = 0; id < sprt_list.size(); id++)
//            {
//                std::cout << sprt_list[id].pt[0] << ", " << sprt_list[id].pt[1]  << std::endl;
//            }

            std::sort(sprt_list.begin(), sprt_list.end(), SprtYiComparatorIncrease);

            for(int id = 1; id < sprt_list.size(); id++)
                sprt_list[id].pt -= sprt_list[0].pt;

//            for(int id = 0; id < sprt_list.size(); id++)
//                std::cout << sprt_list[id].pt[0] << ", " << sprt_list[id].pt[1]  << std::endl;

            std::sort(sprt_list.begin() + 1, sprt_list.end(), SprtCompratorTurnLeft);

            for(int id = 1; id < sprt_list.size(); id++)
                sprt_list[id].pt += sprt_list[0].pt;

            std::vector<int> stack;
            stack.push_back(0);

            for(int id = 1; id < sprt_list.size(); id++)
            {
                Eigen::Vector2d p0, p1, p2 = sprt_list[id].pt;
                int p1_sprt_id;
                double cross = 0;
                do
                {
                    p1 = sprt_list[stack.back()].pt;
                    p1_sprt_id = stack.back();
                    stack.pop_back();
                    if(stack.empty()) break;
                    p0 = sprt_list[stack.back()].pt;
                    cross = (p1 - p0).x() * (p2 - p1).y() - (p1 - p0).y() * (p2 - p1).x();
                }while(cross < 0);
                stack.push_back(p1_sprt_id);
                stack.push_back(id);
            }

            for(int id = 0; id < stack.size(); id++)
                convex_polys[pin_id].push_back(sprt_list[stack[id]].idx);
        }
    }

    return;
}

void Wall_Support_Generator::minimum_spanning_tree(std::vector< Fermat_Level_Set> &fermat,
                                                   std::vector< std::vector<MSTGraphNode>> &T,
                                                   std::vector<int> &roots)
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

    std::vector<bool> is_roots;
    is_roots.resize(fermat.size(), false);
    for(int id = 0; id < T.size(); id++)
    {
        if(!fermat[id].empty()) is_roots[id] = true;
    }

    for(int id = 0; id < T.size(); id++)
    {
        for(int jd = 0; jd < T[id].size(); jd++)
        {
            std::cout << id << "->" << T[id][jd].x << " " <<  T[id][jd].key << std::endl;
            is_roots[T[id][jd].x] = false;
        }
    }

    for(int id = 0; id < T.size(); id++)
    {
        if(is_roots[id]) roots.push_back(id);
    }
    return;
}

void Wall_Support_Generator::level_set( std::vector<Fermat_Level_Set> &fermat)
{
    if(sp_pin.empty())
        setup_platform();

    std::vector<std::vector<int>> convex_polys;
    convex_hull(convex_polys);
    fermat.resize(convex_polys.size());
    for(int pin_id = 0; pin_id < convex_polys.size(); pin_id++)
    {
        if(!convex_polys[pin_id].empty())
        {
            ClipperLib::Path outside;
            for(int id = 0; id < convex_polys[pin_id].size(); id++)
            {
                int p0_id = convex_polys[pin_id][id];
                outside.push_back(ClipperLib::IntPoint(settings.mm2int(sp_pin[pin_id][p0_id].x()),
                                                       settings.mm2int(sp_pin[pin_id][p0_id].y())));
            }

            fermat[pin_id] = Fermat_Level_Set(outside, platform(pin_id / settings.pillar_column, pin_id % settings.pillar_column));

        }
    }

    return;
}

void Wall_Support_Generator::projecting_sp_into_pin()
{
    sp_pin.resize(settings.pillar_row * settings.pillar_column);
    sp_height_pin.resize(settings.pillar_row * settings.pillar_column);

    //valid means the sample points are supported by metal pin
    //non-valid means the sample points are supported by mesh itself
    valid_sp_size = 0;
    sampling(sp_layer);

    //if a sample point intersect with these polygons, this sample point will not be a valid point
    ClipperLib::Paths downward_polys;
    for(int layer_id = 1; layer_id < layer_slices.size(); layer_id ++)
    {

        for(int sp_id = 0; sp_id < sp_layer[layer_id].size(); sp_id++)
        {
            Eigen::Vector2d sp = sp_layer[layer_id][sp_id];
            if(point_in_polygons(sp, downward_polys) == false)
            {
                // projecting the sample points into its corresponding pin
                int pin_id = (int)(sp.y() / settings.pad_size) * settings.pillar_column + (int)(sp.x() / settings.pad_size);
                sp_pin[pin_id].push_back(sp);
                sp_height_pin[pin_id].push_back(layer_id);
                valid_sp_size ++;
            }
        }

        //downward_polys will accumulate durning iteration
        //downward_polys = downward_polys + layer_silces[layer_id]
        ClipperLib::Paths subj, clip;
        subj = downward_polys;
        clip = layer_slices[layer_id];

        ClipperLib::Clipper clipper;
        clipper.AddPaths(subj, ClipperLib::ptSubject, true);
        clipper.AddPaths(clip, ClipperLib::ptClip, true);

        clipper.Execute(ClipperLib::ctUnion, downward_polys, ClipperLib::pftPositive, ClipperLib::pftPositive);
        ClipperLib::SimplifyPolygons(downward_polys);
    }

    return;
}

void Wall_Support_Generator::draw_platform(Eigen::MatrixXd &H) {

    if(platform.isZero())
        setup_platform();

    H = platform;
    return;
}

void Wall_Support_Generator::draw_sp(Eigen::MatrixXd &SP)
{
    if(platform.isZero())
        setup_platform();

    int jd = 0;
    SP.resize(valid_sp_size, 3);
    for (int row = 0; row < settings.pillar_row; row++) {
        for (int col = 0; col < settings.pillar_column; col++) {
            int pin_id = row * settings.pillar_column + col;
            for(int sp_id = 0; sp_id < sp_pin[pin_id].size(); sp_id++)
            {
                SP(jd, 0) = sp_pin[pin_id][sp_id][0];
                SP(jd, 1) = platform(row, col);
                SP(jd, 2) = sp_pin[pin_id][sp_id][1];
                jd++;
            }
        }
    }

    return;
}

void Wall_Support_Generator::setup_platform() {

    if(sp_pin.empty())
        projecting_sp_into_pin();

    //compute each pin's height
    for (int row = 0; row < settings.pillar_row; row++)
    {
        for (int col = 0; col < settings.pillar_column; col++)
        {
            //if there is not sample points on this pin, set pin's height to 0
            if(sp_pin[row * settings.pillar_column + col].empty())
            {
                platform(row, col) = 0;
                continue;
            }

            //create a square which represents a metal pin
            ClipperLib::Path square;
            int sq_width =  settings.mm2int(settings.pad_size);
            square.push_back(ClipperLib::IntPoint(col * sq_width, row * sq_width));
            square.push_back(ClipperLib::IntPoint((col + 1) * sq_width, row * sq_width));
            square.push_back(ClipperLib::IntPoint((col + 1) * sq_width, (row + 1) * sq_width));
            square.push_back(ClipperLib::IntPoint(col * sq_width, (row + 1) * sq_width));

            // for each layer, computing the intersection between the slicing polygons and the square
            // if intersection existed, the height of this pin is P[id - 1]
            // else continue;
            int id = 0;
            for(;id < layer_slices.size(); id++)
            {
                ClipperLib::Paths subj, clip;
                subj.push_back(square);
                clip = layer_slices[id];

                ClipperLib::Clipper clipper;
                clipper.AddPaths(subj, ClipperLib::ptSubject, true);
                clipper.AddPaths(clip, ClipperLib::ptClip, true);

                //computing the intersection
                ClipperLib::Paths intersection;
                clipper.Execute(ClipperLib::ctIntersection, intersection, ClipperLib::pftPositive, ClipperLib::pftPositive);
                ClipperLib::SimplifyPolygons(intersection);

                //check the collision
                if(!intersection.empty())
                    break;
            }

            //!! since the length of each pin can only be the mutiple of settings.pillar_standard_height (present 0.5inch)
            //   draw down the pin to make it satisfy to the fabrication requirement
            if(id < layer_slices.size() && id > 0) {
                int num_standard_height = settings.int2mm(P[id - 1]) / settings.pillar_standard_height;
                platform(row, col) = (double)num_standard_height * settings.pillar_standard_height;
            }
            else
                platform(row, col) = 0;
        }
    }
}

void Wall_Support_Generator::draw_sp_lines(Eigen::MatrixXi &E)
{
    if(sp_pin.empty())
        setup_platform();

    std::vector<Eigen::RowVector2i> sp_lines;
    int prev_sp_nums = 0;
    for (int row = 0; row < settings.pillar_row; row++) {
        for (int col = 0; col < settings.pillar_column; col++) {
            int pin_id = row * settings.pillar_column + col;
            if (!sp_pin[pin_id].empty()) {

                std::list<SprtPoint> support_points;
                for(int sp_id = 0; sp_id < sp_pin[pin_id].size(); sp_id++)
                {
                    SprtPoint sprt;
                    sprt.pt = sp_pin[pin_id][sp_id];
                    sprt.idx = sp_id;
                    sprt.layer = sp_height_pin[pin_id][sp_id];
                    support_points.push_back(sprt);
                }

                connecting_points(support_points);
                std::list<SprtPoint>::iterator it, nit;
                it = support_points.begin();
                nit = std::next(it, 1);
                for(; nit != support_points.end(); ++it, ++nit)
                {
                    sp_lines.push_back(Eigen::RowVector2i(it->idx + prev_sp_nums, nit->idx + prev_sp_nums));
                }
            }
            prev_sp_nums += sp_pin[pin_id].size();
        }
    }

    E.resize(sp_lines.size(), 2);
    for(int id = 0; id < sp_lines.size(); id++)
        E.row(id) = sp_lines[id];

    return;
}

void Wall_Support_Generator::draw_support(Eigen::MatrixXd &V, Eigen::MatrixXi &F) {

    if(sp_pin.empty())
        setup_platform();

    int Vsize = 0;
    for (int row = 0; row < 9; row++) {
        for (int col = 0; col < settings.pillar_column; col++) {
            Vsize += sp_pin[row * settings.pillar_column + col].size();
        }
    }
    V.resize(2 * Vsize, 3);
    std::vector<Eigen::RowVector3i> fs;

    int size = 0;
    for (int row = 0; row < 9; row++) {
        for (int col = 0; col < settings.pillar_column; col++) {
            if (!sp_pin[row * settings.pillar_column + col].empty()) {

                std::list<SprtPoint> support_points;
                for (int id = 0; id < sp_pin[row * settings.pillar_column + col].size(); id++) {
                    SprtPoint sprt;
                    sprt.pt = sp_pin[row * settings.pillar_column + col][id];
                    sprt.idx = id;
                    sprt.layer = sp_height_pin[row * settings.pillar_column + col][id];
                    support_points.push_back(sprt);

                    V(size + id, 0) = sp_pin[row * settings.pillar_column + col][id][0];
                    V(size + id, 1) = platform(row, col);
                    V(size + id, 2) = sp_pin[row * settings.pillar_column + col][id][1];

                    V(size + id + Vsize, 0) = sp_pin[row * settings.pillar_column + col][id][0];
                    V(size + id + Vsize, 1) = settings.int2mm(P[sp_height_pin[row * settings.pillar_column + col][id]]);
                    V(size + id + Vsize, 2) = sp_pin[row * settings.pillar_column + col][id][1];
                }

                connecting_points(support_points);
                std::list<SprtPoint>::iterator it = support_points.begin(), nit = std::next(it, 1);
                for(; nit != support_points.end(); ++nit, it++)
                {
                    fs.push_back(Eigen::RowVector3i(nit -> idx + size, it->idx + size, nit -> idx + size + Vsize));
                    fs.push_back(Eigen::RowVector3i(it->idx + size + Vsize, nit -> idx + size + Vsize, it -> idx + size));
                }
            }
            size += sp_height_pin[row * settings.pillar_column + col].size();
        }
    }

    F.resize(fs.size(), 3);
    for(int id = 0; id < fs.size(); id++)
        F.row(id) = fs[id];

    return;
}

bool Wall_Support_Generator::turnLeft(Eigen::Vector2d pA, Eigen::Vector2d pB, Eigen::Vector2d pt)
{
    Eigen::Vector3d p1((pB - pA)(0), (pB - pA)(1), 0);
    Eigen::Vector3d p2((pt - pA)(0), (pt - pA)(1), 0);
    if(p1.cross(p2)[2] > 0)
        return true;
    else
        return false;
}

void Wall_Support_Generator::connecting_points(std::list<SprtPoint> &sprt)
{
    if(sprt.size() <= 1)
        return;

    sprt.sort(SprtXiComparatorIncrease);

    std::list<SprtPoint> upper, lower;
    std::list<SprtPoint>::iterator it;

    SprtPoint sA = sprt.front();
    SprtPoint sB = sprt.back();
    sprt.pop_front();sprt.pop_back();
    for(it = sprt.begin(); it != sprt.end(); it++)
    {
        Eigen::Vector2d pt = it->pt;
        if(turnLeft(sA.pt, sB.pt, pt))
            upper.push_back(*it);
        else
            lower.push_back(*it);
    }

    upper.sort(SprtXiComparatorDecrease);
    lower.sort(SprtXiComparatorIncrease);

    sprt.clear();
    sprt.push_back(sA);
    sprt.insert(sprt.end(), lower.begin(), lower.end());
    sprt.push_back(sB);
    sprt.insert(sprt.end(), upper.begin(), upper.end());

}


#endif //SUPPORTER_WALL_SUPPORT_GENERATOR_H
