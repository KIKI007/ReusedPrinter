//
// Created by 汪子琦 on 4/2/17.
//

#ifndef SUPPORTER_WALL_SUPPORT_GENERATOR_H
#define SUPPORTER_WALL_SUPPORT_GENERATOR_H

#include "slice_overhang_detector.h"
#include <algorithm>
typedef struct tagSprtPoint
{
    Eigen::Vector2d pt;
    int layer;
    int idx;
}SprtPoint;

typedef struct tagFermatEdge
{
    //edge is in clockwise
public:
    tagFermatEdge(ClipperLib::IntPoint in_p0, ClipperLib::IntPoint in_p1)
    {
        p0 = Eigen::Vector2d(settings.int2mm(in_p0.X), settings.int2mm(in_p0.Y));
        p1 = Eigen::Vector2d(settings.int2mm(in_p1.X), settings.int2mm(in_p1.Y));
    }

public:

    bool on(Eigen::Vector2d p)
    {
        double cross = (p - p0)[0] * (p - p1)[1] - (p - p0)[1] * (p - p1)[0];
        return true;
    }
public:
    Eigen::Vector2d p0, p1;
    Settings settings;
}FermatEdge;

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

private:
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
            if(sp_pin[pin_id].size() <= 1) continue;

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
