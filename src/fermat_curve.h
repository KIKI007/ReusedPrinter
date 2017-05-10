//
// Created by 汪子琦 on 5/3/17.
//

#ifndef SUPPORTER_FERMAT_CURVE_H
#define SUPPORTER_FERMAT_CURVE_H

#include <vector>
#include <cmath>
#include "clipper.hpp"
#include "Settings.h"
#include <Eigen/Dense>

using std::vector;
using Eigen::RowVector2d;
using Eigen::Vector2d;
using ClipperLib::Path;
using ClipperLib::Paths;
using ClipperLib::IntPoint;


class FermatCurve
{
public:

    FermatCurve (Settings s)
    {
        settings = s;
    }

public:/*I/O*/
    bool set_boundary(Path &bdary, IntPoint sta, IntPoint end);

    void output_curve(Path &curve);

    bool is_line();

protected:

    bool set_boundary(int layer, IntPoint sta, IntPoint end);

    void rearrage_origin(int layer, double s);

    IntPoint get_point(int layer, double s);

    double get_arc_length(int layer, IntPoint &p);

    void shrink_boundary(Path &bdary);

    void get_path(Path &path, int layer, double s1, double s2);

    void connect(Path &path1, Path &path2);

    void reverse(Path &path);

protected: /*geometry*/

    double closed_path_length(Path &path);

    double open_path_length(Path &path);

    double Dist(IntPoint p, IntPoint q);

    IntPoint point(IntPoint p, IntPoint q, double ratio);

    IntPoint close_point(IntPoint p, IntPoint q, IntPoint m);

public:

    vector<Path> polygon_curve;

    vector< double> polygon_length;

    vector<Path> upper_curve, lower_curve;

    Settings settings;
};

bool FermatCurve::set_boundary(Path &bdary, IntPoint sta, IntPoint end)
{
    polygon_curve.push_back(bdary);
    polygon_length.push_back(closed_path_length(bdary));

    return set_boundary(0, sta, end);
}

bool FermatCurve::set_boundary(int layer, IntPoint sta, IntPoint end) {

    Path polygon = polygon_curve[layer];
    shrink_boundary(polygon);
    double cut_width = settings.fermat_cut_width;
    double next_level_length = closed_path_length(polygon);

    double s1 = get_arc_length(layer, sta);
    double s2 = get_arc_length(layer, end);
    sta = get_point(layer, s1);
    end = get_point(layer, s2);

    rearrage_origin(layer, s1);
    if(s2 > s1) s2 -= s1;
    else s2 = polygon_length[layer] - s1 + s2;

    double sta_end_length = s2;
    double end_sta_length = polygon_length[layer] - s2;
    Path path;

    if(next_level_length < settings.support_center_factor * cut_width)
    {
        if(layer == 0)
        {
            path = polygon_curve[layer];
            path.push_back(polygon_curve[layer].front());
            upper_curve.push_back(path);
            return false;
        }
        else {

            if (sta_end_length > end_sta_length) {
                get_path(path, layer, 0, sta_end_length);
                upper_curve.push_back(path);
            } else {
                get_path(path, layer, sta_end_length, polygon_length[layer]);
                lower_curve.push_back(path);
            }
            return true;
        }
    }
    else if(sta_end_length < 2 * cut_width)
    {
        sta_end_length += 2 * cut_width;
        end_sta_length -= 2 * cut_width;
    }
    else if(end_sta_length < 2 * cut_width)
    {
        end_sta_length += 2 * cut_width;
        sta_end_length -= 2 * cut_width;
    }

    get_path(path, layer, 0, sta_end_length - cut_width);
    upper_curve.push_back(path);

    get_path(path, layer, sta_end_length,  polygon_length[layer] - cut_width);

    lower_curve.push_back(path);


    polygon_curve.push_back(polygon);
    polygon_length.push_back(next_level_length);

    return set_boundary(layer + 1, sta, end);
}

void FermatCurve::output_curve(Path &curve)
{
    curve.clear();

    int num_circle = polygon_curve.size();
    int orientation = 1;
    Path path;
    for(int id = 0; id < num_circle - 1; id++, orientation ^= 1)
    {
        if(!orientation) connect(curve, lower_curve[id]);
        else connect(curve, upper_curve[id]);
    }

    if(upper_curve.size() > lower_curve.size())
        connect(curve, upper_curve.back());
    else{
        path = lower_curve.back();
        reverse(path);
        connect(curve, path);
    }


    for(int id = num_circle - 2; id >= 0; id--, orientation ^= 1)
    {
        path = lower_curve[id];

        if(!orientation) path = lower_curve[id];
        else  path = upper_curve[id];

        reverse(path);
        connect(curve, path);
    }

    return;
}

void FermatCurve::rearrage_origin(int layer, double s)
{
    assert(0 <= layer && layer < polygon_curve.size());
    assert(0 <= s && s <= polygon_length[layer]);

    Path polygon = polygon_curve[layer];
    Path new_polygon;
    for(int id = 0; id < polygon.size(); id++)
    {
        IntPoint p0 = polygon[id];
        IntPoint p1 = polygon[(id + 1) % polygon.size()];
        double edge_length = Dist(p0, p1);
        if(s < edge_length)
        {
            IntPoint q = point(p0, p1, s / edge_length);
            new_polygon.push_back(q);
            for(int jd = (id + 1) % polygon.size(); jd != id; jd = (jd + 1) % polygon.size())
            {
                new_polygon.push_back(polygon[jd]);
            }
            if(s > settings.ZERO_EPS) new_polygon.push_back(p0);
            polygon_curve[layer] = new_polygon;
            return;
        }
        s -= edge_length;
    }
}

IntPoint FermatCurve::get_point(int layer, double s)
{
    assert(0 <= layer && layer < polygon_curve.size());
    assert(0 <= s && s <= polygon_length[layer]);

    Path polygon = polygon_curve[layer];

    for(int id = 0; id < polygon.size(); id++)
    {
        IntPoint p0 = polygon[id];
        IntPoint p1 = polygon[(id + 1) % polygon.size()];
        double edge_length = Dist(p0, p1);
        if(s <= edge_length)
            return point(p0, p1, s / edge_length);
        s -= edge_length;
    }

    return IntPoint();
}

double FermatCurve::get_arc_length(int layer, IntPoint &p) {
    assert(0 <= layer && layer < polygon_curve.size());

    Path polygon = polygon_curve[layer];

    double s = 0;
    double minimum_length = settings.MAX_DOUBLE;
    double opt_s = 0;
    for(int id = 0; id < polygon.size(); id++)
    {
        IntPoint p0 = polygon[id];
        IntPoint p1 = polygon[(id + 1) % polygon.size()];
        IntPoint q = close_point(p0, p1, p);
        double dist = Dist(p, q);
        if(minimum_length > dist)
        {
            minimum_length = dist;
            opt_s = s + Dist(q, p0);
        }
        s += Dist(p0, p1);
    }

    return opt_s;
}

void FermatCurve::shrink_boundary(Path &bdary) {
    ClipperLib::ClipperOffset co;
    co.AddPath(bdary, ClipperLib::jtRound, ClipperLib::etClosedPolygon);
    ClipperLib::Paths sol;
    co.Execute(sol, -settings.mm2int(settings.extrusion_width * 2));

    ClipperLib::SimplifyPolygons(sol);
    bdary.clear();
    for(int id = 0; id < sol.size(); id++)
    {
        if(bdary.size() < sol[id].size())
            bdary = sol[id];
    }
}

void FermatCurve::get_path(Path &path, int layer, double s1, double s2) {
    assert(0 <= layer && layer < polygon_curve.size());
    assert(0 <= s1 && s1 <= polygon_length[layer]);
    assert(0 <= s2 && s2 <= polygon_length[layer]);
    assert(s1 <= s2);

    path.clear();
    if(s1 == s2)  return;

    Path polygon = polygon_curve[layer];
    for(int id = 0; id < polygon.size(); id++)
    {
        IntPoint p0 = polygon[id];
        IntPoint p1 = polygon[(id + 1) % polygon.size()];
        double dist = Dist(p0, p1);
        if(s1 < dist)
        {
            IntPoint q0 = point(p0, p1, s1 / dist);
            path.push_back(q0);
            if(s2 < dist)
            {
                IntPoint q1 = point(p0, p1, s2 / dist);
                path.push_back(q1);
                return;
            }
            else
            {
                s2 -= dist;
                for(int jd = (id + 1) % polygon.size(); ; jd = (jd + 1) % polygon.size())
                {
                    p0 = polygon[jd];
                    p1 = polygon[(jd + 1) % polygon.size()];
                    dist = Dist(p0, p1);
                    path.push_back(p0);

                    if(s2 <= dist)
                    {
                        IntPoint q1 = point(p0, p1, s2 / dist);
                        path.push_back(q1);
                        return;
                    }
                    s2 -= dist;
                }
            }

        }
        s1 -= dist;
        s2 -= dist;
    }

    return;
}

void FermatCurve::connect(Path &path1, Path &path2) {
    path1.insert(path1.end(), path2.begin(), path2.end());
}

void FermatCurve::reverse(Path &path) {
    Path r_path(path.size());
    for(int id = 0; id < path.size(); id++)
    {
        r_path[id] = path[path.size() - id - 1];
    }
    path = r_path;
    return;
}

double FermatCurve::closed_path_length(Path &path) {
    double length = 0;
    for(int id = 0; id < path.size(); id++)
    {
        length += Dist(path[(id + 1) % path.size()], path[id]);
    }
    return length;
}

double FermatCurve::open_path_length(Path &path) {
    double length = 0;
    for(int id = 0; id < path.size() - 1; id++)
    {
        length += Dist(path[id + 1], path[id]);
    }
    return length;
}

double FermatCurve::Dist(IntPoint p, IntPoint q) {

    double x = settings.int2mm(p.X - q.X);
    double y = settings.int2mm(p.Y - q.Y);
    return std::sqrt(x * x+y * y);
}

IntPoint FermatCurve::point(IntPoint p, IntPoint q, double ratio) {

     IntPoint m ;
     m.X = p.X + (q.X - p.X) * ratio;
     m.Y = p.Y + (q.Y - p.Y) * ratio;
     return m;
}

IntPoint FermatCurve::close_point(IntPoint p, IntPoint q, IntPoint m) {

    RowVector2d pq(settings.int2mm(q.X - p.X), settings.int2mm(q.Y - p.Y));
    RowVector2d pm(settings.int2mm(m.X - p.X), settings.int2mm(m.Y - p.Y));
    double ratio = pq.dot(pm) / Dist(p, q);

    if(ratio < 0) return p;
    else if(ratio > 1) return q;
    else return point(p, q, ratio);
}

bool FermatCurve::is_line() {
    if(polygon_curve.size() <= 1)
        return true;
    else
        return false;
}

#endif //SUPPORTER_FERMAT_CURVE_H
