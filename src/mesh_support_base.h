//
// Created by 汪子琦 on 5/4/17.
//

#ifndef SUPPORTER_MESH_SUPPORT_BASE_H
#define SUPPORTER_MESH_SUPPORT_BASE_H

#include "mesh_slicer_base.h"
#include "fermat_curve.h"
#include <iostream>
#include <queue>

using ClipperLib::IntPoint;
using ClipperLib::Path;
using ClipperLib::Paths;
using Eigen::RowVector2i;
typedef  vector<IntPoint> PointsCluster;
typedef  vector< vector<IntPoint>> VecPointsCluster;
typedef  std::queue<RowVector2i> Queue;
using ClipperLib::Clipper;
using ClipperLib::pftPositive;
using ClipperLib::ctIntersection;
using ClipperLib::ctDifference;
using ClipperLib::ptSubject;
using ClipperLib::ptClip;

typedef long long INT64;

typedef struct tagConvexHullPoint {
public:

    tagConvexHullPoint() {
        pi = IntPoint(0, 0);
    }

    tagConvexHullPoint(IntPoint p) {
        pi = p;
    }

public:
    IntPoint pi;
}ConvexHullPoint;

class CurveGroup
{

public:
    CurveGroup(double lowest, double highest, Path &poly, bool line = false)
    {
        lowest_layer_height = lowest;
        highest_layer_height = highest;
        polygon = poly;
        is_line = line;
        set_center();
    }

public:

    void closet_point_polygon(IntPoint &p0, IntPoint &p1, CurveGroup &A)
    {
        Settings settings;
        p0 = IntPoint(0, 0);
        p1 = IntPoint(0, 0);
        double opt_dist = settings.MAX_DOUBLE;
        for(int id = 0; id < polygon.size(); id++)
        {
            for(int jd = 0; jd < A.polygon.size(); jd++)
            {
                IntPoint close;
                double dist = close_point(A.polygon[jd], A.polygon[(jd + 1) % A.polygon.size()], polygon[id], close);
                if(dist < opt_dist)
                {
                    p0 = polygon[id];
                    p1 = close;
                    opt_dist = dist;
                }
                if(A.is_line) break;
            }
        }
        return;
    }

public:

    double close_point(IntPoint p, IntPoint q, IntPoint m, IntPoint &close)
    {
        Settings settings;
        RowVector2d pq(settings.int2mm(q.X - p.X), settings.int2mm(q.Y - p.Y));
        RowVector2d pm(settings.int2mm(m.X - p.X), settings.int2mm(m.Y - p.Y));
        double ratio = pq.dot(pm) / Dist(p, q);

        if(ratio < 0) close = p;
        else if(ratio > 1) close = q;
        else close = point(p, q, ratio);

        return Dist(close, m);
    }

    IntPoint point(IntPoint p, IntPoint q, double ratio) {
        IntPoint m ;
        m.X = p.X + (q.X - p.X) * ratio;
        m.Y = p.Y + (q.Y - p.Y) * ratio;
        return m;
    }

    double Dist(IntPoint p, IntPoint q) {
        Settings settings;
        double x = settings.int2mm(p.X - q.X);
        double y = settings.int2mm(p.Y - q.Y);
        return std::sqrt(x * x+y * y);
    }

private:

    void set_center()
    {
        int center_x = 0, center_y = 0;
        for(int id = 0; id < polygon.size(); id++)
        {
            center_x += polygon[id].X;
            center_y += polygon[id].Y;
        }
        center_x /= polygon.size();
        center_y /= polygon.size();
    }
public:
    bool is_line;
    double lowest_layer_height;
    double highest_layer_height;
    Path polygon;

    double center_x;
    double center_y;
};

bool CurveGroupIncrease(CurveGroup &A, CurveGroup &B)
{
    if(A.center_x + A.center_y < B.center_x + B.center_y)
        return true;
    else
        return false;
}

bool ConvexHullPoint_Yi_Increase(ConvexHullPoint& p1, ConvexHullPoint& p2)
{
    if(p1.pi.Y < p2.pi.Y)
        return true;
    if(p1.pi.Y == p2.pi.Y && p1.pi.X < p2.pi.X)
        return true;
    return false;
}

IntPoint minus(IntPoint p, IntPoint q)
{
    return IntPoint(p.X - q.X, p.Y - q.Y);
}

IntPoint plus(IntPoint p, IntPoint q)
{
    return IntPoint(p.X + q.X, p.Y + q.Y);
}

INT64 cross(IntPoint p1, IntPoint p2)
{
    return (INT64)p1.X * (INT64)p2.Y - (INT64)p1.Y * (INT64)p2.X;
}

INT64 square_norm(IntPoint p)
{
    return (INT64) p.X * (INT64) p.X + (INT64) p.Y * (INT64) p.Y;
}

bool ConvexHullPoint_TurnLeft(ConvexHullPoint& p1, ConvexHullPoint& p2)
{
    INT64 cross_value = cross(p1.pi, p2.pi);
    if(cross_value > 0)
        return true;
    if(cross_value == 0 && square_norm(p1.pi) < square_norm(p2.pi))
        return true;
    return false;
}

class MeshSupportBase
{
public:
    MeshSupportBase()
    {

    }

public: //fermat curves
    void get_fermat_curves(vector<Paths> &fermat_curves);

    virtual void get_fermat_curves(int layer, Paths &fermat_curves);

    virtual void fermat_curves_construction();

    void fermat_curves_difference_from_model(int layer, Paths &paths);

public: //regroup and get convex hull

    void regroup_support_points(vector<VecPointsCluster> &pin_clusters);

    void regroup_support_points(int ir, int ic, vector<VecPointsCluster> &pin_clusters);

    void curve_group_construction();

    int convex_hull(PointsCluster &cluster, Path &polygon);

public: //utility function

    void set_slicer(MeshSlicerBase &slicer_base);

    void set_map(MatrixXi &hmap, MatrixXi &smap, MatrixXd &plfm);

    void maximum_pin_layer_height_construction(MatrixXd &maximum_height);

    bool is_cluster_line(PointsCluster &cluster);

    void cluster_bounding_box(PointsCluster &cluster, double &width, double &height);

    void singleline_curvegroup_construction(PointsCluster &cluster, Path &polygon);

    void clear();

    void shrink_support_map(MatrixXi &smap);

protected:

    MatrixXi height_map, support_map;

    MatrixXd platform;

    vector<CurveGroup> curve_groups;

    vector<Paths> groups_fermat_curves;

protected:

    MeshSlicerBase slicer;

    vector<Paths> slices;

    Settings settings;
};

void MeshSupportBase::set_slicer(MeshSlicerBase &slicer_base)
{
    slicer = slicer_base;

    slicer_base.get_slices(slices);

    settings = slicer_base.return_settings();
}

void MeshSupportBase::set_map(MatrixXi &hmap, MatrixXi &smap, MatrixXd &plfm) {
    height_map = hmap;
    support_map= smap;
    shrink_support_map(support_map);
    platform = plfm;
}

void MeshSupportBase::get_fermat_curves(vector<Paths> &fermat_curves)
{

    if(curve_groups.empty()) curve_group_construction();
    if(groups_fermat_curves.empty()) fermat_curves_construction();

    fermat_curves.resize(slicer.number_layer());
    for(int id = 1; id < slicer.number_layer(); id++)
    {
        Paths curves;
        get_fermat_curves(id, curves);
        fermat_curves[id] = curves;
    }

    return;
}

void MeshSupportBase::get_fermat_curves(int layer, Paths &fermat_curves) {

    for(int id = 0; id < groups_fermat_curves.size(); id++)
    {
            double layer_height = slicer.layer_height(layer);
            double max_layer_height = curve_groups[id].highest_layer_height;
            double min_layer_height = curve_groups[id].lowest_layer_height;

            if(min_layer_height < layer_height && layer_height <= max_layer_height)
            {
                Paths curves = groups_fermat_curves[id];
                fermat_curves_difference_from_model(layer, curves);
                if(curves.empty())
                {
                    curve_groups[id].highest_layer_height = layer_height;
                } else{
                    fermat_curves.insert(fermat_curves.end(), curves.begin(), curves.end());
                    groups_fermat_curves[id] = curves;
                }

            }
    }

    return;
}

void MeshSupportBase::regroup_support_points(vector<VecPointsCluster> &pin_clusters)
{
    pin_clusters.resize(settings.pillar_column * settings.pillar_row);

    for(int ir = 0; ir < settings.pillar_row; ir++)
    {
        for(int ic = 0; ic < settings.pillar_row; ic++)
        {
            regroup_support_points(ir, ic, pin_clusters);
        }
    }

    return;
}

void MeshSupportBase::regroup_support_points(int ir, int ic, vector<VecPointsCluster> &pin_clusters)
{
    int r1 = ir * settings.xy_sample_num_each_pin;
    int r2 = (ir + 1) * settings.xy_sample_num_each_pin - 1;
    int c1 = ic * settings.xy_sample_num_each_pin;
    int c2 = (ic + 1) * settings.xy_sample_num_each_pin - 1;

    MatrixXi visited = MatrixXi::Zero(r2 - r1 + 1, c2 - c1 + 1);
    MatrixXi group = MatrixXi::Zero(r2 - r1 + 1, c2 - c1 + 1);
    VecPointsCluster clusters;
    int group_index = 0;


    for(int id = r1; id <= r2; id++)
    {
        for(int jd = c1; jd <= c2; jd++)
        {
            if(!visited(id - r1, jd - c1) && support_map(id, jd) > 0)
            {
                Queue Q;
                Q.push(RowVector2i(id, jd));
                group_index ++;
                clusters.push_back(PointsCluster());
                visited(id - r1, jd - c1) = 1;
                while(!Q.empty())
                {
                    RowVector2i  u = Q.front();
                    group(u(0) - r1, u(1) - c1) = group_index;
                    IntPoint point(settings.pin_center_x(u(1)), settings.pin_center_y(u(0)));
                    clusters.back().push_back(point);
                    Q.pop();

                    for(int dr = -settings.group_expand_size; dr <= settings.group_expand_size; dr++)
                    {
                        for(int dc = -settings.group_expand_size; dc <= settings.group_expand_size; dc++)
                        {
                            int vr = u(0) + dr;
                            int vc = u(1) + dc;
                            if(r1 <= vr && vr <= r2 && c1 <= vc && vc <= c2)
                            {
                                if(!visited(vr - r1, vc - c1) && support_map(vr, vc) > 0)
                                {
                                    Q.push(RowVector2i(vr, vc));
                                    visited(vr - r1, vc - c1) = 1;
                                }
                            }
                        }
                    }

                }

            }

        }
    }

    pin_clusters[ir * settings.pillar_column + ic] = clusters;
    return;
}

int MeshSupportBase::convex_hull(PointsCluster &cluster, Path &polygon)
{
    if(cluster.empty()) return 0;

    polygon.clear();
    std::vector<ConvexHullPoint> pointslist;
    for(int id = 0; id < cluster.size(); id++)
    {
        ConvexHullPoint cpoint(cluster[id]);
        pointslist.push_back(cpoint);
        polygon.push_back(cluster[id]);
    }

    if(pointslist.size() <= 2 ) {
        return pointslist.size();
    }

    std::sort(pointslist.begin(), pointslist.end(), ConvexHullPoint_Yi_Increase);
    for(int id = 1; id < pointslist.size(); id++)
        pointslist[id].pi = minus(pointslist[id].pi, pointslist[0].pi);
    std::sort(pointslist.begin() + 1, pointslist.end(), ConvexHullPoint_TurnLeft);
    for(int id = 1; id < pointslist.size(); id++)
        pointslist[id].pi = plus(pointslist[id].pi , pointslist[0].pi);

    std::vector<int> stack;
    stack.push_back(0);

    for(int id = 1; id < pointslist.size(); id++)
    {
        IntPoint p0, p1, p2 = pointslist[id].pi;
        int p1_id = 0;
        INT64 cross_value = 0;
        do
        {
            p1 = pointslist[stack.back()].pi;
            p1_id = stack.back();
            stack.pop_back();
            if(!stack.empty())
            {
                p0 = pointslist[stack.back()].pi;
                cross_value = cross(minus(p1,p0), minus(p2,p1));
            } else
            {
                break;
            }
        }while(cross_value <= 0);
        stack.push_back(p1_id);
        stack.push_back(id);
    }

    polygon.clear();
    for(int id = 0; id < stack.size(); id++)
        polygon.push_back(pointslist[stack[id]].pi);

    Paths simple_polygon;
    ClipperLib::SimplifyPolygon(polygon, simple_polygon);
    polygon.clear();
    for(int id = 0; id < simple_polygon.size(); id++)
    {
        if(polygon.size() < simple_polygon[id].size())
            polygon = simple_polygon[id];
    }

    return polygon.size();
}

void MeshSupportBase::fermat_curves_construction() {

    groups_fermat_curves.resize(curve_groups.size());
    std::sort(curve_groups.begin(), curve_groups.end(), CurveGroupIncrease);

    IntPoint sta(-settings.mm2int(1000), -settings.mm2int(1000));
    Path fermat;
    for(int id = 0; id < curve_groups.size(); id++)
    {
        if(curve_groups[id].is_line)
        {
            fermat.clear();
            fermat = curve_groups[id].polygon;
            fermat.push_back(curve_groups[id].polygon[0]);
            groups_fermat_curves[id].push_back(fermat);
            sta = curve_groups[id].polygon[0];
        }
        else
        {
            IntPoint new_sta, end;
            if(id < curve_groups.size() - 1) {
                curve_groups[id].closet_point_polygon(end, new_sta, curve_groups[id + 1]);
            } else{
                end = IntPoint(settings.mm2int(1000), settings.mm2int(1000));
            }

            fermat.clear();
            Path polygon = curve_groups[id].polygon;
            FermatCurve fermat_curve(settings);
            fermat_curve.set_boundary(curve_groups[id].polygon, sta, end);
            fermat_curve.output_curve(fermat);
            groups_fermat_curves[id].push_back(fermat);

            //!! sometimes only after determine the fermat spiral then we figure this is not a spiral but a line
            if(!fermat_curve.is_line())
                sta = new_sta;
            else
                sta = fermat.back();
        }
    }
    return;
}

void MeshSupportBase::fermat_curves_difference_from_model(int layer, Paths &paths) {

    Clipper clipper;
    clipper.AddPaths(paths, ptSubject, false);
    clipper.AddPaths(slices[layer], ptClip, true);
    ClipperLib::PolyTree polytree;
    clipper.Execute(ctDifference, polytree, pftPositive, pftPositive);
    Paths open;  ClipperLib::OpenPathsFromPolyTree(polytree, open);
    Paths close; ClipperLib::ClosedPathsFromPolyTree(polytree, close);
    for(int id = 0; id < close.size(); id++) close[id].push_back(close[id].front());
    paths = open;
    paths.insert(paths.end(), close.begin(), close.end());
    return;
}

void MeshSupportBase::maximum_pin_layer_height_construction(MatrixXd &maximum_height) {

    maximum_height = MatrixXd::Zero(settings.pillar_row, settings.pillar_column);
    for (int ir = 0; ir < settings.pillar_row; ir++) {
        for (int ic = 0; ic < settings.pillar_column; ic++) {

            int r1 = ir * settings.xy_sample_num_each_pin;
            int r2 = (ir + 1) * settings.xy_sample_num_each_pin - 1;
            int c1 = ic * settings.xy_sample_num_each_pin;
            int c2 = (ic + 1) * settings.xy_sample_num_each_pin - 1;

            int maximum = 0;
            for(int id = r1; id <= r2; id++)
            {
                for(int jd = c1; jd <= c2; jd++)
                {
                    if(support_map(id, jd) && maximum < height_map(id, jd))
                        maximum = height_map(id, jd);
                }
            }

            maximum_height(ir, ic) = slicer.layer_height(maximum);
        }
    }
}

bool MeshSupportBase::is_cluster_line(PointsCluster &cluster) {
    double width, height;
    cluster_bounding_box(cluster, width, height);
    if(width < settings.fermat_cut_width * 2 || height < settings.fermat_cut_width)
        return true;
    return false;
}

void MeshSupportBase::cluster_bounding_box(PointsCluster &cluster, double &width, double &height)
{
    double max_x = -settings.MAX_DOUBLE;
    double min_x = settings.MAX_DOUBLE;
    double max_y = -settings.MAX_DOUBLE;
    double min_y = settings.MAX_DOUBLE;

    for(int id = 0; id < cluster.size(); id++)
    {
        if(max_x < cluster[id].X) max_x = cluster[id].X;
        if(min_x > cluster[id].X) min_x = cluster[id].X;
        if(max_y < cluster[id].Y) max_y = cluster[id].Y;
        if(min_y > cluster[id].Y) min_y = cluster[id].Y;
    }

    width = settings.int2mm(max_x - min_x);
    height = settings.int2mm(max_y - min_y);

    return;
}

void MeshSupportBase::singleline_curvegroup_construction(PointsCluster &cluster, Path &polygon)
{
    double max_x = -settings.MAX_DOUBLE;
    double min_x = settings.MAX_DOUBLE;
    double max_y = -settings.MAX_DOUBLE;
    double min_y = settings.MAX_DOUBLE;

    for(int id = 0; id < cluster.size(); id++)
    {
        if(max_x < cluster[id].X) max_x = cluster[id].X;
        if(min_x > cluster[id].X) min_x = cluster[id].X;
        if(max_y < cluster[id].Y) max_y = cluster[id].Y;
        if(min_y > cluster[id].Y) min_y = cluster[id].Y;
    }

    double width = settings.int2mm(max_x - min_x);
    double height = settings.int2mm(max_y - min_y);

    if(width < settings.extrusion_width)
    {
        min_x -= settings.mm2int(settings.extrusion_width / 2);
        max_x += settings.mm2int(settings.extrusion_width / 2);
    }
    if(height < settings.extrusion_width)
    {
        min_y -= settings.mm2int(settings.extrusion_width / 2);
        max_y += settings.mm2int(settings.extrusion_width / 2);
    }

    polygon.push_back(IntPoint(min_x,min_y));
    polygon.push_back(IntPoint(max_x,min_y));
    polygon.push_back(IntPoint(max_x,max_y));
    polygon.push_back(IntPoint(min_x,max_y));
    polygon.push_back(IntPoint(min_x,min_y));

    return;
}

void MeshSupportBase::curve_group_construction()
{
    MatrixXd maximum_height;
    maximum_pin_layer_height_construction(maximum_height);
    vector<VecPointsCluster> pin_clusters;
    regroup_support_points(pin_clusters);

    for(int ir = 0; ir < settings.pillar_row; ir++)
    {
        for(int ic = 0; ic < settings.pillar_column; ic++)
        {
            int pin_id = ir * settings.pillar_column + ic;
            if(!pin_clusters[pin_id].empty())
            {
                for(int id = 0; id < pin_clusters[pin_id].size(); id++)
                {
                    PointsCluster cluster = pin_clusters[pin_id][id];
                    Path polygon;
                    bool is_line = false;
                    if(is_cluster_line(cluster))
                    {
                        singleline_curvegroup_construction(cluster, polygon);
                        is_line = true;
                    }
                    else
                    {
                        convex_hull(cluster, polygon);
                    }
                    CurveGroup group(platform(ir, ic), maximum_height(ir, ic), polygon, is_line);
                    curve_groups.push_back(group);
                }
            }
        }
    }
}

void MeshSupportBase::clear() {
    platform.setZero();
    support_map.setZero();
    height_map.setZero();

    groups_fermat_curves.clear();
    curve_groups.clear();

    slicer.clear();
}

void MeshSupportBase::shrink_support_map(MatrixXi &smap)
{
    int nr = settings.pillar_row;
    int nc = settings.pillar_column;
    int K = settings.xy_sample_num_each_pin;
    int e = settings.sew_pin_width_factor;
    MatrixXi tmap = MatrixXi::Zero(K * nr, K * nc);
    for(int ir = 0; ir < settings.pillar_row; ir ++)
    {
        for(int ic = 0; ic < settings.pillar_column; ic++)
        {
            int r1 = ir * K + e;
            int r2 = (ir + 1) * K - 1 - e;
            int c1 = ic * K + e;
            int c2 = (ic + 1) * K - 1 - e;
            for(int id = r1; id < r2; id++)
            {
                for(int jd = c1; jd < c2; jd++)
                {
                    tmap(id, jd) = smap(id, jd);
                }
            }
        }
    }

    smap = tmap;
}

#endif //SUPPORTER_MESH_SUPPORT_BASE_H
