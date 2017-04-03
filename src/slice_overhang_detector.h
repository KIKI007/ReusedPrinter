//
// Created by 汪子琦 on 3/28/17.
//

#ifndef SUPPORTER_SLICE_OVERHANG_DETECTOR_H
#define SUPPORTER_SLICE_OVERHANG_DETECTOR_H

#include "slice.h"
#include <list>

typedef struct tagEDGE
{
    int xi;
    double dx;
    int ymax;
    int ymin;
}EDGE;

bool EdgeXiComparator(EDGE& e1, EDGE& e2)
{
    return (e1.xi <= e2.xi);
}

class Slice_Overhang_Detector: public Slice
{
public:
    Slice_Overhang_Detector():Slice() { }

private:
    int num_rows(int ymax, int ymin) {  return (double)(ymax - ymin) / settings.mm2int(settings.sample_distance) + 2; }
    int row_y(int row) { return row * settings.mm2int(settings.sample_distance); }
    int y_rows(int y, int ymax, int ymin);
    double anti_clockwise_angle(Eigen::Vector2d p0, Eigen::Vector2d p1);
private:
    void get_bounding_box(ClipperLib::Paths &poly, int &ymax, int &ymin);
    void init_scanline_new_edgetable(std::vector< std::list<EDGE> >& slNet, ClipperLib::Paths &poly, int ymax, int ymin);
    void scanline_polygon_fill(ClipperLib::Paths &poly, std::vector<Eigen::Vector2d> &sample_vertices);
    void polygon_edge_fill(ClipperLib::Paths &poly, std::vector<Eigen::Vector2d> &sample_vertices);
    void insert_net_into_aet(std::list<EDGE> &net, std::list<EDGE> &aet, int y);
    void fill_ate_scanline(std::list<EDGE> &aet, int y, std::vector<Eigen::Vector2d> &sample_vertices);
    void remove_non_active_edge_from_aet(std::list<EDGE> &aet, int y);
    void update_and_resort_aet(std::list<EDGE>& aet, int dy);
    void process_scanline_fill(std::vector< std::list<EDGE> >& slNet, std::vector<Eigen::Vector2d> &sample_vertices, int ymax, int ymin);
    int  get_winding_number(double x, double y, ClipperLib::Paths poly);

public:
    void removing_overlap(std::vector<ClipperLib::Paths> &slices);
    void sampling(Eigen::MatrixXd &SP);

private:
};

int Slice_Overhang_Detector::y_rows(int y, int ymax, int ymin)
{
    int dy =  settings.mm2int(settings.sample_distance);
    if((y - ymin) % dy == 0)
        return (y - ymin) / dy;
    else
        return (y - ymin) / dy + 1;
}

void Slice_Overhang_Detector::get_bounding_box(ClipperLib::Paths &poly, int &ymax, int &ymin)
{
    assert(!poly.empty());
    ymin = ymax = poly[0][0].Y;
    for(size_t id = 0; id < poly.size(); id++)
        for(size_t jd = 0; jd < poly[id].size(); jd++)
        {
            if(ymax < poly[id][jd].Y)
                ymax = poly[id][jd].Y;
            if(ymin > poly[id][jd].Y)
                ymin = poly[id][jd].Y;
        }
}

double Slice_Overhang_Detector::anti_clockwise_angle(Eigen::Vector2d p0, Eigen::Vector2d p1)
{
    if(p0.norm() < settings.ZERO_EPS || p1.norm() < settings.ZERO_EPS) return 0;
    p0 /= p0.norm();
    p1 /= p1.norm();
    double sign = p0.x() * p1.y() - p0.y() * p1.x();
    if(sign > 0)
        return std::acos(p0.dot(p1));
    else
        return -std::acos(p0.dot(p1));
}


void Slice_Overhang_Detector::init_scanline_new_edgetable(std::vector< std::list<EDGE> >& slNet, ClipperLib::Paths &poly, int ymax, int ymin)
{
    std::cout << num_rows(ymax, ymin) << std::endl;
    slNet.resize(num_rows(ymax, ymin));
    for(int id = 0; id < poly.size(); id++)
    {
        int size_poly_id = poly[id].size();
        for (int jd = 0; jd < poly[id].size(); jd++)
        {
            ClipperLib::IntPoint pe, ps, pss, pee;
            ps = poly[id][jd];
            pe = poly[id][(jd + 1) % size_poly_id];
            pee = poly[id][(jd + 2) % size_poly_id];
            pss = poly[id][(jd - 1 + size_poly_id) % size_poly_id];
            EDGE e;
            if (ps.Y != pe.Y) {
                e.dx = (double) (pe.X - ps.X) / (double) (pe.Y - ps.Y);
                if (pe.Y > ps.Y) {
                    e.xi = ps.X;
                    e.ymin = ps.Y;
                    if (pee.Y >= pe.Y)
                        e.ymax = pe.Y - 1;
                    else
                        e.ymax = pe.Y;
                    slNet[y_rows(ps.Y, ymax, ymin)].push_front(e);
                } else {
                    e.xi = pe.X;
                    e.ymin = pe.Y;
                    if (pss.Y >= ps.Y)
                        e.ymax = ps.Y - 1;
                    else
                        e.ymax = ps.Y;
                    slNet[y_rows(pe.Y, ymax, ymin)].push_front(e);
                }
            }
        }
    }
    return;
}

void Slice_Overhang_Detector::scanline_polygon_fill(ClipperLib::Paths &poly, std::vector<Eigen::Vector2d> &sample_vertices)
{
    assert(!poly.empty());

    int ymin, ymax;
    get_bounding_box(poly, ymax, ymin);
    std::vector< std::list<EDGE>> slNet;
    init_scanline_new_edgetable(slNet, poly, ymax, ymin);
    //polygon_edge_fill(poly, sample_vertices);
    process_scanline_fill(slNet, sample_vertices ,ymax, ymin);
}

void Slice_Overhang_Detector::polygon_edge_fill(ClipperLib::Paths &poly, std::vector<Eigen::Vector2d> &sample_vertices)
{
    int dx = settings.mm2int(settings.sample_distance);
    for(int id = 0; id < poly.size(); id++)
    {
        int size_poly_id = poly[id].size();
        for (int jd = 0; jd < poly[id].size(); jd++) {
            Eigen::Vector2d ps, pe;
            ps = Eigen::Vector2d(poly[id][jd].X,
                                 poly[id][jd].Y);
            pe = Eigen::Vector2d(poly[id][(jd + 1) % size_poly_id].X,
                                 poly[id][(jd + 1) % size_poly_id].Y);
            int num_sample = (ps - pe).norm() / (double)dx;
            for(int kd = 0; kd < num_sample; kd++)
            {
                Eigen::Vector2d sample = (ps - pe) * (double) (kd + 1) / (double)(num_sample + 1) + pe;
                sample_vertices.push_back(sample);
            }
            sample_vertices.push_back(ps);
        }
    }
}

void Slice_Overhang_Detector::insert_net_into_aet(std::list<EDGE> &net, std::list<EDGE> &aet, int y)
{
    std::list<EDGE>::iterator nit, ait;
    for(nit = net.begin(); nit != net.end(); ++nit)
    {
        nit->xi += (double)(y - nit->ymin) * nit->dx;
        for(ait = aet.begin(); ait != aet.end(); ++ait)
        {
            if(ait->xi > nit->xi)
                break;
        }
        aet.insert(ait, *nit);
    }

    return;
}

void Slice_Overhang_Detector::fill_ate_scanline(std::list<EDGE> &aet, int y, std::vector<Eigen::Vector2d> &sample_vertices)
{
    std::list<EDGE>::iterator ait, it0, it1;
    int dx = settings.mm2int(settings.sample_distance), sx;
    for(ait = aet.begin(); ait != aet.end();)
    {
        it0 = ait++; if(ait == aet.end()) break;
        it1 = ait++;
        int partition_size = std::abs((double)(it1->xi - it0->xi)) / (double) dx;
        for(int kd = 1; kd <= partition_size; kd++)
        {
            sx = (double)(it1->xi - it0->xi) * (double) kd / (double) (partition_size + 1) + it0->xi;
            sample_vertices.push_back(Eigen::Vector2d(sx, y));
        }
    }
}

void Slice_Overhang_Detector::remove_non_active_edge_from_aet(std::list<EDGE> &aet, int y)
{
    std::list<EDGE>::iterator ait;
    for(ait = aet.begin(); ait != aet.end();)
    {
        if(ait->ymax < y)
        {
            ait = aet.erase(ait);
        }
        else
        {
            ait++;
        }
    }
    return;
}

void Slice_Overhang_Detector::update_and_resort_aet(std::list<EDGE>& aet, int dy)
{
    std::list<EDGE>::iterator ait;
    for(ait = aet.begin(); ait != aet.end(); ait++)
    {
        ait->xi += ait->dx * dy;
    }

    aet.sort(EdgeXiComparator);
}

void Slice_Overhang_Detector::process_scanline_fill(std::vector< std::list<EDGE> >& slNet, std::vector<Eigen::Vector2d> &sample_vertices, int ymax, int ymin)
{
    std::list<EDGE> aet;
    int dy = settings.mm2int(settings.sample_distance);
    for(int y = ymin; y < ymax; y += dy)
    {
        std::cout << "y, " << y << std::endl;
        std::cout << "insert" << std::endl;
        insert_net_into_aet(slNet[y_rows(y, ymax, ymin)], aet, y);
        remove_non_active_edge_from_aet(aet, y);
        std::cout << "fill" << std::endl;
        if(y != ymin)
            fill_ate_scanline(aet, y, sample_vertices);
        std::cout << "remove" << std::endl;
        std::cout << "resort" << std::endl;
        update_and_resort_aet(aet, dy);
    }
}

int Slice_Overhang_Detector::get_winding_number(double x, double y, ClipperLib::Paths poly)
{
    Eigen::Vector3d center(x, y, 0);
    double wing_angle = 0;
    for(int id = 0; id < poly.size(); id++)
    {
        for(int jd = 0; jd < poly[id].size(); jd++)
        {
            Eigen::Vector3d p0(settings.int2mm(poly[id][jd].X),
                               settings.int2mm(poly[id][jd].Y), 0);
            Eigen::Vector3d p1(settings.int2mm(poly[id][(jd + 1) % poly[id].size()].X),
                               settings.int2mm(poly[id][(jd + 1) % poly[id].size()].Y), 0);
            if((center - p0).norm() < settings.ZERO_EPS || (center - p1).norm() < settings.ZERO_EPS)
                return 1;

            if((center - p0).cross(center - p1)[2] < settings.ZERO_EPS)
            {
                if((center - p0).dot(center - p1) < 0)
                    return 1;
            }
        }
    }
//    std::cout << std::round(wing_angle / (2 * settings.PI)) << std::endl;
//    return std::round(wing_angle / (2 * settings.PI));
    return 0;
}

void Slice_Overhang_Detector::sampling(Eigen::MatrixXd &SP)
{
    std::vector<Eigen::RowVector3d> sample_vertices;

    //inner points which require support
    std::vector<ClipperLib::Paths> overlap_removed_slices;
    removing_overlap(overlap_removed_slices);
    for(int id = 0; id < overlap_removed_slices.size(); id++)
    {
        ClipperLib::Paths poly = overlap_removed_slices[id];
        if(!poly.empty())
        {
            std::vector<Eigen::Vector2d> vlist;
            scanline_polygon_fill(poly, vlist);
            for(int kd = 0; kd < vlist.size(); kd++)
                sample_vertices.push_back(Eigen::RowVector3d(vlist[kd][0] * settings.UNIT,settings.int2mm(P[id]), vlist[kd][1] * settings.UNIT));
        }
    }

    //outer points(perimeter) which require support
    std::vector<ClipperLib::Paths> intensive_slices;
    intensive_slices.resize(layer_slices.size());
    for(int id = 0; id < layer_slices.size(); id++)
    {
        intensive_slices[id].resize(layer_slices[id].size());
        for(int jd = 0; jd < layer_slices[id].size(); jd++)
        {
            int poly_size = layer_slices[id][jd].size();
            for(int kd = 0; kd < layer_slices[id][jd].size(); kd++)
            {
                ClipperLib::IntPoint ps, pe;
                ps = layer_slices[id][jd][kd];
                pe = layer_slices[id][jd][(kd + 1) % poly_size];

                double dist_x = settings.int2mm(ps.X - pe.X);
                double dist_y = settings.int2mm(ps.Y - pe.Y);
                double dist = std::sqrt(dist_x * dist_x + dist_y * dist_y);
                int sequences = dist / settings.sample_distance + 2;
                for(int ld = 0; ld < sequences - 1; ld++)
                {
                    intensive_slices[id][jd].push_back(ClipperLib::IntPoint(
                            (pe.X - ps.X) * (double) ld / (sequences - 1) + ps.X,
                            (pe.Y - ps.Y) * (double) ld / (sequences - 1) + ps.Y
                    ));
                }
            }
        }
    }

    for(int id = 1; id < intensive_slices.size(); id++)
    {
        if(overlap_removed_slices[id].size() != 0)
        {
            ClipperLib::Paths poly = overlap_removed_slices[id];
            for(int jd = 0; jd < intensive_slices[id].size(); jd++)
            {
                std::vector<bool> supported;
                int poly_size = intensive_slices[id][jd].size();
                supported.resize(poly_size, false);

                for (int kd = 0; kd < intensive_slices[id][jd].size(); kd++)
                {
                    ClipperLib::IntPoint p0 = intensive_slices[id][jd][kd];
                    int winding_number = get_winding_number(
                            settings.int2mm(p0.X),
                            settings.int2mm(p0.Y),
                            poly
                    );
                    if(winding_number == 0)
                        supported[kd] = true;
                }

                for (int kd = 0; kd < intensive_slices[id][jd].size(); kd++)
                {
                    ClipperLib::IntPoint p0 = intensive_slices[id][jd][kd];
                    if(supported[kd] == false)
                    {
                        if(supported[(kd + 1) % poly_size] == false
                           && supported[(kd - 1 + poly_size) % poly_size] == false)
                        {
                            sample_vertices.push_back(Eigen::Vector3d(settings.int2mm(p0.X), settings.int2mm(P[id]), settings.int2mm(p0.Y)));
                        }
                    }
                }
            }
        }
    }

    SP.resize(sample_vertices.size(), 3);
    for(int id = 0; id < sample_vertices.size(); id++)
        SP.row(id) = sample_vertices[id];
    return;
}

void Slice_Overhang_Detector::removing_overlap(std::vector<ClipperLib::Paths> &slices)
{
    assert(!V.isZero());

    if(layer_slices.empty())
        contour_construction();

    slices.clear();
    slices.resize(layer_slices.size());

    for(int id = 1; id < layer_slices.size(); id++)
    {
        ClipperLib::Clipper clipper;
        ClipperLib::ClipperOffset co;
        ClipperLib::Paths upper = layer_slices[id], lower = layer_slices[id - 1];

        co.AddPaths(lower, ClipperLib::jtRound, ClipperLib::etClosedPolygon);
        co.Execute(lower, settings.mm2int(settings.overhang_offset));
        ClipperLib::SimplifyPolygons(lower);

        clipper.StrictlySimple(true);
        clipper.AddPaths(upper, ClipperLib::ptSubject, true);
        clipper.AddPaths(lower, ClipperLib::ptClip, true);
        clipper.Execute(ClipperLib::ctDifference, slices[id], ClipperLib::pftPositive, ClipperLib::pftPositive);

        ClipperLib::SimplifyPolygons(slices[id]);
    }

    slices[0].clear();
    return;
}

#endif //SUPPORTER_SLICE_OVERHANG_DETECTOR_H
