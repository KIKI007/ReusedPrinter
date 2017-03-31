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



public:

    Slice_Overhang_Detector():Slice()
    {

    }

private:

    int num_rows(int ymax, int ymin) {  return (double)(ymax - ymin) / settings.mm2int(settings.sample_distance) + 2; }

    int row_y(int row) { return row * settings.mm2int(settings.sample_distance); }

    int y_rows(int y, int ymax, int ymin)
    {
        int dy =  settings.mm2int(settings.sample_distance);
        if((y - ymin) % dy == 0)
            return (y - ymin) / dy;
        else
            return (y - ymin) / dy + 1;
    }

public:

    void get_bounding_box(ClipperLib::Paths &poly, int &ymax, int &ymin)
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

    void init_scanline_new_edgetable(std::vector< std::list<EDGE> >& slNet, ClipperLib::Paths &poly, int ymax, int ymin)
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

    void scanline_polygon_fill(ClipperLib::Paths &poly, std::vector<Eigen::Vector2d> &sample_vertices)
    {
        assert(!poly.empty());

        int ymin, ymax;
        get_bounding_box(poly, ymax, ymin);
        std::vector< std::list<EDGE>> slNet;
        init_scanline_new_edgetable(slNet, poly, ymax, ymin);
        polygon_edge_fill(poly, sample_vertices);
        process_scanline_fill(slNet, sample_vertices ,ymax, ymin);
    }

    void polygon_edge_fill(ClipperLib::Paths &poly, std::vector<Eigen::Vector2d> &sample_vertices)
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

    void insert_net_into_aet(std::list<EDGE> &net, std::list<EDGE> &aet, int y)
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

    void fill_ate_scanline(std::list<EDGE> &aet, int y, std::vector<Eigen::Vector2d> &sample_vertices)
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

    void remove_non_active_edge_from_aet(std::list<EDGE> &aet, int y)
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
    void update_and_resort_aet(std::list<EDGE>& aet, int dy)
    {
        std::list<EDGE>::iterator ait;
        for(ait = aet.begin(); ait != aet.end(); ait++)
        {
            ait->xi += ait->dx * dy;
        }

        aet.sort(EdgeXiComparator);
    }

    void process_scanline_fill(std::vector< std::list<EDGE> >& slNet, std::vector<Eigen::Vector2d> &sample_vertices, int ymax, int ymin)
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

    void sampling(Eigen::MatrixXd &SP)
    {
        removing_overlap();
        std::vector<Eigen::RowVector3d> sample_vertices;
        for(int id = 0; id < number_layer(); id++)
        {
            ClipperLib::Paths poly = layer_slices[id];
            if(!poly.empty())
            {
                std::vector<Eigen::Vector2d> vlist;
                scanline_polygon_fill(poly, vlist);
                for(int kd = 0; kd < vlist.size(); kd++)
                    sample_vertices.push_back(Eigen::RowVector3d(vlist[kd][0] * settings.UNIT,settings.int2mm(P[id]), vlist[kd][1] * settings.UNIT));
            }
        }

        SP.resize(sample_vertices.size(), 3);
        for(int id = 0; id < sample_vertices.size(); id++)
            SP.row(id) = sample_vertices[id];
        return;
    }

    void removing_overlap()
    {
        assert(!V.isZero());

        if(layer_slices.empty())
            contour_construction();

        std::vector<ClipperLib::Paths> solutions;
        solutions.resize(layer_slices.size());

        for(int id = 1; id < layer_slices.size(); id++)
        {
            ClipperLib::Clipper clipper;
            ClipperLib::ClipperOffset co;
            ClipperLib::Paths upper = layer_slices[id], lower = layer_slices[id - 1];

            co.AddPaths(lower, ClipperLib::jtRound, ClipperLib::etClosedPolygon);
            co.Execute(lower, 400);

            clipper.StrictlySimple(true);
            clipper.AddPaths(upper, ClipperLib::ptSubject, true);
            clipper.AddPaths(lower, ClipperLib::ptClip, true);
            clipper.Execute(ClipperLib::ctDifference, solutions[id], ClipperLib::pftPositive, ClipperLib::pftPositive);

            ClipperLib::SimplifyPolygons(solutions[id]);
        }
        for(int id = 1; id < layer_slices.size(); id++)
            layer_slices[id] = solutions[id];
        return;
    }

private:
};

#endif //SUPPORTER_SLICE_OVERHANG_DETECTOR_H
