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

public:

    void get_bounding_box(ClipperLib::Paths &poly, int &ymax, int &ymin)
    {
        assert(!poly.empty());
        ymin = ymax = poly[0][0].Y;
        for(size_t id = 0; id < poly.size(); id++)
            for(size_t jd = 0; jd < poly[id].size(); jd++)
            {
                if(ymax > poly[id][jd].Y)
                    ymax = poly[id][jd].Y;
                if(ymin < poly[id][jd].Y)
                    ymin = poly[id][jd].Y;
            }
    }

    void init_scanline_new_edgetable(std::vector< std::list<EDGE> >& slNet, ClipperLib::Paths &poly, int ymax, int ymin)
    {
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
                        if (pee.Y >= pe.Y)
                            e.ymax = pe.Y - 1;
                        else
                            e.ymax = pe.Y;
                        slNet[ps.Y - ymin].push_front(e);
                    } else {
                        e.xi = pe.X;
                        if (pss.Y >= ps.Y)
                            e.ymax = ps.Y - 1;
                        else
                            e.ymax = ps.Y;
                        slNet[pe.Y - ymin].push_front(e);
                    }
                }
            }
        }
        return;
    }

    void ScanLinePolygonFill(ClipperLib::Paths &poly, std::vector<Eigen::Vector2d> &sample_vertices)
    {
        assert(!poly.empty());

        int ymin, ymax;
        get_bounding_box(poly, ymax, ymin);
        std::vector< std::list<EDGE>> slNet;
        init_scanline_new_edgetable(slNet, poly, ymax, ymin);
        horizon_edge_fill(poly, sample_vertices);
        process_scanline_fill(slNet, sample_vertices ,ymax, ymin);
    }

    void horizon_edge_fill(ClipperLib::Paths &poly, std::vector<Eigen::Vector2d> &sample_vertices)
    {
        int dx = settings.mm2int(settings.sample_distance);
        for(int id = 0; id < poly.size(); id++)
        {
            int size_poly_id = poly[id].size();
            for (int jd = 0; jd < poly[id].size(); jd++) {
                ClipperLib::IntPoint pe, ps, pss, pee;
                ps = poly[id][jd];
                pe = poly[id][(jd + 1) % size_poly_id];
                if (ps.Y == pe.Y)
                {
                    int sx;
                    int partition_size = std::abs((double)(pe.X - ps.X)) / (double) dx + 1;
                    for(int kd = 0; kd <= partition_size; kd++)
                    {
                        sx = (double)(pe.X - ps.X) * kd / partition_size + ps.X;
                        sample_vertices.push_back(Eigen::Vector2d(sx, ps.Y));
                    }
                }
            }
        }
    }

    void insert_net_into_aet(std::list<EDGE> &net, std::list<EDGE> &aet)
    {
        std::list<EDGE>::iterator nit, ait;
        for(nit = net.begin(); nit != net.end(); ++nit)
        {
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
            it0 = ait++;
            it1 = ait++;
            int partition_size = std::abs((double)(it1->xi - it0->xi)) / (double) dx + 1;
            for(int kd = 0; kd <= partition_size; kd++)
            {
                sx = (double)(it1->xi - it0->xi) * kd / partition_size + it0->xi;
                sample_vertices.push_back(Eigen::Vector2d(sx, y));
            }
        }
    }

    void remove_non_active_edge_from_aet(std::list<EDGE> &aet, int y)
    {
        std::list<EDGE>::iterator ait;
        for(ait = aet.begin(); ait != aet.end();)
        {
            if(ait->ymax == y)
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
    void update_and_resort_aet(std::list<EDGE>& aet)
    {
        std::list<EDGE>::iterator ait;
        for(ait = aet.begin(); ait != aet.end();)
        {
            ait->xi += ait->dx;
        }

        aet.sort(EdgeXiComparator);
    }

    void process_scanline_fill(std::vector< std::list<EDGE> >& slNet, std::vector<Eigen::Vector2d> &sample_vertices, int ymax, int ymin)
    {
        std::list<EDGE> aet;
        for(int y = ymin; y <= ymax; y++)
        {
            insert_net_into_aet(slNet[y - ymin], aet);
            fill_ate_scanline(aet, y, sample_vertices);
            remove_non_active_edge_from_aet(aet, y);
            update_and_resort_aet(aet);
        }
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
            co.Execute(lower, 100);

            clipper.StrictlySimple(true);
            clipper.AddPaths(upper, ClipperLib::ptSubject, true);
            clipper.AddPaths(lower, ClipperLib::ptClip, true);
            clipper.Execute(ClipperLib::ctDifference, solutions[id], ClipperLib::pftNonZero, ClipperLib::pftNonZero);

            ClipperLib::SimplifyPolygons(solutions[id]);
        }

        layer_slices = solutions;
        return;
    }

private:
};

#endif //SUPPORTER_SLICE_OVERHANG_DETECTOR_H
