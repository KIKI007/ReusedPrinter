//
// Created by 汪子琦 on 4/27/17.
//

#ifndef SUPPORTER_SCAN_LINE_FILL_H
#define SUPPORTER_SCAN_LINE_FILL_H

#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <iostream>
#include <igl/png/writePNG.h>
#include <igl/png/readPNG.h>

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

class ScanLineFill
{

public:

    ScanLineFill(bool expand_fill = false)
    {
        N = settings.xy_sample_num_each_pin;
        sample_width = settings.sample_width;
        expand = expand_fill;
    }

public:

    void polygon_fill(ClipperLib::Paths &poly, Eigen::MatrixXi &imag);

protected:

    void bounding_box(ClipperLib::Paths &poly);

    void init_scanline_new_edgetable(std::vector< std::list<EDGE>>& slNet, ClipperLib::Paths &poly);

    void process_scanline_fill(std::vector< std::list<EDGE> >& slNet,  Eigen::MatrixXi &image);

    void insert_net_into_aet(std::list<EDGE> &net, std::list<EDGE> &aet, int y);

    void remove_non_active_edge_from_aet(std::list<EDGE> &aet, int y);

    void fill_ate_scanline(std::list<EDGE> aet, int y, Eigen::MatrixXi &image);

    void update_and_resort_aet(std::list<EDGE>& aet, int dy);

    void process_horizontal_line(ClipperLib::Paths poly, Eigen::MatrixXi &imag);

protected:

    //return the layer before y
    //if y just in that layer, return that layer
    int before_layer(int y) { y-= sample_width / 2; return std::floor((double)y / sample_width); }

    //return the layer after y
    //if y just in that layer, return that layer
    int after_layer(int y)  { y-= sample_width / 2; return std::ceil((double)y / sample_width); }

    int layer(int y){ assert(on_layer(y)); return (y - sample_width / 2) / sample_width; }

    bool on_layer(int y) {
        y-= sample_width / 2;
        if(y % sample_width == 0)
            return true;
        else
            return false;
    }

    int layer_height(int layer)
    {
        return layer * sample_width + sample_width / 2;
    }

private:
    Settings settings;
    bool expand;

private:
    int sta_layer;
    int end_layer;
    int N;
    int sample_width;
};

void ScanLineFill::polygon_fill(ClipperLib::Paths &poly, Eigen::MatrixXi &imag)
{
    imag = Eigen::MatrixXi::Zero(settings.pillar_row * N, settings.pillar_column * N);
    if(poly.empty()) return;
    bounding_box(poly);
    std::vector< std::list<EDGE>> slNet;
    process_horizontal_line(poly, imag);
    init_scanline_new_edgetable(slNet, poly);
    process_scanline_fill(slNet, imag);
}

void ScanLineFill::bounding_box(ClipperLib::Paths &poly)
{
    assert(!poly.empty());
    int ymin = poly[0][0].Y, ymax = poly[0][0].Y;
    for(size_t id = 0; id < poly.size(); id++)
        for(size_t jd = 0; jd < poly[id].size(); jd++)
        {
            if(ymax < poly[id][jd].Y)
                ymax = poly[id][jd].Y;
            if(ymin > poly[id][jd].Y)
                ymin = poly[id][jd].Y;
        }
    sta_layer = before_layer(ymin);
    end_layer = after_layer(ymax);
}

void ScanLineFill::init_scanline_new_edgetable(std::vector< std::list<EDGE> >& slNet, ClipperLib::Paths &poly)
{
    slNet.resize(end_layer - sta_layer + 1);
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
                    slNet[after_layer(ps.Y) - sta_layer].push_front(e);
                } else {
                    e.xi = pe.X;
                    e.ymin = pe.Y;
                    if (pss.Y >= ps.Y)
                        e.ymax = ps.Y - 1;
                    else
                        e.ymax = ps.Y;
                    slNet[after_layer(pe.Y) - sta_layer].push_front(e);
                }
            }
        }
    }
    return;
}

void ScanLineFill::process_horizontal_line(ClipperLib::Paths poly, Eigen::MatrixXi &imag)
{
    for(int id = 0; id < poly.size(); id++)
    {
        int size_poly_id = poly[id].size();
        for (int jd = 0; jd < poly[id].size(); jd++) {
            ClipperLib::IntPoint pe, ps;
            ps = poly[id][jd];
            pe = poly[id][(jd + 1) % size_poly_id];
            if (ps.Y == pe.Y && on_layer(pe.Y))
            {
                int X0, X1, row = before_layer(pe.Y), c0, c1;
                if(ps.X < pe.X){X0 = ps.X, X1 = pe.X;}
                else           {X1 = ps.X, X0 = pe.X;}

                if(expand) c0 = before_layer(X0), c1 = after_layer(X1);
                else c0 = after_layer(X0), c1 = before_layer(X1);

                for(int col = c0; col <= c1; col++)
                {
                    imag(row, col) = 1;
                }
            }
        }
    }
    return;
}


void ScanLineFill::process_scanline_fill(std::vector< std::list<EDGE> >& slNet,  Eigen::MatrixXi &image)
{
    std::list<EDGE> aet;
    for(int layer = sta_layer; layer < end_layer; layer += 1)
    {
        int y = layer_height(layer);
        insert_net_into_aet(slNet[layer - sta_layer], aet, y);
        remove_non_active_edge_from_aet(aet, y);
        fill_ate_scanline(aet, y, image);
        update_and_resort_aet(aet, sample_width);
    }
}

void ScanLineFill::insert_net_into_aet(std::list<EDGE> &net, std::list<EDGE> &aet, int y)
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

void ScanLineFill::remove_non_active_edge_from_aet(std::list<EDGE> &aet, int y)
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

void ScanLineFill::fill_ate_scanline(std::list<EDGE> aet, int y, Eigen::MatrixXi &imag)
{
    std::list<EDGE>::iterator ait, it0, it1;
    int X0 = 0, X1 = 0;
    int c0, c1;
    int row = layer(y);
    for(ait = aet.begin(); ait != aet.end();)
    {
        it0 = ait++; if(ait == aet.end()) break;
        it1 = ait++;
        if(it0->xi < it1->xi) {X0 = it0->xi, X1 = it1->xi;}
        else {X1 = it0->xi, X0 = it1->xi;}

        if(expand) c0 = before_layer(X0), c1 = after_layer(X1);
        else c0 = after_layer(X0), c1 = before_layer(X1);

        for(int col = c0; col <= c1; col++)
        {
            imag(row, col) = 1;
        }
    }
}

void ScanLineFill::update_and_resort_aet(std::list<EDGE>& aet, int dy)
{
    std::list<EDGE>::iterator ait;
    for(ait = aet.begin(); ait != aet.end(); ait++)
    {
        ait->xi += ait->dx * dy;
    }

    aet.sort(EdgeXiComparator);
}

#endif