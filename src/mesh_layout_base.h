//
// Created by 汪子琦 on 5/2/17.
//

#ifndef SUPPORTER_MESH_LAYOUT_BASE_H
#define SUPPORTER_MESH_LAYOUT_BASE_H

#include "mesh_slicer_overhang.h"
#include "scan_line_fill.h"
#include "scene_organizer.h"

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::Vector2d;

typedef std::vector< Eigen::MatrixXi, Eigen::aligned_allocator<Eigen::MatrixXi>> VecMatrixXi;
typedef std::vector< Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd>> VecMatrixXd;

class MeshLayoutBase
{

public:
    MeshLayoutBase()
    {

    }

public: /* Input & Utility */

    //set slicer for layout optimization
    void set_slicer(MeshSlicerOverhang &slicer_overhang);

    void clear();

    //move all map a vector of (x, z)
    void transform_map_xz(double x, double z);

    //rotate all map in y axis with the center and the angle
    void rotate_map_yaxis(double angle, Vector2d &center);

public: /* Output */

    void get_height_map(MatrixXi &hmap);

    void get_height_map(MatrixXd &hmap);

    void get_support_map(MatrixXi &smap);

    void get_platform(MatrixXi &hmap, MatrixXi &smap, MatrixXd &platform);

    void get_volume_support(MatrixXi &hmap, MatrixXi &smap, MatrixXd &platform, MatrixXd &V, MatrixXi &F);

protected:

    void height_map_construction();

    void support_map_construction();

    void transform_map_xz(MatrixXi &map, double x, double y, int fill);

    void rotate_map_yaxis(MatrixXi &map, double angle, Vector2d &center, int fill);

protected:
    MatrixXi height_map;   //value is the layer id

    MatrixXi support_map; // value is 0 or 1(repesent this matrix's element is a support point)

    MeshSlicerOverhang slicer;

    Settings settings;
};

void MeshLayoutBase::set_slicer(MeshSlicerOverhang &slicer_overhang) {
    slicer = slicer_overhang;
    settings = slicer.return_settings();
}

void MeshLayoutBase::clear() {
    height_map.setZero();
    support_map.setZero();
    slicer.clear();
}

void MeshLayoutBase::transform_map_xz(double x, double z) {

    if(!height_map.isZero() && !support_map.isZero())
    {
        transform_map_xz(height_map, x, z, slicer.number_layer());
        transform_map_xz(height_map, x, z, slicer.number_layer());
    }
}

void MeshLayoutBase::rotate_map_yaxis(double angle, Vector2d &center)
{
    if(!height_map.isZero() && !support_map.isZero())
    {
        rotate_map_yaxis(height_map, angle, center, slicer.number_layer());
        rotate_map_yaxis(height_map, angle, center, slicer.number_layer());
    }
}

void MeshLayoutBase::get_height_map(MatrixXi &hmap) {
    if(height_map.isZero())
        height_map_construction();
    hmap = height_map;
}

void MeshLayoutBase::get_height_map(MatrixXd &hmap) {
    if(height_map.isZero())
        height_map_construction();

    int nr = height_map.rows();
    int nc = height_map.cols();

    hmap.resize(nr, nc);

    for(int ir = 0; ir < nr; ir++)
    {
        for(int ic = 0; ic < nc; ic++)
        {
            if(height_map(ir, ic)  < slicer.number_layer())
                hmap(ir, ic) = slicer.layer_height(height_map(ir, ic));
            else
                hmap(ir, ic) = settings.maximum_height_map;
        }
    }

    return;
}

void MeshLayoutBase::get_support_map(MatrixXi &smap) {
    if(support_map.isZero())
        support_map_construction();
    smap = height_map;
}

void MeshLayoutBase::height_map_construction()
{
    int nc = settings.pillar_column  * settings.xy_sample_num_each_pin;
    int nr = settings.pillar_row  * settings.xy_sample_num_each_pin;
    height_map = MatrixXi::Ones(nr, nc) * slicer.number_layer();

    std::vector<Paths> slices;
    slicer.get_slices(slices);

    VecMatrixXi slices_map;
    slices_map.resize(slicer.number_layer());

    for(int layer = 0; layer < slicer.number_layer(); layer++)
    {
        ScanLineFill fill(false);
        fill.polygon_fill(slices[layer], slices_map[layer]);
    }

    for(int layer = 0; layer < slicer.number_layer(); layer++)
    {
        for(int ir = 0; ir < nr; ir++)
        {
            for(int ic = 0; ic < nc; ic++)
            {
                if(slices_map[layer](ir, ic) && height_map(ir, ic) > layer)
                {
                    height_map(ir, ic) = layer;
                }
            }
        }
    }
}

void MeshLayoutBase::support_map_construction()
{
    int nc = settings.pillar_column  * settings.xy_sample_num_each_pin;
    int nr = settings.pillar_row  * settings.xy_sample_num_each_pin;
    support_map = MatrixXi::Zero(nr, nc);

    std::vector<Paths> bottom_half;
    slicer.get_bottom_half_overhang(bottom_half);

    ClipperLib::Paths downside;
    for(int layer = 1; layer < slicer.number_layer(); layer++)
    {
        ClipperLib::Clipper clipper;
        clipper.AddPaths(downside, ClipperLib::ptSubject, true);
        clipper.AddPaths(bottom_half[layer], ClipperLib::ptClip, true);
        clipper.Execute(ClipperLib::ctUnion, downside, ClipperLib::pftPositive, ClipperLib::pftPositive);
    }

    ScanLineFill fill(false);
    fill.polygon_fill(downside, support_map);
}

void MeshLayoutBase::get_platform(MatrixXi &hmap, MatrixXi &smap, MatrixXd &platform)
{
    platform = MatrixXd::Zero(9, 11);

    for(int ir = 0; ir < settings.pillar_row; ir++)
    {
        for(int ic = 0; ic < settings.pillar_column; ic++)
        {
            int r1 = ir * settings.xy_sample_num_each_pin  ;
            int r2 = (ir + 1) * settings.xy_sample_num_each_pin - 1;
            int c1 = ic * settings.xy_sample_num_each_pin ;
            int c2 = (ic + 1) * settings.xy_sample_num_each_pin - 1;

            int num_red = 0;
            int minimum_layer = slicer.number_layer();
            for(int id = r1; id <= r2; id++)
            {
                for(int jd = c1; jd <= c2; jd++)
                {
                    if(minimum_layer > height_map(id, jd)) minimum_layer = height_map(id, jd);
                    num_red += support_map(id, jd);
                }
            }

            if(num_red > 0)
                platform(ir, ic) = slicer.layer_pin_height(minimum_layer);
        }
    }

    return;
}

void MeshLayoutBase::transform_map_xz(MatrixXi &map, double x, double z, int fill) {
    Eigen::MatrixXi new_map;
    int nc = settings.pillar_column  * settings.xy_sample_num_each_pin;
    int nr = settings.pillar_row  * settings.xy_sample_num_each_pin;

    //init
    new_map =  Eigen::MatrixXi::Ones(nr, nc) * fill;

    for(int ir = 0; ir < nr; ir++)
    {
        for (int ic = 0; ic < nc; ic++)
        {
            Eigen::Vector2d p(settings.pin_center_x(ic), settings.pin_center_y(ir));
            p+= Eigen::Vector2d(settings.mm2int(x), settings.mm2int(z));

            int new_ir = std::floor(p(1) / settings.sample_width);
            int new_ic = std::floor(p(0) / settings.sample_width);
            if(0 <= new_ir && new_ir < nr && 0 <= new_ic && new_ic < nc)
                new_map(new_ir, new_ic) = map(ir, ic);
        }
    }

    map = new_map;
    return;
}

void MeshLayoutBase::rotate_map_yaxis(MatrixXi &map, double angle, Vector2d &center, int fill) {
    Eigen::MatrixXi new_map;
    int nc = settings.pillar_column  * settings.xy_sample_num_each_pin;
    int nr = settings.pillar_row  * settings.xy_sample_num_each_pin;

    //init
    new_map =  Eigen::MatrixXi::Ones(nr, nc) * fill;

    //center
    Eigen::Matrix2d rot_mat;
    rot_mat <<  cos(angle), -sin(angle),
                sin(angle),  cos(angle);
    for(int ir = 0; ir < nr; ir++)
    {
        for (int ic = 0; ic < nc; ic++)
        {
            Eigen::Vector2d p(settings.pin_center_x(ic), settings.pin_center_y(ir));
            p -= center;
            p = rot_mat * p;
            p += center;

            int new_ir = std::floor(p(1) / settings.sample_width);
            int new_ic = std::floor(p(0) / settings.sample_width);
            if(0 <= new_ir && new_ir < nr && 0 <= new_ic && new_ic < nc)
                new_map(new_ir, new_ic) = map(ir, ic);
        }
    }

    map = new_map;
    return;
}

void MeshLayoutBase::get_volume_support(MatrixXi &hmap, MatrixXi &smap, MatrixXd &platform, MatrixXd &V, MatrixXi &F)
{
    SceneOrganizer organizer;
    VecMatrixXi Fs;
    VecMatrixXd Vs;
    for(int ir = 0; ir < smap.rows(); ir++)
    {
        for(int ic = 0; ic < smap.cols(); ic++)
        {
            int pin_r = ir / settings.xy_sample_num_each_pin;
            int pin_c = ic / settings.xy_sample_num_each_pin;
            if(smap(ir, ic) > 0)
            {
                Eigen::MatrixXd tV;
                Eigen::MatrixXi tF;
                Eigen::Vector3d bottom(settings.int2mm(settings.pin_center_x(ic)),
                                       platform(pin_r, pin_c),
                                       settings.int2mm(settings.pin_center_y(ir)));
                double height = slicer.layer_height(hmap(ir, ic)) - platform(pin_r, pin_c) - settings.layer_height;

                if(height >= settings.layer_height) {
                    organizer.draw_pillar(bottom, settings.pad_size / settings.xy_sample_num_each_pin, height, tV, tF);
                    Vs.push_back(tV);
                    Fs.push_back(tF);
                }
            }
        }
    }
    V.setZero();
    F.setZero();
    int nv = 0, nf = 0;
    for(int id = 0; id < Vs.size(); id ++)
    {
        nv += Vs[id].rows(); nf += Fs[id].rows();
    }
    V.resize(nv, 3); F.resize(nf, 3);
    nv = 0; nf = 0;
    for(int id = 0; id < Vs.size(); id ++)
    {
        //F
        for(int jd = 0; jd < Fs[id].rows(); jd++)
        {
            F.row(nf ++) = Fs[id].row(jd) + Eigen::RowVector3i(nv, nv, nv);
        }

        //V
        for(int jd = 0; jd < Vs[id].rows(); jd++)
        {
            V.row(nv ++) = Vs[id].row(jd);
        }
    }

    return;
}

#endif //SUPPORTER_MESH_LAYOUT_BASE_H
