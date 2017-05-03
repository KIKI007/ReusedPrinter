//
// Created by 汪子琦 on 5/2/17.
//

#ifndef SUPPORTER_MENU_FUNCTION_H
#define SUPPORTER_MENU_FUNCTION_H

//Rendering
#include <igl/viewer/Viewer.h>
#include <igl/triangle/triangulate.h>

//I/O
#include "igl/read_triangle_mesh.h"
#include "igl/readSTL.h"
#include "igl/writeSTL.h"
#include "igl/write_triangle_mesh.h"
#include <igl/png/writePNG.h>
#include <igl/png/readPNG.h>
#include <igl/jet.h>

//Std
#include <iostream>

//Self Class
#include "mesh_slicer_base.h"
#include "mesh_slicer_overhang.h"
#include "mesh_layout_base.h"
#include "mesh_layout_opt.h"
#include "scene_organizer.h"
#include "clipper.hpp"

//Model path
#include "testing_models_path.h"

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using ClipperLib::Clipper;
using ClipperLib::Paths;
using ClipperLib::Path;
using ClipperLib::ClipperOffset;
using ClipperLib::pftPositive;
using ClipperLib::ptClip;
using ClipperLib::ptSubject;
using Eigen::RowVector3d;

class NormalizingModel
{
public:
    NormalizingModel()
    {

    }
public:

    void size_normalize(Eigen::MatrixXd &V)
    {
        Eigen::Vector3d m = V.colwise().minCoeff();
        Eigen::Vector3d M = V.colwise().maxCoeff();

        //scale
        //V  *= (settings_.maximum_model_heigh / (M(1) - m(1)));

        m = V.colwise().minCoeff();
        M = V.colwise().maxCoeff();
        for(size_t kd = 0; kd < V.rows(); kd++)
        {

            //transform
            V.row(kd) -=  Eigen::RowVector3d((m(0) + M(0)) / 2 , m(1), (m(2) + M(2)) / 2 );
            V.row(kd) +=  Eigen::RowVector3d(11.0 / 2 * settings_.pad_size, 0, 9.0 / 2  * settings_.pad_size);

        }
    }

private:

    Settings settings_;
};


string file(string path, string name, string ext)
{
    return path + '/' + name + '/' + name + '.' + ext;
}

void loadModel(MatrixXd &V, MatrixXi &F, string file_path)
{
    Eigen::MatrixXd temp_V, N;
    bool success = igl::readSTL(file_path, V, F, N);
    if(success)
    {
        NormalizingModel normaler;
        normaler.size_normalize(V);
    }
    return;
}

bool slicing_for_one_layer(MeshSlicerBase &slicer, int &layer, MatrixXd &V, MatrixXi &F)
{
    if(slicer.empty()) return false;

    if(layer > slicer.number_layer() - 1)
        layer = slicer.number_layer() - 1;

    slicer.get_intersecting_surface(layer, V, F);

    if(V.rows() >= 3 && F.rows() >= 1)
        return  true;
    return false;
}

bool slicing_in_step(MeshSlicerBase &slicer, int layer_step, MatrixXd &V, MatrixXi &F)
{
    if(slicer.empty()) return false;

    SceneOrganizer organizer;
    for(int layer_id = 0; layer_id < slicer.number_layer(); layer_id+= layer_step)
    {
        slicer.get_intersecting_surface(layer_id, V, F);
        if(V.rows() >= 3) organizer.add_mesh(V, F, Eigen::RowVector3d(1, 1, 0));
    }

    MatrixXd C;
    organizer.get_mesh(V, F, C);

    if(V.rows() >= 3 && F.rows() >= 1)
        return  true;
    return false;
}

bool slicing_bottom_half(MeshSlicerOverhang &slicer, MatrixXd &V, MatrixXi &F)
{
    if(slicer.empty()) return false;

    vector<Paths> bottom_half;
    slicer.get_bottom_half_overhang(bottom_half);

    SceneOrganizer organizer;
    for(int layer_id = 0; layer_id < slicer.number_layer(); layer_id++)
    {
        slicer.get_intersecting_surface(bottom_half, layer_id, V, F);
        if(V.rows() >= 3) organizer.add_mesh(V, F, Eigen::RowVector3d(1, 1, 0));
    }

    MatrixXd C;
    organizer.get_mesh(V, F, C);

    if(V.rows() >= 3 && F.rows() >= 1)
        return  true;
    return false;
}

bool slicing_upper_half(MeshSlicerOverhang &slicer, MatrixXd &V, MatrixXi &F)
{
    if(slicer.empty()) return false;

    vector<Paths> upper_half;
    slicer.get_upper_half_overhang(upper_half);

    SceneOrganizer organizer;
    for(int layer_id = 0; layer_id < slicer.number_layer(); layer_id++)
    {
        slicer.get_intersecting_surface(upper_half, layer_id, V, F);
        if(V.rows() >= 3) organizer.add_mesh(V, F, Eigen::RowVector3d(1, 1, 0));
    }

    MatrixXd C;
    organizer.get_mesh(V, F, C);

    if(V.rows() >= 3 && F.rows() >= 1)
        return  true;
    return false;
}

bool move_xy_opt(MeshSlicerOverhang &slicer, MatrixXd &V, MatrixXi &F, MatrixXd &C)
{
    if(slicer.empty()) return false;

    MeshLayoutOpt layouter;
    layouter.set_slicer(slicer);
    SceneOrganizer organizer;

    LayoutOptResult result;
    layouter.request_layout_xz_opt(result);

    MatrixXd spV, plV;
    MatrixXi spF, plF;

    //platform
    MatrixXd platform; layouter.get_opt_platform(platform);
    GeneratingPlatform platform_builder;
    platform_builder.draw_platform(plV, plF, platform);
    organizer.add_mesh(plV, plF, RowVector3d(0, 0 , 1));

    //support
    MatrixXi hmap, smap;
    layouter.get_opt_hmap(hmap);
    layouter.get_opt_smap(smap);
    layouter.get_volume_support(hmap, smap, platform, spV, spF);
    organizer.add_mesh(spV, spF, RowVector3d(1, 0 , 0));

    //mesh
    slicer.get_vertices_mat(V);
    slicer.get_faces_mat(F);
    for(int id = 0; id < V.rows(); id++)
    {
        V(id, 0) += result.dx;
        V(id, 2) += result.dz;
    }
    organizer.add_mesh(V, F, RowVector3d(1, 1, 0));

    organizer.get_mesh(V, F, C);

    return true;
}

#else
void loadModel()
{
    // Load a mesh in OBJ format
    Eigen::MatrixXd temp_V, N;
    //std::cout << TESTING_MODELS_PATH "/" + pathname + "/" + pathname + ".stl";
    bool success = igl::readSTL(TESTING_MODELS_PATH "/" + pathname + "/" + pathname + ".stl",V,F,N);
    if(success)
    {
        NormalizingModel normaler;
        normaler.size_normalize(V);
        viewer.data.clear();
        viewer.data.set_mesh(V, F);

        slicer.clear();
        slicer.set_mesh(V, F);
        slicer.contour_construction();
        platform.setZero();
    }
}

void slicing()
{
    slicer.set_mesh(V, F);
    slicer.contour_construction();
}

void image()
{
    for(int id = 0; id < slicer.number_layer(); id+= menu_input.layer_step)
    {
        std::vector<ClipperLib::Paths> slices;
        slicer.get_slices(slices);
        Eigen::MatrixXi imag;
        ScanLineFill filler;
        filler.polygon_fill(slices[id], imag);
        int nr = imag.rows();
        int nc = imag.cols();
        Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R(nc, nr);
        Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> G(nc ,nr);
        Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> B(nc ,nr);
        Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> A(nc ,nr);
        for(int ir = 0; ir < nr; ir++)
        {
            for(int ic = 0; ic < nc; ic++)
            {
                A(ic, ir) = 255;
                if(imag(ir, ic))
                    R(ic, ir) = G(ic, ir) = B(ic, ir) = 0;
                else
                    R(ic, ir) = G(ic, ir) = B(ic, ir) = 255;
            }
        }

        std::string output_path = TESTING_MODELS_PATH "/" + pathname + "/image/";
        char num[100];
        sprintf(num, "%d", id);
        output_path = output_path + num + ".png";

        igl::png::writePNG(R,G,B,A,output_path);
    }
}

void draw_image(Eigen::MatrixXd map)
{
    int nr = map.rows();
    int nc = map.cols();
    typedef Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> ColorMatrix;
    ColorMatrix R,G,B,A;
    A = ColorMatrix::Zero(nc, nr);
    R = ColorMatrix::Zero(nc, nr);
    G = ColorMatrix::Zero(nc, nr);
    B = ColorMatrix::Zero(nc, nr);
    for(int id = 0; id < nr; id++)
    {
        for(int jd = 0; jd < nc; jd++)
        {
            double h = map(id, jd) / settings.maximum_height_map;
            double r, g, b;
            igl::jet(h, r, g, b);
            R(jd, id) = r * 255;
            B(jd, id) = g * 255;
            G(jd, id) = b * 255;
            A(jd, id) = 255;
        }
    }

    std::string output_path = TESTING_MODELS_PATH "/" + pathname + "/image/";
    char num[100];
    output_path = output_path + "height.png";

    igl::png::writePNG(R,G,B,A,output_path);
}

void show_slice()
{
    if(menu_input.layer > slicer.number_layer() - 1)
        menu_input.layer = slicer.number_layer() - 1;
    slicer.get_intersecting_surface(menu_input.layer, V, F);
//    std::cout << V << F << std::endl;
    if(V.rows() >= 3)
    {
        viewer.data.clear();
        viewer.data.set_mesh(V, F);
    }
    return;
}

bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier)
{
    std::cout<<"Key: "<<key<<" "<<(unsigned int)key<<std::endl;
    if (key == 10)
    {
        show_slice();
    }
    return false;
}

void recover()
{
    loadModel();
    viewer.data.clear();
    viewer.data.set_mesh(V, F);
    return;
}

void series_slicing()
{
    SceneOrganizer organizer;
    for(int id = 0; id < slicer.number_layer(); id+= menu_input.layer_step)
    {
        slicer.get_intersecting_surface(id, V, F);
        if(V.rows() >= 3)
            organizer.add_mesh(V, F, Eigen::RowVector3d(1, 1, 0));
    }
    Eigen::MatrixXd C;
    organizer.get_mesh(V, F, C);
    viewer.data.clear();
    viewer.data.set_mesh(V, F);
    return;
}

void overhang_bottom_slicing()
{
    SceneOrganizer organizer;
    std::vector<ClipperLib::Paths> bslice,uslice;
    slicer.get_bottom_half(bslice);
    for(int id = 0; id < slicer.number_layer(); id++)
    {
        slicer.get_intersecting_surface(bslice, id, V, F);
        if(V.rows() >= 3)
            organizer.add_mesh(V, F, Eigen::RowVector3d(1, 1, 0));
    }
    Eigen::MatrixXd C;
    organizer.get_mesh(V, F, C);
    viewer.data.clear();
    if(V.rows() >= 3)
        viewer.data.set_mesh(V, F);
    return;
}

void output_screen()
{
    igl::writeSTL(TESTING_MODELS_PATH "/" + pathname + "/" + pathname + "_tmp.stl",V, F);
    return;
}

void output_model(Eigen::MatrixXd oV, Eigen::MatrixXi oF, std::string aprex)
{
    for(int id = 0; id < oV.rows(); id++)
    {
        double y = oV(id, 1);
        double z = oV(id, 2);
        oV(id, 0) *= 0.03;
        oV(id, 1) = -z * 0.03;
        oV(id, 2) = y * 0.03;
    }
    igl::writeSTL(TESTING_MODELS_PATH "/" + pathname + "/" + pathname + "_" + aprex + ".stl",oV, oF);
    return;
}

void output_printing()
{
    loadModel();
    SceneOrganizer organizer;

    //layout
    MeshLayout layout;
    LayoutOptOutput opt;
    layout.set_slicer(slicer);
    opt = layout.xy_layout();
    slicer.move_XY(opt.dx, opt.dy);
    platform = opt.platform;

    //output mesh for rendering
    igl::writeSTL(TESTING_MODELS_PATH "/" + pathname + "/" + pathname + "_render.stl",V, F);

    //output mesh
    for(int id = 0; id < V.rows(); id++)
    {
        double y = V(id, 1);
        double z = V(id, 2);
        V(id, 0) += settings.platform_zero_x;
        V(id, 1) = -z + settings.platform_zero_y;
        V(id, 2) = y;
    }
    igl::writeSTL(TESTING_MODELS_PATH "/" + pathname + "/" + pathname + "_print.stl",V, F);

    //ouput platform
    GeneratingPlatform platform_builder;
    platform_builder.draw_platform(V, F, platform);
    igl::writeSTL(TESTING_MODELS_PATH "/" + pathname + "/" + pathname + "_platform.stl",V, F);

}

void convex_hull()
{
    loadModel();
    viewer.data.clear();
    SceneOrganizer organizer;

    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(9, 11);
    Eigen::MatrixXd SP;
    std::vector<std::vector<ConvexHullPoint>> convex;

    MeshSupport support;
    support.sp_pin_construction(slicer, platform);
    std::cout << "done" << std::endl;
    support.convex_hull_construction(convex);

    slicer.get_vertices_mat(V);
    organizer.add_mesh(V, F, Eigen::RowVector3d(1, 1, 0));

    GeneratingPlatform platform_builder;
    platform_builder.draw_platform(V, F, platform);
    organizer.add_mesh(V, F, Eigen::RowVector3d(0, 0 ,1));

    Eigen::MatrixXd C;
    organizer.get_mesh(V, F, C);
    viewer.data.set_face_based(true);
    viewer.data.set_mesh(V, F);

    viewer.data.set_colors(C);
    viewer.core.point_size = 3;
    viewer.data.add_points(SP, Eigen::RowVector3d(1, 0 ,0));
    viewer.core.line_width = 2;

    for(int id = 0; id < convex.size(); id++)
    {
        for(int jd = 0; jd < convex[id].size(); jd++)
        {
            viewer.data.add_edges(
                    convex[id][jd].PT(),
                    convex[id][(jd + 1) % convex[id].size()].PT(),
                    Eigen::RowVector3d(0, 1, 0));
        }
    }

    return;
}

void level_set()
{
    loadModel();
    viewer.data.clear();
    SceneOrganizer organizer;
    std::vector<Fermat_Level_Set> fermat;

    MeshSupport support;
    support.sp_pin_construction(slicer, platform);
    support.level_set(fermat);
    slicer.get_vertices_mat(V);
    organizer.add_mesh(V, F, Eigen::RowVector3d(1, 1 ,0));

    GeneratingPlatform platform_builder;
    platform_builder.draw_platform(V, F, platform);
    organizer.add_mesh(V, F, Eigen::RowVector3d(0, 0 ,1));

    Eigen::MatrixXd C;
    organizer.get_mesh(V, F, C);
    viewer.data.set_face_based(true);
    viewer.data.set_mesh(V, F);
    viewer.data.set_colors(C);

    for(int id = 0; id < fermat.size(); id++)
        for(int jd = 0; jd < fermat[id].num_level(); jd++)
        {
            std::vector<FermatEdge> polygon;
            fermat[id].get_level(polygon, jd);
            for(int kd = 0; kd < polygon.size(); kd++)
                viewer.data.add_edges(
                        polygon[kd].P0(),
                        polygon[kd].P1(),
                        Eigen::RowVector3d(0, 1, 0));
        }
}

void fermat_spiral()
{
    settings.tic("Time");
    loadModel();
    viewer.data.clear();
    SceneOrganizer organizer;
    std::vector<std::list<FermatEdge>> path_layer;
    MeshSupport support;
    support.fermat_spiral(path_layer, slicer, platform);
    settings.toc();

    slicer.get_vertices_mat(V);
    organizer.add_mesh(V, F, Eigen::RowVector3d(1, 1 ,0));

    GeneratingPlatform platform_builder;
    platform_builder.draw_platform(V, F, platform);
    organizer.add_mesh(V, F, Eigen::RowVector3d(204.0f/255, 204.0/255 ,255.0f/255));

    Eigen::MatrixXd C;
    organizer.get_mesh(V, F, C);
    viewer.data.set_face_based(true);
    viewer.data.set_mesh(V, F);
    viewer.data.set_colors(C);
//
    int id = 1;
    for(int id = 0; id < path_layer.size(); id+= 10) {
        std::list<FermatEdge> path = path_layer[id];
        for (std::list<FermatEdge>::iterator it = path.begin(); it != path.end(); ++it) {
            viewer.data.add_edges(
                    it->P0(),
                    it->P1(),
                    Eigen::RowVector3d(1, 0, 0));
        }
    }
    std::cout << "finished" << std::endl;
}
//
void heightmap()
{
//    Eigen::MatrixXd hmap;
//    Eigen::MatrixXi smap;
//    MeshLayout layout;
//    double dx, dy;
//    layout.set_slicer(slicer);
//    layout.get_height_map(hmap);
//
//    Eigen::MatrixXd V;
//    V.resize(hmap.size(), 3);
//
//    //std::cout << hmap << std::endl;
//
//    int kd = 0;
//    for(int id = 0; id < hmap.rows(); id++)
//    {
//        for(int jd = 0; jd < hmap.cols(); jd++)
//        {
//            double sq_width = settings.pad_size / settings.xy_sample_num_each_pin;
//            V(kd, 0) = jd * sq_width + sq_width / 2;
//            V(kd, 2) = id * sq_width + sq_width / 2;
//            V(kd, 1) = hmap(id, jd);
//            //V(kd, 1) = smap(id, jd) ? settings.pillar_standard_height * 5 : 0;
//            kd++;
//        }
//    }
//    viewer.core.point_size = 2;
//    viewer.data.add_points(V, Eigen::RowVector3d(1, 0 ,0));

    Eigen::MatrixXd hmap;
    MeshLayout layout;
    layout.set_slicer(slicer);
    layout.get_height_map(hmap);
    draw_image(hmap);
    return;
}

void redmap()
{
    Eigen::MatrixXi smap;
    MeshLayout layout;
    double dx, dy;
    layout.set_slicer(slicer);
    layout.get_red_map(smap);

    Eigen::MatrixXd V;
    V.resize(smap.size(), 3);

    int kd = 0;
    for(int id = 0; id < smap.rows(); id++)
    {
        for(int jd = 0; jd < smap.cols(); jd++)
        {
            double sq_width = settings.pad_size / settings.xy_sample_num_each_pin;
            V(kd, 0) = jd * sq_width + sq_width / 2;
            V(kd, 2) = id * sq_width + sq_width / 2;
            //V(kd, 1) = hmap(id, jd) < settings.maximum_height_map ? hmap(id, jd) : 0;
            V(kd, 1) = smap(id, jd) > 0 ? settings.pillar_standard_height * 5 : 0;
            kd++;
        }
    }
    viewer.core.point_size = 2;
    viewer.data.add_points(V, Eigen::RowVector3d(1, 0 ,0));
}

void xy_move()
{
    //loadModel();
    SceneOrganizer organizer;
    viewer.data.clear();
    MeshLayout layout;
    LayoutOptOutput opt;
    layout.set_slicer(slicer);
    opt = layout.xy_layout();

    platform = opt.platform;
    slicer.move_XY(opt.dx, opt.dy);

    slicer.get_vertices_mat(V);
    organizer.add_mesh(V, F, Eigen::RowVector3d(1, 1 ,0));

    GeneratingPlatform platform_builder;
    platform_builder.draw_platform(V, F, platform);
    organizer.add_mesh(V, F, Eigen::RowVector3d(0, 0 ,1));

    Eigen::MatrixXd C;
    organizer.get_mesh(V, F, C);
    viewer.data.set_face_based(true);
    viewer.data.set_mesh(V, F);
    viewer.data.set_colors(C);
}

void rotate_move()
{
    //loadModel();
    SceneOrganizer organizer;
    viewer.data.clear();
    MeshLayout layout;
    LayoutOptOutput opt;
    layout.set_slicer(slicer);
    opt = layout.rotate_layout();

    platform = opt.platform;
    slicer.rotate(opt.center, opt.angle);
    slicer.move_XY(opt.dx, opt.dy);

    slicer.get_vertices_mat(V);
    organizer.add_mesh(V, F, Eigen::RowVector3d(1, 1 ,0));

    GeneratingPlatform platform_builder;
    platform_builder.draw_platform(V, F, platform);
    organizer.add_mesh(V, F, Eigen::RowVector3d(0, 0 ,1));

    Eigen::MatrixXd C;
    organizer.get_mesh(V, F, C);
    viewer.data.set_face_based(true);
    viewer.data.set_mesh(V, F);
    viewer.data.set_colors(C);
}

void rotate()
{
    SceneOrganizer organizer;
    viewer.data.clear();

    Eigen::Vector2d center(0, 0);
    slicer.rotate(center, menu_input.angle / 180 * settings.PI);

    slicer.get_vertices_mat(V);
    organizer.add_mesh(V, F, Eigen::RowVector3d(1, 1 ,0));

    GeneratingPlatform platform_builder;
    platform_builder.draw_platform(V, F, platform);
    organizer.add_mesh(V, F, Eigen::RowVector3d(0, 0 ,1));

    Eigen::MatrixXd C;
    organizer.get_mesh(V, F, C);
    viewer.data.set_face_based(true);
    viewer.data.set_mesh(V, F);
    viewer.data.set_colors(C);
}

void gcode()
{
    std::cout << TESTING_MODELS_PATH  "/" + pathname + "/" + pathname + ".gcode" << std::endl;

    loadModel();
    viewer.data.clear();
    SceneOrganizer organizer;
    std::vector<std::list<FermatEdge>> path_layer;
    MeshSupport support;
    support.fermat_spiral(path_layer, slicer, platform);
    settings.toc();

    slicer.get_vertices_mat(V);
    organizer.add_mesh(V, F, Eigen::RowVector3d(1, 1 ,0));

    GeneratingPlatform platform_builder;
    platform_builder.draw_platform(V, F, platform);
    organizer.add_mesh(V, F, Eigen::RowVector3d(204.0f/255, 204.0/255 ,255.0f/255));

    Eigen::MatrixXd C;
    organizer.get_mesh(V, F, C);
    viewer.data.set_face_based(true);
    viewer.data.set_mesh(V, F);
    viewer.data.set_colors(C);
    Gcode gcode(TESTING_MODELS_PATH  "/"  + pathname + "/" + pathname + ".gcode");
    gcode.get_minXY();
    slicer.get_vertices_mat(V);
    std::cout << "Gcode: MinX " << gcode.min_X  << ",\tMinY " << gcode.min_Y << std::endl;
    std::cout << "Support: MinX "    << V.colwise().minCoeff()[0]  + settings.platform_zero_x
              << ",\tMinY "          << -V.colwise().maxCoeff()[2] + settings.platform_zero_y;
    gcode.move_XY( V.colwise().minCoeff()[0] + settings.platform_zero_x - gcode.min_X,
                   -V.colwise().maxCoeff()[2] + settings.platform_zero_y - gcode.min_Y);
    gcode.add_support(path_layer);
    gcode.print(TESTING_MODELS_PATH  "/"  + pathname + "/" + pathname + "_new.gcode");
}

void virtual_support(SupportType type)
{
    SceneOrganizer organizer;
    viewer.data.clear();
    MeshLayout layout;
    LayoutOptOutput opt;
    layout.set_slicer(slicer);
    Eigen::MatrixXd supportV;
    Eigen::MatrixXi supportF;
    if(type == Non)
        opt = layout.non_pin_support(supportV, supportF);
    else if(type == Pin)
        opt = layout.stand_support(supportV, supportF);
    else if(type == Pin_XY)
        opt = layout.xy_opt_support(supportV, supportF);

    output_model(supportV, supportF, "support");
    organizer.add_mesh(supportV, supportF, Eigen::RowVector3d(1, 0 ,0));

    platform = opt.platform;

    if(opt.rotate) slicer.rotate(opt.center, opt.angle);
    slicer.move_XY(opt.dx, opt.dy);

    slicer.get_vertices_mat(V);
    organizer.add_mesh(V, F, Eigen::RowVector3d(1, 1 ,0));
    output_model(V, F, "model");

    GeneratingPlatform platform_builder;
    platform_builder.draw_platform(V, F, platform);
    output_model(V, F, "platform");
    organizer.add_mesh(V, F, Eigen::RowVector3d(0, 0 ,1));

    Eigen::MatrixXd C;
    organizer.get_mesh(V, F, C);
    viewer.data.set_face_based(true);
    viewer.data.set_mesh(V, F);
    viewer.data.set_colors(C);
}

void none_pin_support()
{
    loadModel();
    slicer.move_XY(menu_input.layout_dx, menu_input.layout_dy);
    virtual_support(Non);
}

void pin_xy_support()
{
    loadModel();
    virtual_support(Pin_XY);
}

void pin_support()
{
    loadModel();
    std::cout << menu_input.layout_dx << ", " << menu_input.layout_dy << std::endl;
    slicer.move_XY(menu_input.layout_dx, menu_input.layout_dy);
    virtual_support(Pin);
}

#endif //SUPPORTER_MENU_FUNCTION_H