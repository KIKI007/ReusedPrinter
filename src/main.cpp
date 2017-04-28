#include <igl/viewer/Viewer.h>
#include <igl/triangle/triangulate.h>
#include <nanogui/formhelper.h>
#include <nanogui/screen.h>
#include <iostream>
#include "igl/read_triangle_mesh.h"
#include "igl/readSTL.h"
#include "igl/writeSTL.h"
#include "igl/write_triangle_mesh.h"

#include "testing_models_path.h"
#include "mesh_slicer.h"
#include "mesh_layout.h"
#include "mesh_support.h"
#include "clipper.hpp"
#include "normalizing_model.h"
#include "scene_organizer.h"
#include "gcode.h"
#include "scan_line_fill.h"

#include <igl/png/writePNG.h>
#include <igl/png/readPNG.h>
using namespace ClipperLib;
using std::cout;
using std::endl;
// Input polygon
Eigen::MatrixXd V;
Eigen::MatrixXi E, F;

// Triangulated interior
Eigen::MatrixXd V2;
Eigen::MatrixXi F2;

int layer;
igl::viewer::Viewer viewer;
MeshSlicer slicer;
Settings settings;
std::string pathname;
Eigen::MatrixXd platform;

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
    if(layer > slicer.number_layer() - 1)
        layer = slicer.number_layer() - 1;
    std::vector<ClipperLib::Paths> slices;
    slicer.get_slices(slices);
    Eigen::MatrixXd imag;
    ScanLineFill filler;
    filler.polygon_fill(slices[layer], imag);
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
    igl::png::writePNG(R,G,B,A,"/Users/wangziqi/Desktop/USC Spring/Code/Supporter/out.png");
}

void show_slice()
{
    if(layer > slicer.number_layer() - 1)
        layer = slicer.number_layer() - 1;
    slicer.get_intersecting_surface(layer, V, F);
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
    for(int id = 0; id < slicer.number_layer(); id+= 2)
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

void overhang_upper_slicing()
{
//    SceneOrganizer organizer;
//    std::vector<ClipperLib::Paths> bslice,uslice;
//    slicer.get_overhang_slices(bslice, uslice);
//    for(int id = 0; id < slicer.number_layer(); id++)
//    {
//        slicer.get_intersecting_surface(uslice, id, V, F);
//        if(V.rows() >= 3)
//            organizer.add_mesh(V, F, Eigen::RowVector3d(1, 1, 0));
//    }
//    Eigen::MatrixXd C;
//    organizer.get_mesh(V, F, C);
//    viewer.data.clear();
//    if(V.rows() >= 3)
//        viewer.data.set_mesh(V, F);
//    return;
}

//void platform_n_mesh()
//{
//    loadModel();
//    viewer.data.clear();
//    SceneOrganizer organizer;
//    organizer.add_mesh(V, F, Eigen::RowVector3d(1, 1 ,0));
//
//    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(9, 11);
//    Eigen::MatrixXi E;
//    Eigen::MatrixXd SP;
//    slicer.draw_platform(H);
//    slicer.draw_sp(SP);
//    slicer.draw_sp_lines(E);
//
//    GeneratingPlatform platform_builder;
//    platform_builder.draw_platform(V, F, H);
//
//    organizer.add_mesh(V, F, Eigen::RowVector3d(0, 0 ,1));
//
//    Eigen::MatrixXd C;
//    organizer.get_mesh(V, F, C);
//    viewer.data.set_face_based(true);
//    viewer.data.set_mesh(V, F);
//
//    viewer.data.set_colors(C);
//    viewer.core.point_size = 3;
//    viewer.data.add_points(SP, Eigen::RowVector3d(1, 0 ,0));
//
//    for (int id = 0;id < E.rows(); id++)
//        viewer.data.add_edges(SP.row(E(id, 0)), SP.row(E(id, 1)), Eigen::RowVector3d(1, 0, 0));
//}
//
void output()
{
    loadModel();
    SceneOrganizer organizer;

    //layout
    MeshLayout layout;
    double dx = 0, dz = 0;
    layout.set_slicer(slicer);
    layout.xy_layout(dx, dz, platform);
    slicer.move_XY(dx,dz);
    slicer.getV(V);

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

//void generating_support()
//{
//    loadModel();
//    viewer.data.clear();
//    SceneOrganizer organizer;
//    organizer.add_mesh(V, F, Eigen::RowVector3d(1, 1, 0));
//
//    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(9, 11);
//    Eigen::MatrixXd SP;
//    Eigen::MatrixXi E;
//    slicer.draw_platform(H);
//    slicer.draw_sp(SP);
//    slicer.draw_sp_lines(E);
//    slicer.draw_support(V, F);
//    organizer.add_mesh(V, F, Eigen::RowVector3d(1, 0, 1));
//
//    GeneratingPlatform platform_builder;
//    platform_builder.draw_platform(V, F, H);
//
//    organizer.add_mesh(V, F, Eigen::RowVector3d(0, 0 ,1));
//
//    Eigen::MatrixXd C;
//    organizer.get_mesh(V, F, C);
//    viewer.data.set_face_based(true);
//    viewer.data.set_mesh(V, F);
//
//    viewer.data.set_colors(C);
//    viewer.core.point_size = 3;
//    viewer.data.add_points(SP, Eigen::RowVector3d(1, 0 ,0));
//
//}
//
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

    slicer.getV(V);
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
    slicer.getV(V);
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

    slicer.getV(V);
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
    Eigen::MatrixXd hmap;
    Eigen::MatrixXi smap;
    MeshLayout layout;
    double dx, dy;
    layout.set_slicer(slicer);
    layout.get_height_map(hmap);

    Eigen::MatrixXd V;
    V.resize(hmap.size(), 3);

    //std::cout << hmap << std::endl;

    int kd = 0;
    for(int id = 0; id < hmap.rows(); id++)
    {
        for(int jd = 0; jd < hmap.cols(); jd++)
        {
            double sq_width = settings.pad_size / settings.xy_sample_num_each_pin;
            V(kd, 0) = jd * sq_width + sq_width / 2;
            V(kd, 2) = id * sq_width + sq_width / 2;
            V(kd, 1) = hmap(id, jd);
            //V(kd, 1) = smap(id, jd) ? settings.pillar_standard_height * 5 : 0;
            kd++;
        }
    }
    viewer.core.point_size = 2;
    viewer.data.add_points(V, Eigen::RowVector3d(1, 0 ,0));
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
    loadModel();
    SceneOrganizer organizer;
    viewer.data.clear();
    MeshLayout layout;
    double dx = 0, dz = 0;
    layout.set_slicer(slicer);
    layout.xy_layout(dx, dz, platform);
    slicer.move_XY(dx, dz);
    slicer.getV(V);
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

    slicer.getV(V);
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
    slicer.getV(V);
    std::cout << "Gcode: MinX " << gcode.min_X  << ",\tMinY " << gcode.min_Y << std::endl;
    std::cout << "Support: MinX "    << V.colwise().minCoeff()[0]  + settings.platform_zero_x
              << ",\tMinY "          << -V.colwise().maxCoeff()[2] + settings.platform_zero_y;
    gcode.move_XY( V.colwise().minCoeff()[0] + settings.platform_zero_x - gcode.min_X,
                  -V.colwise().maxCoeff()[2] + settings.platform_zero_y - gcode.min_Y);
    gcode.add_support(path_layer);
    gcode.print(TESTING_MODELS_PATH  "/"  + pathname + "/" + pathname + "_new.gcode");
}

int main(int argc, char *argv[])
{
    pathname = "arch";
    loadModel();
    viewer.callback_init = [&](igl::viewer::Viewer& viewer)
    {
        viewer.ngui->addGroup("Model Loading");
        viewer.ngui->addVariable("Model Name", pathname);
        viewer.ngui->addButton("Load Model", loadModel);
        viewer.ngui->addButton("Output", output);
        viewer.ngui->addButton("Gcode", gcode);

        viewer.ngui->addWindow(Eigen::Vector2i(220,10),"Slicing & Overhang");
        viewer.ngui->addGroup("Slicing Settings");
        viewer.ngui->addVariable("layer", layer);
        viewer.ngui->addButton("slice", show_slice);
        viewer.ngui->addButton("recovery", recover);
        viewer.ngui->addButton("series slicing", series_slicing);
        viewer.ngui->addButton("overhang bottom", overhang_bottom_slicing);
        viewer.ngui->addButton("overhang upper", overhang_upper_slicing);
        viewer.ngui->addButton("image", image);

        viewer.ngui->addWindow(Eigen::Vector2i(300,10),"Layout Optimization");
        viewer.ngui->addButton("Height Map", heightmap);
        viewer.ngui->addButton("Red Map", redmap);
        viewer.ngui->addButton("XY Layout", xy_move);

        viewer.ngui->addWindow(Eigen::Vector2i(400,10),"Support");
//        viewer.ngui->addButton("platform & mesh", platform_n_mesh);
//        viewer.ngui->addButton("support", generating_support);
        viewer.ngui->addButton("convex hull", convex_hull);
        viewer.ngui->addButton("level set", level_set);
        viewer.ngui->addButton("fermat", fermat_spiral);

//        viewer.ngui->addWindow(Eigen::Vector2i(600,10),"Layout");
//        viewer.ngui->addButton("Height Map", heightmap);
//        viewer.ngui->addButton("Layout Optimization", layout_optimization);
        viewer.screen->performLayout();
        return false;
    };

    viewer.callback_key_down = &key_down;

    viewer.launch();
}
