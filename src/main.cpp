#include <igl/viewer/Viewer.h>
#include <igl/triangle/triangulate.h>
#include <nanogui/formhelper.h>
#include <nanogui/screen.h>
#include <iostream>
#include "igl/readOBJ.h"
#include "igl/readSTL.h"

#include "testing_models_path.h"
#include "slice.h"
#include "clipper.hpp"
#include "normalizing_model.h"
#include "scene_organizer.h"
#include "slice_overhang_detector.h"
#include "wall_support_generator.h"
using namespace ClipperLib;
using std::cout;
using std::endl;
// Input polygon
Eigen::MatrixXd V;
Eigen::MatrixXi E, F;
Eigen::MatrixXd H;

// Triangulated interior
Eigen::MatrixXd V2;
Eigen::MatrixXi F2;

int layer;
igl::viewer::Viewer viewer;
Wall_Support_Generator slicer;
Settings settings;

void loadModel()
{
    // Load a mesh in OBJ format
    Eigen::MatrixXd temp_V, N;
    bool success = igl::readSTL(TESTING_MODELS_PATH "/timber/timber.stl",V,F,N);
    NormalizingModel normaler;
    normaler.size_normalize(V);
}

void slicing()
{
    slicer.set_mesh(V, F);
    slicer.contour_construction();
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
    for(int id = 0; id < slicer.number_layer(); id+= 5)
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

void overhang_slicing()
{
    SceneOrganizer organizer;
    std::vector<ClipperLib::Paths> slices;
    slicer.removing_overlap(slices);
    for(int id = 0; id < slicer.number_layer(); id++)
    {
        slicer.get_intersecting_surface(slices, id, V, F);
        if(V.rows() >= 3)
            organizer.add_mesh(V, F, Eigen::RowVector3d(1, 1, 0));
    }
    Eigen::MatrixXd C;
    organizer.get_mesh(V, F, C);
    viewer.data.clear();
    viewer.data.set_mesh(V, F);
    return;
}

void sampling()
{
    loadModel();
    viewer.data.clear();
    viewer.core.point_size = 5;
    viewer.data.set_mesh(V, F);
    Eigen::MatrixXd SP;
    slicer.sampling(SP);
    viewer.data.add_points(SP, Eigen::RowVector3d(1, 0 ,0));
}

void platform_n_mesh()
{
    loadModel();
    viewer.data.clear();
    SceneOrganizer organizer;
    organizer.add_mesh(V, F, Eigen::RowVector3d(1, 1 ,0));

    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(9, 11);
    Eigen::MatrixXi E;
    Eigen::MatrixXd SP;
    slicer.draw_platform(H);
    slicer.draw_sp(SP);
    slicer.draw_sp_lines(E);

    GeneratingPlatform platform_builder;
    platform_builder.draw_platform(V, F, H);

    organizer.add_mesh(V, F, Eigen::RowVector3d(0, 0 ,1));

    Eigen::MatrixXd C;
    organizer.get_mesh(V, F, C);
    viewer.data.set_face_based(true);
    viewer.data.set_mesh(V, F);

    viewer.data.set_colors(C);
    viewer.core.point_size = 3;
    viewer.data.add_points(SP, Eigen::RowVector3d(1, 0 ,0));

    for (int id = 0;id < E.rows(); id++)
        viewer.data.add_edges(SP.row(E(id, 0)), SP.row(E(id, 1)), Eigen::RowVector3d(1, 0, 0));
}

void platform()
{
    loadModel();
    viewer.data.clear();
    SceneOrganizer organizer;

    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(9, 11);
    Eigen::MatrixXi E;
    Eigen::MatrixXd SP;
    slicer.draw_platform(H);
    slicer.draw_sp(SP);
    slicer.draw_sp_lines(E);

    GeneratingPlatform platform_builder;
    platform_builder.draw_platform(V, F, H);

    organizer.add_mesh(V, F, Eigen::RowVector3d(0, 0 ,1));

    Eigen::MatrixXd C;
    organizer.get_mesh(V, F, C);
    viewer.data.set_face_based(true);
    viewer.data.set_mesh(V, F);

    viewer.data.set_colors(C);
    viewer.core.point_size = 3;
    viewer.data.add_points(SP, Eigen::RowVector3d(1, 0 ,0));

//    for (int id = 0;id < E.rows(); id++)
//        viewer.data.add_edges(SP.row(E(id, 0)), SP.row(E(id, 1)), Eigen::RowVector3d(1, 0, 0));
}

void generating_support()
{
    loadModel();
    viewer.data.clear();
    SceneOrganizer organizer;
    organizer.add_mesh(V, F, Eigen::RowVector3d(1, 1, 0));

    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(9, 11);
    Eigen::MatrixXd SP;
    Eigen::MatrixXi E;
    slicer.draw_platform(H);
    slicer.draw_sp(SP);
    slicer.draw_sp_lines(E);
    slicer.draw_support(V, F);
    organizer.add_mesh(V, F, Eigen::RowVector3d(1, 0, 1));

    GeneratingPlatform platform_builder;
    platform_builder.draw_platform(V, F, H);

    organizer.add_mesh(V, F, Eigen::RowVector3d(0, 0 ,1));

    Eigen::MatrixXd C;
    organizer.get_mesh(V, F, C);
    viewer.data.set_face_based(true);
    viewer.data.set_mesh(V, F);

    viewer.data.set_colors(C);
    viewer.core.point_size = 3;
    viewer.data.add_points(SP, Eigen::RowVector3d(1, 0 ,0));

}

void convex_hull()
{
    loadModel();
    viewer.data.clear();
    SceneOrganizer organizer;
    //organizer.add_mesh(V, F, Eigen::RowVector3d(1, 1, 0));

    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(9, 11);
    Eigen::MatrixXd SP;
    std::vector<std::vector<int>> convex;
    slicer.draw_platform(H);
    slicer.draw_sp(SP);
    slicer.convex_hull(convex);

    GeneratingPlatform platform_builder;
    platform_builder.draw_platform(V, F, H);

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
                    SP.row(convex[id][jd] + slicer.get_num_vertices_before(id)),
                    SP.row(convex[id][(jd + 1) % convex[id].size()] + slicer.get_num_vertices_before(id)),
                    Eigen::RowVector3d(0, 1, 0));
        }
    }
}

void level_set()
{
    loadModel();
    viewer.data.clear();
    SceneOrganizer organizer;
    std::vector<Fermat_Level_Set> fermat;
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(9, 11);
    slicer.level_set(fermat);
    slicer.draw_platform(H);

    GeneratingPlatform platform_builder;
    platform_builder.draw_platform(V, F, H);

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
    loadModel();
    viewer.data.clear();
    SceneOrganizer organizer;
    //organizer.add_mesh(V, F, Eigen::RowVector3d(1, 1 ,0));
    std::vector<std::list<FermatEdge>> path_layer;
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(9, 11);
    slicer.fermat_spiral(path_layer);
    slicer.draw_platform(H);

    GeneratingPlatform platform_builder;
    platform_builder.draw_platform(V, F, H);

    organizer.add_mesh(V, F, Eigen::RowVector3d(204.0f/255, 204.0/255 ,255.0f/255));

    Eigen::MatrixXd C;
    organizer.get_mesh(V, F, C);
    viewer.data.set_face_based(true);
    viewer.data.set_mesh(V, F);
    viewer.data.set_colors(C);

    for(int id = 0; id < slicer.number_layer(); id+= 2)
    {
        std::list<FermatEdge> path = path_layer[id];
        for(std::list<FermatEdge>::iterator it = path.begin(); it != path.end(); ++it)
        {
            viewer.data.add_edges(
                    it->P0(),
                    it->P1(),
                    Eigen::RowVector3d(1, 0, 0));
        }
    }

}

void heightmap()
{
    viewer.data.clear();
    Eigen::MatrixXd hmap;
    Eigen::MatrixXi smap;
    slicer.mesh_height_map(hmap, smap);

    Eigen::MatrixXd V;
    V.resize(hmap.size(), 3);

    std::cout << smap << std::endl;

    int kd = 0;
    for(int id = 0; id < hmap.rows(); id++)
    {
        for(int jd = 0; jd < hmap.cols(); jd++)
        {
            double sq_width = settings.pad_size / settings.xy_sample_num_each_pin;
            V(kd, 0) = jd * sq_width + sq_width / 2;
            V(kd, 2) = id * sq_width + sq_width / 2;
            V(kd, 1) = hmap(id, jd) < settings.maximum_height_map ? hmap(id, jd) : 0;
            //V(kd, 1) = smap(id, jd) ? settings.pillar_standard_height * 5 : 0;
            kd++;
        }
    }
    viewer.core.point_size = 1;
    viewer.data.add_points(V, Eigen::RowVector3d(1, 0 ,0));
}

void layout_optimization()
{
    loadModel();
    SceneOrganizer organizer;
    viewer.data.clear();
    Eigen::MatrixXd hmap;
    Eigen::MatrixXi smap;
    slicer.mesh_height_map(hmap, smap);
    double dx = 0, dz = 0;
    slicer.layout_optimization_xy(dx, dz, hmap, smap);


    for(int id = 0; id < V.rows(); id++)
    {
        V(id, 0) = V(id, 0) + dx;
        V(id, 2) = V(id, 2) + dz;
    }

    organizer.add_mesh(V, F, Eigen::RowVector3d(1, 1 ,0));

    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(9, 11);

    slicer = Wall_Support_Generator();
    slicer.set_mesh(V, F);
    slicer.contour_construction();
    slicer.draw_platform(H);

    GeneratingPlatform platform_builder;
    platform_builder.draw_platform(V, F, H);
    organizer.add_mesh(V, F, Eigen::RowVector3d(0, 0 ,1));

    Eigen::MatrixXd C;
    organizer.get_mesh(V, F, C);
    viewer.data.set_face_based(true);
    viewer.data.set_mesh(V, F);
    viewer.data.set_colors(C);
}

int main(int argc, char *argv[])
{
    loadModel();
    slicing();
    viewer.callback_init = [&](igl::viewer::Viewer& viewer)
    {
        viewer.ngui->addWindow(Eigen::Vector2i(220,10),"Sampling Control");
        viewer.ngui->addGroup("Slicing Settings");
        viewer.ngui->addVariable("layer", layer);
        viewer.ngui->addButton("slice", show_slice);
        viewer.ngui->addButton("recovery", recover);
        viewer.ngui->addButton("series slicing", series_slicing);
        viewer.ngui->addButton("overhang slicing", overhang_slicing);
        viewer.ngui->addButton("sampling", sampling);

        viewer.ngui->addWindow(Eigen::Vector2i(400,10),"Support Generation");
        viewer.ngui->addButton("platform & mesh", platform_n_mesh);
        viewer.ngui->addButton("platform", platform);
        viewer.ngui->addButton("support", generating_support);
        viewer.ngui->addButton("convex hull", convex_hull);
        viewer.ngui->addButton("level set", level_set);
        viewer.ngui->addButton("fermat", fermat_spiral);

        viewer.ngui->addWindow(Eigen::Vector2i(600,10),"Layout");
        viewer.ngui->addButton("Height Map", heightmap);
        viewer.ngui->addButton("Layout Optimization", layout_optimization);
        viewer.screen->performLayout();
        return false;
    };

    viewer.callback_key_down = &key_down;

    //draw();
    viewer.data.set_mesh(V,F);
    viewer.launch();
}
