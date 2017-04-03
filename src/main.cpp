//#include <igl/readOFF.h>
//#include <igl/viewer/Viewer.h>
//#include <sstream>
//#include "testing_models_path.h"
//#include "detecting_overhangs.h"
//#include "rendering_tree_support.h"
//#include "scene_organizer.h"
//#include "generating_support.h"
//#include "normalizing_model.h"
//
//Eigen::Vector3d m, M;
//Eigen::MatrixXi E_box(12,2);
//Eigen::MatrixXd V, V_Box;
//Eigen::MatrixXi F;
//
//void loadModel(Eigen::MatrixXd &V, Eigen::MatrixXi &F)
//{
//    // Load a mesh in OBJ format
//    igl::readOFF(TESTING_MODELS_PATH "/arch.off", V, F);
//    NormalizingModel normalize;
//    normalize.size_normalize(V);
//}
//
//void getBoundingBox(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &V_box)
//{
//    // Find the bounding box
//    m = V.colwise().minCoeff();
//    M = V.colwise().maxCoeff();
//
//    // Corners of the bounding box
//    V_box = Eigen::MatrixXd::Zero(8,3);
//
//    V_box <<
//          m(0), m(1), m(2),
//            M(0), m(1), m(2),
//            M(0), M(1), m(2),
//            m(0), M(1), m(2),
//            m(0), m(1), M(2),
//            M(0), m(1), M(2),
//            M(0), M(1), M(2),
//            m(0), M(1), M(2);
//
//    // Edges of the bounding box
//
//    E_box <<
//          0, 1,
//            1, 2,
//            2, 3,
//            3, 0,
//            4, 5,
//            5, 6,
//            6, 7,
//            7, 4,
//            0, 4,
//            1, 5,
//            2, 6,
//            7 ,3;
//}
//
//
//
//void setViewer(igl::viewer::Viewer &viewer, Eigen::MatrixXd V, Eigen::MatrixXi &F, Eigen::MatrixXd V_box) {
//    // Plot the mesh
//
//    viewer.data.set_mesh(V, F);
//
//    // Plot the corners of the bounding box as points
//    viewer.data.add_points(V_box,Eigen::RowVector3d(1,0,0));
//
//    // Plot the edges of the bounding box
//    for (unsigned i=0;i<E_box.rows(); ++i)
//        viewer.data.add_edges
//                (
//                        V_box.row(E_box(i,0)),
//                        V_box.row(E_box(i,1)),
//                        Eigen::RowVector3d(1,0,0)
//                );
//
//    // Plot labels with the coordinates of bounding box vertices
//    std::stringstream l1;
//    l1 << m(0) << ", " << m(1) << ", " << m(2);
//    viewer.data.add_label(m,l1.str());
//    std::stringstream l2;
//    l2 << M(0) << ", " << M(1) << ", " << M(2);
//    viewer.data.add_label(M,l2.str());
//}
//
//// This function is called every time a keyboard button is pressed
//bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier)
//{
//    std::cout<<"Key: "<<key<<" "<<(unsigned int)key<<std::endl;
//    if (key == '1')
//    {
//        loadModel(V, F);
//        // Clear should be called before drawing the mesh
//        viewer.data.clear();
//        // Draw_mesh creates or updates the vertices and faces of the displayed mesh.
//        // If a mesh is already displayed, draw_mesh returns an error if the given V and
//        // F have size different than the current ones
//        viewer.data.set_mesh(V, F);
//        getBoundingBox(V, F, V_Box);
//        viewer.core.align_camera_center(V, F);
//        viewer.data.set_face_based(false);
//        setViewer(viewer, V, F, V_Box);
//    }
//    else if (key == '2')
//    {
//        viewer.data.clear();
//        viewer.data.set_mesh(V, F);
//        viewer.core.align_camera_center(V,F);
//        DetectingOverhangs overhangs_detector;
//        Eigen::MatrixXd C;
//        overhangs_detector.colorOverhangingFaces(V, F, C);
//        viewer.data.set_face_based(true);
//        viewer.data.set_colors(C);
//    }
//    else if (key == '3')
//    {
//        viewer.data.clear();
//        viewer.core.point_size = 5;
//        viewer.data.set_mesh(V, F);
//        Eigen::MatrixXd SP;
//        DetectingOverhangs overhangs_detector;
//        overhangs_detector.sample(V, F, SP);
//        viewer.data.add_points(SP, Eigen::RowVector3d(1, 0 ,0));
//    }
//    else if (key == '4')
//    {
//        viewer.data.clear();
//        SceneOrganizer organizer;
//        organizer.add_mesh(V, F, Eigen::RowVector3d(1, 1 ,0));
//        RenderingTreeSupport render;
//
//
//        render.add_root(Eigen::Vector3d(0,0,0));
//        render.add_node(Eigen::Vector3d(0,0.1,0));
//        render.add_node(Eigen::Vector3d(-0.05,0.15,0));
//        render.add_node(Eigen::Vector3d(-0.05,0.16,0));
//        render.add_node(Eigen::Vector3d(0.05,0.15,0));
//        render.add_node(Eigen::Vector3d(0.05,0.16,0));
//
//        render.add_edge(0, 1);
//        render.add_edge(1, 2);
//        render.add_edge(2, 3);
//        render.add_edge(1, 4);
//        render.add_edge(4, 5);
//        render.draw(V, F);
//
//        organizer.add_mesh(V, F, Eigen::RowVector3d(1, 0 ,0));
//        Eigen::MatrixXd C;
//        organizer.get_mesh(V, F, C);
//
//        viewer.data.set_face_based(true);
//        viewer.data.set_mesh(V, F);
//        viewer.data.set_colors(C);
//    }
//    else if (key == '5')
//    {
//
//        viewer.data.clear();
//        SceneOrganizer organizer;
//        organizer.add_mesh(V, F, Eigen::RowVector3d(1, 1 ,0));
//
//        GeneratingSupport supporter;
//        Eigen::MatrixXd sV;
//        Eigen::MatrixXi sF;
//        supporter.generate_support(V, F, sV, sF);
//        organizer.add_mesh(sV, sF, Eigen::RowVector3d(1, 0 ,0));
//
//        Eigen::MatrixXd C;
//        organizer.get_mesh(V, F, C);
//        viewer.data.set_face_based(true);
//        viewer.data.set_mesh(V, F);
//        viewer.data.set_colors(C);
//    }
//    else if( key=='6')
//    {
//        loadModel(V, F);
//
//        viewer.data.clear();
//        SceneOrganizer organizer;
//        organizer.add_mesh(V, F, Eigen::RowVector3d(1, 1 ,0));
//
//        GeneratingSupport supporter;
//        Eigen::MatrixXd sV;
//        Eigen::MatrixXi sF;
//        supporter.generate_support(V, F, sV, sF);
//        organizer.add_mesh(sV, sF, Eigen::RowVector3d(1, 0 ,0));
//
//        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(9, 11);
//        GeneratingPlatform platform_builder;
//        platform_builder.draw_platform(V, F, H);
//        organizer.add_mesh(V, F, Eigen::RowVector3d(0, 0 ,1));
//
//        Eigen::MatrixXd C;
//        organizer.get_mesh(V, F, C);
//        viewer.data.set_face_based(true);
//        viewer.data.set_mesh(V, F);
//        viewer.data.set_colors(C);
//    }
//    else if(key == '7')
//    {
//        loadModel(V, F);
//        viewer.data.clear();
//        SceneOrganizer organizer;
//        organizer.add_mesh(V, F, Eigen::RowVector3d(1, 1 ,0));
//
//        GeneratingSupport supporter;
//        supporter.platform_button_ = true;
//        Eigen::MatrixXd sV;
//        Eigen::MatrixXi sF;
//        supporter.generate_support(V, F, sV, sF);
//        organizer.add_mesh(sV, sF, Eigen::RowVector3d(1, 0 ,0));
//
//        Eigen::MatrixXd H = supporter.H_;
//        std::cout << H << std::endl;
//        GeneratingPlatform platform_builder;
//        platform_builder.draw_platform(V, F, H);
//        organizer.add_mesh(V, F, Eigen::RowVector3d(0, 0 ,1));
//
//        Eigen::MatrixXd C;
//        organizer.get_mesh(V, F, C);
//        viewer.data.set_face_based(true);
//        viewer.data.set_mesh(V, F);
//        viewer.data.set_colors(C);
//    }
//    return false;
//}
//
//int main(int argc, char *argv[])
//{
//    void loadModel(Eigen::MatrixXd &V, Eigen::MatrixXi &F);
//    void getBoundingBox(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &V_nox);
//    void setViewer(igl::viewer::Viewer &viewer, Eigen::MatrixXd V, Eigen::MatrixXi &F, Eigen::MatrixXd V_box);
//
//    igl::viewer::Viewer viewer;
//
//    loadModel(V,F);
//
//    viewer.core.show_lines = false;
//
//    getBoundingBox(V, F, V_Box);
//
//    setViewer(viewer, V, F, V_Box);
//
//    // Set Keyboard map
//    viewer.callback_key_down = &key_down;
//
//    // Launch the viewer
//    viewer.launch();
//}
//void slicing()
//{
//    Eigen::MatrixXi tF(std::min(present_layer,(int)F.rows()), 3);
//    for(size_t kd = 0; kd < F.rows() && kd < present_layer; kd++)
//        tF.row(kd) = F.row(kd);
//
//    viewer.data.clear();
//    viewer.data.set_mesh(V, tF);
//}
#include <igl/viewer/Viewer.h>
#include <igl/triangle/triangulate.h>
#include <nanogui/formhelper.h>
#include <nanogui/screen.h>
#include <iostream>
#include "igl/readOBJ.h"

#include "testing_models_path.h"
#include "slice.h"
#include "clipper.hpp"
#include "normalizing_model.h"
#include "scene_organizer.h"
#include "slice_overhang_detector.h"
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
Slice_Overhang_Detector slicer;

void loadModel()
{
    // Load a mesh in OBJ format
    igl::readOFF(TESTING_MODELS_PATH "/arch2.off", V, F);
    NormalizingModel normaler;
    normaler.size_normalize(V);
//        V.resize(4, 3);®
//    F.resize(4, 3);
//    V << 0.0, 0.0, 0.0,
//         1.0, 0.0, 0.0,
//         0.0, 1.0, 0.0,
//         0.0, 0.0, 1.0;
//    F << 0, 2, 1,
//         0, 1, 3,
//         0, 3, 2,
//         1, 2, 3;
//    NormalizingModel normaler;
//    normaler.size_normalize(V);

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
    std::cout << V << F << std::endl;
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
        std::cout <<"id:\t" << id << std::endl;
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
        std::cout <<"id:\t" << id << std::endl;
        slicer.get_intersecting_surface(slices, id, V, F);
        if(V.rows() >= 3)
            organizer.add_mesh(V, F, Eigen::RowVector3d(1, 1, 0));
    }
    Eigen::MatrixXd C;
    organizer.get_mesh(V, F, C);
    std::cout << V.rows() << ", " << F.rows() << std::endl;
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

int main(int argc, char *argv[])
{
    loadModel();
    slicing();
    viewer.callback_init = [&](igl::viewer::Viewer& viewer)
    {
        viewer.ngui->addGroup("Slicing Settings");
        viewer.ngui->addVariable("layer", layer);
        viewer.ngui->addButton("slice", show_slice);
        viewer.ngui->addButton("recovery", recover);
        viewer.ngui->addButton("series slicing", series_slicing);
        viewer.ngui->addButton("overhang slicing", overhang_slicing);
        viewer.ngui->addButton("sampling", sampling);
        viewer.screen->performLayout();
        return false;
    };

    viewer.callback_key_down = &key_down;

    //draw();
    viewer.data.set_mesh(V,F);
    viewer.launch();
}
