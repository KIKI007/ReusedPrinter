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
//using namespace ClipperLib;
//// Input polygon
//Eigen::MatrixXd V;
//Eigen::MatrixXi E;
//Eigen::MatrixXd H;
//
//// Triangulated interior
//Eigen::MatrixXd V2;
//Eigen::MatrixXi F2;
//void draw()
//{
//    Paths subj(2), solution;
//    subj[0] <<
//            IntPoint(-100,-100) << IntPoint(100,-100) <<
//            IntPoint(100,100) << IntPoint(-100,100);
//    subj[1] <<
//            IntPoint(-50,-50) << IntPoint(0,50) << IntPoint(50,-50);
//    Clipper c;
//    c.AddPaths(subj, ptSubject, true);
//    c.Execute(ctIntersection, solution, pftEvenOdd, pftEvenOdd);
////    V.resize(7, 2);
////    E.resize(7, 2);
////    H.resize(1, 2);
////    size_t jd = 0;
////    for(size_t kd = 0; kd < solution.size(); kd++)
////    {
////        for(size_t id = 0;id < solution[kd].size(); id++)
////        {
////            IntPoint point = solution[kd][id];
////            V.row(jd++) = Eigen::RowVector2d(point.X,point.Y);
////            E.row(jd++) = Eigen::RowVector2i(point.X,point.Y);
////        }
////    }
////    H << 0 , 0;
////    igl::triangle::triangulate(V,E,H,"a0.005q",V2,F2);
////    igl::viewer::Viewer viewer;
////    viewer.data.set_mesh(V2,F2);
//    V.resize(8,2);
//    E.resize(8,2);
//    H.resize(1,2);
//
//    V << -1,-1, 1,-1, 1,1, -1, 1,
//            -2,-2, 2,-2, 2,2, -2, 2;
//
//    E << 0,1, 1,2, 2,3, 3,0,
//            4,5, 5,6, 6,7, 7,4;
//
//    H << 0,0;
//
//    // Triangulate the interior
//    igl::triangle::triangulate(V,E,H,"a0.005q",V2,F2);
//
//    // Plot the generated mesh
//    igl::viewer::Viewer viewer;
//    viewer.data.set_mesh(V2,F2);
//    viewer.launch();
//}
//
//int main(int argc, char *argv[])
//{
//    Settings settings;
//
//    // Load a mesh in OFF format
//    //igl::readOFF(TESTING_MODELS_PATH "/bunny.off", V, F);
//    //present_layer = F.rows();
//
//    // Extend viewer menu
////    viewer.callback_init = [&](igl::viewer::Viewer& viewer)
////    {
////        // Add an additional menu window
////        viewer.ngui->addWindow(Eigen::Vector2i(220,10),"Print Settings");
////
////        // Add new group
////        viewer.ngui->addGroup("Slicing");
////
////        // Expose variable directly ...
////        viewer.ngui->addVariable("Layer Height",settings.layer_height);
////
////        viewer.ngui->addVariable("Present Layer", present_layer);
////
////        viewer.ngui->addButton("Slicing",slicing);
////
////        // Generate menu
////        viewer.screen->performLayout();
////
////        return false;
////    };
//    draw();
//    viewer.launch();
//}

#include <igl/viewer/Viewer.h>
#include <igl/triangle/triangulate.h>
#include <Eigen/Core>

// Input polygon
Eigen::MatrixXd V;
Eigen::MatrixXi E;
Eigen::MatrixXd H;

// Triangulated interior
Eigen::MatrixXd V2;
Eigen::MatrixXi F2;

int main(int argc, char *argv[])
{
    using namespace Eigen;
    using namespace std;

    // Create the boundary of a square
    V.resize(8,2);
    E.resize(8,2);
    H.resize(1,2);

    V << -1,-1, 1,-1, 1,1, -1, 1,
            -2,-2, 2,-2, 2,2, -2, 2;

    E << 0,1, 1,2, 2,3, 3,0,
            4,5, 5,6, 6,7, 7,4;

    H << 0,0;

    // Triangulate the interior
    igl::triangle::triangulate<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd, Eigen::MatrixXd,Eigen::MatrixXi>(V,E,H,"a0.005q",V2,F2);

    // Plot the generated mesh
    igl::viewer::Viewer viewer;
    viewer.data.set_mesh(V2,F2);
    viewer.launch();
}