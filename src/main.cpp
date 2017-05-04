
//Menu
#include <nanogui/formhelper.h>
#include <nanogui/screen.h>
#include "menu_function.h"

//GUI
igl::viewer::Viewer viewer;

struct MenuInput
{
    int layer;
    int layer_step;
    double angle;
    double layout_dx;
    double layout_dy;
    string model_name;
    int layout_opt_type;
}menu_input;

struct Scene_Data
{
    MatrixXd V; //original model vertices
    MatrixXi F; //original model faces
    MeshSlicerShift slicer;
}scene_data;

struct Scene_Color
{
public:
    Scene_Color()
    {
        yellow = RowVector3d (1, 1, 0);
    }

    RowVector3d yellow;
}scene_color;


void init()
{
    menu_input.model_name = "arch";
    menu_input.layer = 0;
    menu_input.layer_step = 5;
    menu_input.angle = 0;
    menu_input.layout_dy = 0;
    menu_input.layout_dx = 0;
}

void load_model()
{
    void render_mesh(MatrixXd &V, MatrixXi &F);
    string file_path = file(TESTING_MODELS_PATH, menu_input.model_name, "stl");
    loadModel(scene_data.V, scene_data.F, file_path);

    scene_data.slicer.clear();
    scene_data.slicer.set_mesh(scene_data.V, scene_data.F);

    render_mesh(scene_data.V, scene_data.F);
}

void render_mesh(MatrixXd &V, MatrixXi &F)
{
    viewer.data.clear();
    viewer.data.set_mesh(V, F);
    return;
}

void render_mesh(MatrixXd &V, MatrixXi &F, MatrixXd C)
{
    viewer.data.clear();
    viewer.data.set_mesh(V, F);
    viewer.data.set_colors(C);
    return;
}

void single_layer_slicing()
{
    MatrixXd V;
    MatrixXi F;
    if(slicing_for_one_layer(scene_data.slicer, menu_input.layer, V, F))
        render_mesh(V, F);
}

void multi_layers_slicing()
{
    MatrixXd V;
    MatrixXi F;
    if(slicing_in_step(scene_data.slicer, menu_input.layer_step, V, F))
        render_mesh(V, F);
}

void bottom_half()
{
    MatrixXd V;
    MatrixXi F;
    if(slicing_bottom_half(scene_data.slicer, V, F))
        render_mesh(V, F);
}

void upper_half()
{
    MatrixXd V;
    MatrixXi F;
    if(slicing_upper_half(scene_data.slicer, V, F))
        render_mesh(V, F);
}

void xz_opt()
{
    MatrixXd V, C;
    MatrixXi F;
    if(move_xy_opt(scene_data.slicer, V, F, C))
        render_mesh(V, F, C);
}

void rotate()
{
    MatrixXd V, C;
    MatrixXi F;
    if(rotate_opt(scene_data.slicer, V, F, C))
        render_mesh(V, F, C);
}

void test()
{
    test_fermat_curve(viewer);
}

int main(int argc, char *argv[])
{
    init();

    viewer.callback_init = [&](igl::viewer::Viewer& viewer)
    {
        viewer.ngui->addGroup("Model Loading");
        viewer.ngui->addVariable("Model Name", menu_input.model_name);
        viewer.ngui->addButton("Load Model", load_model);

        viewer.ngui->addWindow(Eigen::Vector2i(220,10),"Mesh Slicer Class");
        viewer.ngui->addGroup("Basic Slicing");
        viewer.ngui->addVariable("layer", menu_input.layer);
        viewer.ngui->addButton("Single Layer", single_layer_slicing);
        viewer.ngui->addVariable("layer step", menu_input.layer_step);
        viewer.ngui->addButton("Multiple Layer", multi_layers_slicing);
        viewer.ngui->addGroup("Overhang");
        viewer.ngui->addButton("Bottom half", bottom_half);
        viewer.ngui->addButton("Upper half", upper_half);

        viewer.ngui->addWindow(Eigen::Vector2i(220,10),"Mesh Layout Class");
        viewer.ngui->addGroup("Layout Opt");
        viewer.ngui->addButton("XZ Opt", xz_opt);
        viewer.ngui->addButton("Rotate Opt", rotate);

        viewer.ngui->addWindow(Eigen::Vector2i(220,10),"Test");
        viewer.ngui->addGroup("Test");
        viewer.ngui->addButton("Test", test);
        viewer.screen->performLayout();
        return false;
    };

    viewer.launch();
    return 0;
}
