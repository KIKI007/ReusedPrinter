
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
    double layout_dz;
    string model_name;
    int layout_opt_type;
    bool layout_output;
    bool metal_pin;
    bool load_platform;
    bool single_layer;
}menu_input;

struct Scene_Data
{
    MatrixXd V; //original model vertices
    MatrixXi F; //original model faces
    MeshSlicerShift slicer;
    vector<Paths> tool_path;
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
    menu_input.model_name = "Arch";
    menu_input.layer = 0;
    menu_input.layer_step = 5;
    menu_input.angle = 0;
    menu_input.layout_dy = 0;
    menu_input.layout_dx = 0;
    menu_input.layout_dz = 0;
    menu_input.layout_output = false;
    menu_input.single_layer = false;
}

void load_model()
{
	void render_mesh(MatrixXd &V, MatrixXi &F, MatrixXd &C);
#ifdef _WIN32
	string file_path = file(TESTING_MODELS_PATH, menu_input.model_name, "", "OBJ");
#elif __APPLE__
    string file_path = file(TESTING_MODELS_PATH, menu_input.model_name, "", "STL");
#endif
    if(loadModel(scene_data.V, scene_data.F, file_path))
    {
        scene_data.slicer.clear();
        scene_data.slicer.set_mesh(scene_data.V, scene_data.F);
        SceneOrganizer organizer;
        organizer.add_mesh(scene_data.V, scene_data.F, RowVector3d(1, 1, 0));
        GeneratingPlatform platform_builder;
        MatrixXd V, plV, C, platform; MatrixXi F, plF;
        Settings settings;
        platform = MatrixXd::Zero(settings.pillar_row, settings.pillar_column);
        platform_builder.draw_platform(plV, plF, platform);
        organizer.add_mesh(plV, plF, RowVector3d(0, 0, 1));
        organizer.get_mesh(V, F, C);
        render_mesh(V, F, C);
        viewer.core.align_camera_center(scene_data.V ,scene_data.F);
    }
    return;
}

void render_mesh(MatrixXd &V, MatrixXi &F)
{
    viewer.data.clear();
    viewer.data.set_mesh(V, F);


    return;
}

void render_mesh(MatrixXd &V, MatrixXi &F, MatrixXd &C)
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
    if(move_xy_opt(scene_data.slicer, V, F, C, menu_input.layout_output, menu_input.model_name))
        render_mesh(V, F, C);
}

void rotate()
{
    MatrixXd V, C;
    MatrixXi F;
    if(rotate_opt(scene_data.slicer, V, F, C, menu_input.layout_output, menu_input.model_name))
        render_mesh(V, F, C);
}

void none_opt()
{
    MatrixXd V, C;
    MatrixXi F;
    if(non_move_support(scene_data.slicer, V, F, C, menu_input.layout_output, menu_input.model_name))
        render_mesh(V, F, C);
}

void test()
{
    test_convex_hull(scene_data.slicer, viewer);
    //test_fermat_curve(viewer);
}

void gcode_viewer()
{
    support_previwer(viewer, scene_data.slicer,
                     menu_input.layer_step,
                     menu_input.metal_pin,
                     scene_data.tool_path,
                     menu_input.load_platform,
                     menu_input.model_name,
                     menu_input.single_layer);
}

void move_model()
{
    MatrixXd V, C;
    MatrixXi F;
    if(model_move(scene_data.slicer, V, F, C, menu_input.layout_dx, menu_input.layout_dy, menu_input.layout_dz))
        render_mesh(V, F, C);
}

void gcode_model_output()
{
    write_gcode_model(scene_data.slicer, menu_input.model_name);
}

void gcode_output()
{
    if(!scene_data.tool_path.empty())
        write_gcode(scene_data.slicer, scene_data.tool_path, menu_input.model_name);
}

void height_map()
{
    draw_height_map(scene_data.slicer, menu_input.model_name);
}

void support_map()
{
    draw_support_map(scene_data.slicer, menu_input.model_name);
}

void load_info()
{
	MatrixXd V, C;
	MatrixXi F;
	if (load_file_move_model(scene_data.slicer, V, F, C, menu_input.model_name))
		render_mesh(V, F, C);
}

int main(int argc, char *argv[])
{
    init();

    viewer.callback_init = [&](igl::viewer::Viewer& viewer)
    {
        viewer.ngui->addGroup("Model Loading");
        viewer.ngui->addVariable("Model Name", menu_input.model_name);
        viewer.ngui->addVariable("Overhang Width", scene_data.slicer.settings.overhang_offset);
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

        viewer.ngui->addWindow(Eigen::Vector2i(220,295),"Mesh Layout Class");
        viewer.ngui->addGroup("Move Model");
        viewer.ngui->addVariable("+X", menu_input.layout_dx);
        viewer.ngui->addVariable("+Y", menu_input.layout_dy);
        viewer.ngui->addVariable("+Z", menu_input.layout_dz);
        viewer.ngui->addButton("move", move_model);
        viewer.ngui->addGroup("Virtual Support");
        viewer.ngui->addVariable("Output", menu_input.layout_output);
		viewer.ngui->addButton("Load Info", load_info);
        viewer.ngui->addButton("None Layout Opt", none_opt);
        viewer.ngui->addButton("XZ Layout Opt", xz_opt);
        viewer.ngui->addButton("Rotate Layout Opt", rotate);
        viewer.ngui->addGroup("Color Map");
        viewer.ngui->addButton("Height Map", height_map);
        viewer.ngui->addButton("Support Map", support_map);

        viewer.ngui->addWindow(Eigen::Vector2i(390,10),"Mesh Support Class");
        viewer.ngui->addGroup("Parameter");
        viewer.ngui->addVariable("Gap width", scene_data.slicer.settings.fermat_cut_width);
        viewer.ngui->addVariable("Center size",  scene_data.slicer.settings.support_center_factor);
        viewer.ngui->addVariable("Group size", scene_data.slicer.settings.group_expand_size);

        viewer.ngui->addGroup("Previewer");
        viewer.ngui->addVariable("Metal Pin", menu_input.metal_pin);
        viewer.ngui->addVariable("Load .plf", menu_input.load_platform);
        viewer.ngui->addVariable("single layer", menu_input.single_layer);
        viewer.ngui->addButton("Gcode Preview", gcode_viewer);

        viewer.ngui->addGroup("Generate Gcode");
        viewer.ngui->addButton("Model Output", gcode_model_output);
        viewer.ngui->addButton("Gcode Output", gcode_output);

        viewer.screen->performLayout();
        return false;
    };

    viewer.launch();
    return 0;
}
