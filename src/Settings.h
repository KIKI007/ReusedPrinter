//
// Created by 汪子琦 on 3/19/17.
//

#ifndef SUPPORTER_SETTINGS_H
#define SUPPORTER_SETTINGS_H
#include <limits>
#include <string>
#include <ctime>
#include <cstdio>

class Settings
{
public:
    Settings()
    {
        //global
        PI = 3.1415926;
        MAX_DOUBLE = std::numeric_limits<double>::max();
        MIN_DOUBLE = std::numeric_limits<double>::min();
        ZERO_EPS = 1e-5;
        UNIT = 1e-3;
        circle_vertices_num = 32;

        //slice
        layer_height = 0.2;

        //rendering
        maximum_model_heigh = 40; //mm
        blender_scale_factor = 0.03; //mm

        //support generation
        overhang_offset = 0.2; //mm
        extrusion_width = 0.5; //mm
        support_center_factor = 4; //mm
        fermat_cut_width = 0.5; //mm
        sew_pin_width_factor = 2;

        //platform
        pad_size = 12.7; //mm 0.5inch
        pad_thickness = 1.4224; //mm 0.056inch
        pillar_length = 4 * 25.4; //mm 12inch
        pillar_radius = 3.175; //mm 0.125 inch
        pillar_standard_height = 12.7; //mm 0.5inch
        pillar_row = 9;
        pillar_column = 11;

        //layout
        xy_sample_num_each_pin = 50;
        sample_width = mm2int(pad_size / xy_sample_num_each_pin);
        maximum_height_map = pillar_standard_height * 100;
        angle_sample_num = 360;
        angle_step = PI / angle_sample_num * 2;
        edge_region_num = xy_sample_num_each_pin / 10;

        //gcode
        nF_printing_rest = 1800;
        nF_moving = 1800;
        nF_printing_first = 540;
        nF_reversing = 4800;
        platform_zero_x = -75 + 3.5 - 0.2 + 0.8;
        platform_zero_y = 60 - 2 + 1;

        group_expand_size = 2;
    }

public:

    inline long long mm2int_Even(double mm)
    {
        long long INT = (long long)(round(mm / UNIT));
        INT *= 2;
        return INT;
    }

    inline long long mm2int_Odd(double mm)
    {
        long long INT = (long long)(round(mm / UNIT));
        INT = INT * 2 + 1;
        return INT;
    }

    inline long long mm2int(double mm)
    {
        return (long long)(round(mm / UNIT * 2));
    }

    inline double int2mm(long long INT)
    {
        return (double)INT * UNIT / 2;
    }

    inline double int2mm(double INT)
    {
        return INT * UNIT / 2;
    }

    inline double int2mm(int INT)
    {
        return (double)INT * UNIT / 2;
    }

public:

    int pin_x(int ir){return ir * sample_width;}

    int pin_y(int ic){return ic * sample_width;}

    int pin_center_x(int ir) {return ir * sample_width + sample_width / 2;}

    int pin_center_y(int ic) {return ic * sample_width + sample_width / 2;}

public:

    void print_N()
    {
//        std::cout << std::endl;
    }


    void print_TsN(std::string script, int tab = 0)
    {
//       while(tab --)std::cout << "\t";
//        std::cout << script << std::endl;
    }

    void print_Ts(std::string script, int tab = 0)
    {
//        while(tab --)std::cout << "\t";
//        std::cout << script;
    }


public:

    void tic(std::string script)
    {
        std::cout<< script << "...";
        start = std::clock();
    }

    void toc()
    {
        double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        std::cout<< duration << "s" <<'\n';
    }

public:
    //global
    double PI;
    double MAX_DOUBLE;
    double MIN_DOUBLE;
    double ZERO_EPS;
    double UNIT;
    int circle_vertices_num;

    //slice
    double layer_height;

    //platform
    double pad_size;
    double pad_thickness;
    double pillar_radius;
    double pillar_length;
    double pillar_standard_height;
    int pillar_row;
    int pillar_column;
    double platform_zero_x;
    double platform_zero_y;

    //rendering
    double maximum_model_heigh; //to shrink the model inside the rendering space
    double blender_scale_factor;

    //support generation
    double overhang_offset;
    double extrusion_width;
    double fermat_cut_width;
    int support_center_factor;
    int group_expand_size;
    int sew_pin_width_factor;

    //layout
    int xy_sample_num_each_pin;
    int sample_width;
    int angle_sample_num;
    double angle_step;
    double maximum_height_map;
    int edge_region_num;

    //gcode
    double nF_printing_first;
    double nF_printing_rest;
    double nF_moving;
    double nF_reversing;


public:
    char tmp_str[1024];

    std::clock_t start;
};




#endif //SUPPORTER_SETTINGS_H
