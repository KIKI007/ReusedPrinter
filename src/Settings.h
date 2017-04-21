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

        //sampling
        sample_distance = 2;

        //rendering
        maximum_model_heigh = 40; //mm

        //support generation
        overhang_offset = 0.2; //mm
        support_width = 0.5; //mm
        support_center_area = 0.25; //mm
        fermat_cut_width = 0.5; //mm

        face_overhang_angle = 50 * PI / 180;
        expected_sample_num = 2000;
        printing_max_angle = 45 * PI / 180;

        //platform
        pad_size = 12.7; //mm 0.5inch
        pad_thickness = 1.4224; //mm 0.056inch
        pillar_length = 304.8; //mm 12inch
        pillar_radius = 3.175; //mm 0.125 inch
        pillar_standard_height = 12.7; //mm 0.5inch
        pillar_row = 9;
        pillar_column = 11;

        //layout
        xy_sample_num_each_pin = 1 << 4;
        maximum_height_map = pillar_standard_height * 100;
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

    //sampling
    double sample_distance;

    //rendering
    double maximum_model_heigh; //to shrink the model inside the rendering space

    //support generation
    double overhang_offset;
    double support_width;
    double support_center_area;
    double fermat_cut_width;

    double face_overhang_angle;
    double expected_sample_num;
    double printing_max_angle;

    //layout
    int xy_sample_num_each_pin;
    double maximum_height_map;

public:
    char tmp_str[1024];

    std::clock_t start;
};




#endif //SUPPORTER_SETTINGS_H
