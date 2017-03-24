//
// Created by 汪子琦 on 3/19/17.
//

#ifndef SUPPORTER_SETTINGS_H
#define SUPPORTER_SETTINGS_H
#include <limits>

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
        maximum_model_heigh = 50; //mm

        //support generation
        face_overhang_angle = 50 * PI / 180;
        expected_sample_num = 2000;
        printing_max_angle = 45 * PI / 180;

        //platform
        pad_size = 12.7; //mm 0.5inch
        pad_thickness = 1.4224; //mm 0.056inch
        pillar_length = 304.8; //mm 12inch
        pillar_radius = 3.175; //mm 0.125 inch
    }

public:

    inline long long mm2int_Even(double mm)
    {
        long long INT = (long long)(mm / UNIT);
        if(INT % 2 != 0) INT -= 1;
        return INT;
    }

    inline long long mm2int_Odd(double mm)
    {
        long long INT = (long long)(mm / UNIT);
        if(INT % 2 != 1) INT -= 1;
        return INT;
    }

    inline long long mm2int(double mm)
    {
        return (long long)(mm / UNIT);
    }

    inline double int2mm(long long INT)
    {
        return INT * UNIT;
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

    //rendering
    double maximum_model_heigh; //to shrink the model inside the rendering space

    //support generation
    double face_overhang_angle;
    double expected_sample_num;
    double printing_max_angle;
};

#endif //SUPPORTER_SETTINGS_H
