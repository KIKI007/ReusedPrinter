//
// Created by 汪子琦 on 4/23/17.
//

#ifndef EXAMPLE_GCODE_H
#define EXAMPLE_GCODE_H

#include <string>
#include <cstring>
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <assert.h>
#include <fstream>
#include "settings.h"
#include "clipper.hpp"

#define MAX_CHAR_NUM_LINE 1024

const unsigned G_Type = 1;
const unsigned M_Type = 2;
const unsigned S_Type = 4;
const unsigned E_Type = 8;
const unsigned F_Type = 16;
const unsigned X_Type = 32;
const unsigned Y_Type = 64;
const unsigned Z_Type = 128;
const unsigned Comment_Type = 256;

using std::vector;
using ClipperLib::Paths;
using ClipperLib::Path;

class GcodeLine
{
public:

    GcodeLine()
    {
        zero_eps = 1e-6;
        Empty = 0;
        Comment = Comment_Type;
        M_Code = M_Type;
        M_S_Code = M_Type | S_Type;
        G_Code = G_Type;
        G_E_Code = G_Type | E_Type;
        G_E_F_Code = G_Type | E_Type | F_Type;
        G_Z_Code = G_Type | Z_Type;
        G_Z_E_Code = G_Type | Z_Type | E_Type;
        G_Z_F_Code = G_Type | Z_Type | F_Type;
        G_X_Y_F_Code = G_Type | X_Type | Y_Type | F_Type;
        G_X_Y_E_F_Code = G_Type | X_Type | Y_Type | E_Type | F_Type;
        G_X_Y_E_Code = G_Type | X_Type | Y_Type | E_Type;
    }

    void clear()
    {
        type = Empty;
        nG = nM = 0;
        nE = nS = nF = nX = nY = nZ = 0;
        comment = "";
        zero_eps = 1e-6;
    }

public:

    void SplitString(const std::string& s, std::vector<std::string>& v, const std::string& c)
    {
        std::string::size_type pos1, pos2;
        pos2 = s.find(c);
        pos1 = 0;
        while(std::string::npos != pos2)
        {
            v.push_back(s.substr(pos1, pos2-pos1));

            pos1 = pos2 + c.size();
            pos2 = s.find(c, pos1);
        }
        if(pos1 != s.length())
            v.push_back(s.substr(pos1));
    }

    void read(std::string &str)
    {
        if(str.length() == 0)
        {
            type = Empty;
            return;
        }
        else
        {
            if(remove_comment(str)) return;
            std::vector<std::string> split_strs;
            std::string c = " ";
            SplitString(str, split_strs, c);
            std::vector<char> signs(split_strs.size());
            std::vector<double> nums(split_strs.size());
            std::string name;
            for(int id = 0; id < split_strs.size(); id++)
            {
                Char_Double(split_strs[id], signs[id], nums[id]);
                if(isalpha(signs[id]))
                    name.push_back(signs[id]);
                if(signs[id] == 'G')
                {
                    nG = (int) std::round(nums[id]);
                }
                else if(signs[id] == 'M')
                {
                    nM = (int) std::round(nums[id]);
                }
                else if(signs[id] == 'S')
                {
                    nS = nums[id];
                }
                else if(signs[id] == 'E')
                {
                    nE = nums[id];
                }
                else if(signs[id] == 'F')
                {
                    nF = nums[id];
                }
                else if(signs[id] == 'X')
                {
                    nX = nums[id];
                }
                else if(signs[id] == 'Y')
                {
                    nY = nums[id];
                }
                else if(signs[id] == 'Z')
                {
                    nZ = nums[id];
                }
            }
            /*
                M_Code,
                M_S_Code,
                G_Code,
                G_E_Code,
                G_E_F_Code,
                G_Z_Code,
                G_Z_F_Code,
                G_X_Y_F_Code,
                G_X_Y_E_F_Code,
                G_X_Y_E_Code
             */
            if(name == "M")
            {
                type = M_Code;
            }
            else if(name =="MS")
            {
                type = M_S_Code;
            }
            else if(name == "G")
            {
                type = G_Code;
            }
            else if(name == "GE")
            {
                type = G_E_Code;
            }
            else if(name == "GEF")
            {
                type = G_E_F_Code;
            }
            else if(name == "GZ")
            {
                type = G_Z_Code;
            }
            else if(name == "GZF")
            {
                type = G_Z_F_Code;
            }
            else if(name == "GXYF")
            {
                type = G_X_Y_F_Code;
            }
            else if(name == "GXYEF")
            {
                type = G_X_Y_E_F_Code;
            }
            else if(name == "GXYE")
            {
                type = G_X_Y_E_Code;
            }
            else if(name == "GZE")
            {
                type = G_Z_E_Code;
            }
        }
    }

public:

    std::string print()
    {
        std::string str;
        char char_str[MAX_CHAR_NUM_LINE];

        if(type == M_Code)
        {
            sprintf(char_str, "M%d", nM);
            str = std::string(char_str);
        }
        else if(type == M_S_Code)
        {
            sprintf(char_str, "M%d S%g", nM, nS);
            str = std::string(char_str);
        }
        else if(type == G_Code)
        {
            sprintf(char_str, "G%d", nG);
            str = std::string(char_str);
        }
        else if(type == G_E_Code)
        {
            sprintf(char_str, "G%d E%g", nG, nE);
            str = std::string(char_str);
        }
        else if(type == G_E_F_Code)
        {
            sprintf(char_str, "G%d E%.5lf F%.5lf", nG, nE, nF);
            str = std::string(char_str);
        }
        else if(type == G_Z_Code)
        {
            sprintf(char_str, "G%d Z%.3lf", nG, nZ);
            str = std::string(char_str);
        }
        else if(type == G_Z_F_Code)
        {
            sprintf(char_str, "G%d Z%.3lf F%.3lf", nG, nZ, nF);
            str = std::string(char_str);
        }
        else if(type == G_X_Y_F_Code)
        {
            sprintf(char_str, "G%d X%.3lf Y%.3lf F%.3lf", nG, nX, nY, nF);
            str = std::string(char_str);
        }
        else if(type == G_X_Y_E_F_Code)
        {
            sprintf(char_str, "G%d X%.3lf Y%.3lf E%.5lf F%.3lf", nG, nX, nY, nE, nF);
            str = std::string(char_str);
        }
        else if(type == G_X_Y_E_Code)
        {
            sprintf(char_str, "G%d X%.3lf Y%.3lf E%.5lf", nG, nX, nY, nE);
            str = std::string(char_str);
        }
        else if(type == G_Z_E_Code)
        {
            sprintf(char_str, "G%d Z%g E%.1lf", nG, nZ, nE);
            str = std::string(char_str);
        }

        if(comment.length() > 0)
        {
            if(type != Comment) str.push_back(' ');
            str.insert(str.end(), comment.begin(), comment.end());
        }

        return str;
    }

public:

    bool remove_comment(std::string &str)
    {
        int id = 0;
        for(id = 0; id < str.length(); id++)
        {
            if(str[id] == ';')
            {
                comment = str.substr(id);
                str = str.substr(0, id);
                break;
            }
        }
        if(id == 0)
        {
            type = Comment;
            return true;
        }
        else
        {
            return false;
        }
    }

    void Char_Double(std::string str, char &sgn, double &num)
    {
        char *char_str = new char[str.length()];
        std::strcpy(char_str, str.c_str());
        sscanf(char_str, "%c%lf",&sgn, &num);
        return;
    }

public:

    bool is_GZ()
    {
        if((type & G_Type) && (type & Z_Type))
            return true;
        else
            return false;
    }

    bool is_M107()
    {
        if(type == M_Code && nM == 107)
            return true;
        else
            return false;
    }

    bool is_G1XY()
    {
        if(nG != 1) return false;
        if((type & G_Type) && (type & X_Type) && (type & Y_Type))
            return true;
        else
            return false;
    }

    bool is_G92E0()
    {
        if(nG != 92 || std::abs(nE) < zero_eps)return false;
        if((type & G_Type) && (type & E_Type))
            return true;
        else
            return false;
    }

    bool is_G1E()
    {
        if(nG != 1)return false;
        if((type & G_Type) && (type & E_Type))
            return true;
        else
            return false;
    }


    bool is_G1EF()
    {
        if(nG != 1)return false;
        if((type & G_Type) && (type & E_Type) && (type & F_Type))
            return true;
        else
            return false;
    }

public:

    void set_G1_XYF(double X, double Y, double F)
    {
        clear();
        type = G_X_Y_F_Code;
        nG = 1;
        nX = X;
        nY = Y;
        nF = F;
    }

    void set_G1_XYE(double X, double Y, double E)
    {
        clear();
        type = G_X_Y_E_Code;
        nG = 1;
        nX = X;
        nY = Y;
        nE = E;
    }

    void set_G1_XYEF(double X, double Y, double E, double F)
    {
        clear();
        type = G_X_Y_E_F_Code;
        nG = 1;
        nX = X;
        nY = Y;
        nE = E;
        nF = F;
    }

    void set_G92E0()
    {
        clear();
        type = G_E_Code;
        nG = 92;
        nE = 0;
    }

    void set_G1EF(double E, double F)
    {
        clear();
        type = G_E_F_Code;
        nG = 1;
        nE = E;
        nF = F;
        return;
    }

    void set_G1ZF(double Z, double F)
    {
        clear();
        type = G_Z_F_Code;
        nG = 1;
        nZ = Z;
        nF = F;
    }

    void set_G5()
    {
        clear();
        type = G_Code;
        nG = 5;
    }

public:
    unsigned Empty;
    unsigned Comment;
    unsigned M_Code;
    unsigned M_S_Code;
    unsigned G_Code;
    unsigned G_E_Code;
    unsigned G_E_F_Code;
    unsigned G_Z_Code;
    unsigned G_Z_E_Code;
    unsigned G_Z_F_Code;
    unsigned G_X_Y_F_Code;
    unsigned G_X_Y_E_F_Code;
    unsigned G_X_Y_E_Code;

public:

    unsigned type;
    std::string comment;
    int nG;
    int nM;
    double nF;
    double nS;
    double nE;
    double nZ;
    double nX;
    double nY;

    double zero_eps;
};

class Gcode
{
public:

    Gcode(std::string pathname)
    {
        std::ifstream fin;
        fin.open(pathname.c_str());
        std::vector<GcodeLine> lines;
        if(fin.is_open())
        {
            char str_line[MAX_CHAR_NUM_LINE];
            while(fin.getline(str_line, MAX_CHAR_NUM_LINE))
            {
                std::string str(str_line);
                GcodeLine line;
                line.read(str);
                lines.push_back(line);
            }
            layering(lines);
            //remove_G1Z_code_inside_layer();
        }
    }

public:

    void print(std::string pathname)
    {
        std::ofstream fout;
        fout.open(pathname.c_str());
        if(fout.is_open())
        {
            for(int id = 0; id < lines_layer.size(); id++)
            {
                for(int jd = 0; jd < lines_layer[id].size(); jd++)
                {
                    fout << lines_layer[id][jd].print() << std::endl;
                }
            }
        }
        fout.close();
        return;
    }

    void layering(std::vector<GcodeLine> &lines)
    {
        int num_lines = lines.size();
        int id = 0, jd = 0;
        int sta_line = 0, end_line = 0;
        std::vector<GcodeLine> head;
        std::vector<GcodeLine> tail;

        //layer head
        for(id = 0; id < num_lines; id++)
        {
            if(lines[id].is_GZ())
                break;
            else
                head.push_back(lines[id]);
        }
        sta_line = id;

        //layer tail
        for(id = num_lines - 1; id >= 0; id--)
        {
            if(lines[id].is_M107())
            {
                end_line = id - 1;
                break;
            }
        }
        for(id = end_line + 1; id < num_lines; id++)
        {
            tail.push_back(lines[id]);
        }

        //rest
        int total_layer = 0;
        std::vector<int> G1_Z;
        for(id = sta_line; id <= end_line; id++)
        {
            if(lines[id].is_GZ())
            {
                G1_Z.push_back(id);
                if(total_layer < get_layer(lines[id]))
                    total_layer = get_layer(lines[id]);
            }
        }

        // layer end
        std::vector<int> layer_end;
        layer_end.resize(total_layer + 2, 0);
        for(id = 0; id < G1_Z.size(); id++)
        {
            int layer = get_layer(lines[G1_Z[id]]);
            layer_end[layer] = G1_Z[id];
        }
        for(id = 1; id <= total_layer; id++)
        {
            for(jd = layer_end[id] + 1; jd <= end_line; jd++)
            {
                if(lines[jd].is_GZ())
                    break;
            }
            layer_end[id] = jd - 1;
        }
        layer_end[0] = sta_line - 1;
        layer_end[total_layer + 1] = num_lines;

        // layer
        lines_layer.resize(layer_end.size());
        for(id = 0; id < layer_end.size(); id++)
        {
            //std::cout << layer_end[id] + 1 << std::endl;
            if(id == 0) lines_layer[id] = head;
            else if(id == total_layer + 1) lines_layer[id] = tail;
            else
            {
                for(jd = layer_end[id - 1] + 1; jd <= layer_end[id]; jd++)
                    lines_layer[id].push_back(lines[jd]);
            }
        }
    }

    int get_layer(GcodeLine &code)
    {
        return std::round(code.nZ / settings.layer_height);
    }

    void get_minXY()
    {
        min_X = settings.MAX_DOUBLE;
        min_Y = settings.MAX_DOUBLE;
        double X = 0, Y = 0;
        for(int layer = 1; layer < lines_layer.size() - 1; layer++) {
            for(int id = 0; id < lines_layer[layer].size(); id++)
            if (lines_layer[layer][id].is_G1XY() && lines_layer[layer][id].is_G1E())
            {
                min_X = std::min(min_X, X);
                min_X = std::min(min_X, lines_layer[layer][id].nX);
                min_Y = std::min(min_Y, Y);
                min_Y = std::min(min_Y, lines_layer[layer][id].nY);
            }
        }
        min_Y -= settings.extrusion_width / 2;
        min_X -= settings.extrusion_width / 2;
    }

    void move_XY(double dX, double dY)
    {
        for(int layer = 1; layer < lines_layer.size() - 1; layer++) {
                for(int id = 0; id < lines_layer[layer].size(); id++)
                {
                    if(lines_layer[layer][id].is_G1XY())
                    {
                        lines_layer[layer][id].nX += dX;
                        lines_layer[layer][id].nY += dY;
                    }
                }
        }
        return;
    }


    void remove_G1Z_code_inside_layer()
    {
        for(int id = 1; id < lines_layer.size() - 1; id++)
        {
            if(lines_layer[id].empty()) continue;
            std::vector<GcodeLine> new_layer;
            new_layer.push_back(lines_layer[id][0]);
            for(int jd = 1; jd < lines_layer[id].size(); jd++)
            {
                if(!lines_layer[id][jd].is_GZ())
                    new_layer.push_back(lines_layer[id][jd]);
            }
            lines_layer[id] = new_layer;
        }
    }

    void add_support(vector<Paths> &fermat_curves)
    {
        get_dE_D_dL();                      //get the ratio between dE and dL
        double X = 0, Y = 0, F = 0, E = 0, dL;  //present status value
        GcodeLine code;                     //temporary gcode term
        bool first_move_touch, new_speed;
        for(int layer = 1; layer < fermat_curves.size(); layer++)
        {
            first_move_touch = true;        //to address problem of the filament real in and off
            for(int id = 0; id < lines_layer[layer].size(); id++)
            {
                code = lines_layer[layer][id];

                //The first G1E code is to real in the filament
                //Once we add the support, this code no longer meaning real in 7mm filament
                //We need to find the previous layer's E and minus 7
                if(code.is_G1E())
                {
                    if(first_move_touch && layer > 1)
                    {
                        lines_layer[layer][id].nE = E - 7;
                        first_move_touch = false;
                    }
                    E = code.nE;
                }

                //To locate the position of the printing nozzel and find the travel speed;
                if(code.is_G1XY()) {
                    X = code.nX;
                    Y = code.nY;
                    if(code.is_G1EF())
                    {
                        F = code.nF;
                    }
                }

            }

            vector<GcodeLine> support;


            //real in 7mm material
            code.set_G1EF(E - 7, 4800);
            support.push_back(code);

            //set E to the orgin;
            code.set_G92E0();
            support.push_back(code);

            //move Z a layer height to avoid collision
            code.set_G1ZF((layer + 1) * settings.layer_height, 1800);
            support.push_back(code);

            int curve_id = 0;
            new_speed = false;
            first_move_touch = true;

            for(curve_id = 0; curve_id != fermat_curves[layer].size(); curve_id++)
            {
                Path path = fermat_curves[layer][curve_id];
                for(int id = 1; id < path.size(); id++)
                {
                    double X1 = settings.int2mm(path[id - 1].X) + settings.platform_zero_x;
                    double Y1 = -settings.int2mm(path[id - 1].Y) + settings.platform_zero_y;
                    double X2 = settings.int2mm(path[id].X) + settings.platform_zero_x;
                    double Y2 = -settings.int2mm(path[id].Y) + settings.platform_zero_y;

                    //if(X1 != X, Y1 != Y), we have to move the printing nozzel
                    if(std::abs(X1 - X) > settings.ZERO_EPS || std::abs(Y1 - Y) > settings.ZERO_EPS)
                    {
                        //we first move the nozzel
                        code.set_G1_XYF(X1, Y1, settings.nF_moving);
                        support.push_back(code);

                        //if this move from mesh to the support
                        //we have to reel off 7.1mm material to make sure the new filament can be attached to the support
                        if(first_move_touch)
                        {
                            //move down nozzel
                            code.set_G1ZF(layer * settings.layer_height, 1800);
                            support.push_back(code);

                            //reel off material
                            code.set_G1EF(7.1, 4800);
                            support.push_back(code);

                            first_move_touch = false;
                            E = 7.1;
                        }

                        new_speed = true;
                    }


                    dL = std::sqrt((X2 - X1) * (X2 - X1) + (Y2 - Y1) * (Y2 - Y1)); //compute the travel distance
                    E += dL * dEdL[layer];

                    if(new_speed) code.set_G1_XYEF(X2, Y2, E, F); //the new_speed means to set the travel speed as same as the mesh
                    else code.set_G1_XYE(X2, Y2, E);
                    support.push_back(code);

                    X = X2; Y = Y2;
                    new_speed = false;
                }
            }
            lines_layer[layer].insert(lines_layer[layer].end(), support.begin(), support.end());

            code.set_G5();
            lines_layer[layer].push_back(code);
        }
    }

    void get_dE_D_dL()
    {
        dEdL.resize(lines_layer.size());
        double X = 0, Y = 0, E = 0;
        for(int layer = 1; layer < lines_layer.size() - 1; layer++)
        {
            std::vector<double> dLs;
            std::vector<double> dEs;
            for(int id = 0; id < lines_layer[layer].size(); id++)
            {
                GcodeLine code = lines_layer[layer][id];
                double X1 = 0, Y1 = 0, E1 = 0;
                if(code.is_G92E0())
                {
                    E = 0;
                }
                else if(code.is_G1XY())
                {
                    if(code.is_G1E())
                    {
                        X1 = code.nX;
                        Y1 = code.nY;
                        E1 = code.nE;
                        dLs.push_back(std::sqrt((X - X1) * (X - X1) + (Y - Y1) * (Y - Y1) ));
                        dEs.push_back(E1 - E);
                        E = E1;
                    }
                    X = code.nX;
                    Y = code.nY;
                }
                else if(code.is_G1E())
                {
                    E = code.nE;
                }
            }

            double average_dE_D_dL = 0;
            for(int id = 0; id < dLs.size(); id++) {
                //std::cout << dEs[id] << " " << dLs[id] << " " << dEs[id] / dLs[id] << std::endl;
                average_dE_D_dL += dEs[id] / dLs[id];
            }
            if(!dLs.empty())average_dE_D_dL /= dLs.size();
            dEdL[layer] = average_dE_D_dL;
        }
        return;
    }

public:

    std::vector<std::vector<GcodeLine>> lines_layer;

    Settings settings;

    std::vector<double> dEdL;

    double min_X;

    double min_Y;
};

#endif //EXAMPLE_GCODE_H
