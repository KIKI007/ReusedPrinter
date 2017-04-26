//
// Created by 汪子琦 on 4/20/17.
//

#ifndef EXAMPLE_MESH_LAYOUT_H
#define EXAMPLE_MESH_LAYOUT_H

#include "mesh_slicer.h"
#include <cmath>
#include <queue>
#include <Eigen/Dense>



typedef struct tagHeightMapNode
{
    int R;
    int C;
    int XYsize;
    int layer;
    bool new_layer;
    ClipperLib::Paths polys;
}HeightMapNode;

class MeshLayout
{
public:
    MeshLayout()
    {
        clear();
    }

    MeshLayout(Settings s)
    {
        settings = s;
    }

public:

    void clear()
    {
        red_map_.setZero();
        height_map_.setZero();
        slicer_.clear();
    }

    void set_slicer(MeshSlicer &slicer)
    {
        clear();
        slicer_ = slicer;
        settings = slicer.return_settings();
    }

public:

    void layer_projecting(Eigen::MatrixXi &map, std::vector<ClipperLib::Paths> &slices);

    void get_height_map(Eigen::MatrixXd &height_map);

    void get_red_map(Eigen::MatrixXi &red_map);

    void xy_layout(double &dx, double &dz, Eigen::MatrixXd &platform);

private:

    double min(double a, double b)
    {
        if(a > b) return b;
        else return a;
    }

public:

    MeshSlicer slicer_;

    Eigen::MatrixXd height_map_;

    Eigen::MatrixXi red_map_;

    Settings settings;
};


void MeshLayout::layer_projecting(Eigen::MatrixXi &map, std::vector<ClipperLib::Paths> &slices)
{
    std::queue<HeightMapNode> Queue;
    map = Eigen::MatrixXi::Zero(settings.pillar_row  * settings.xy_sample_num_each_pin,
                                settings.pillar_column * settings.xy_sample_num_each_pin);

    for(int row_id = 0; row_id < settings.pillar_row; row_id++)
    {
        for(int col_id = 0; col_id < settings.pillar_column; col_id++)
        {
            HeightMapNode node;
            node.R = row_id * settings.xy_sample_num_each_pin;
            node.C = col_id * settings.xy_sample_num_each_pin;
            node.XYsize = settings.xy_sample_num_each_pin;
            node.layer = 0;
            node.new_layer = true;
            node.polys.clear();
            Queue.push(node);
        }
    }

    while(!Queue.empty()) {
        HeightMapNode u = Queue.front();
        Queue.pop();

        if (u.layer >= slices.size()) {
            map.block(u.R, u.C, u.XYsize, u.XYsize) = Eigen::MatrixXi::Ones(u.XYsize, u.XYsize) * -1;
            continue;
        }

        ClipperLib::Clipper clipper;
        if (u.new_layer) {
            clipper.AddPaths(slices[u.layer], ClipperLib::ptClip, true);
        } else {
            clipper.AddPaths(u.polys, ClipperLib::ptClip, true);
        }

        ClipperLib::Path square;
        int L =  settings.mm2int(settings.pad_size / settings.xy_sample_num_each_pin * u.C);
        int R =  settings.mm2int(settings.pad_size / settings.xy_sample_num_each_pin * (u.C + u.XYsize));
        int T =  settings.mm2int(settings.pad_size / settings.xy_sample_num_each_pin * u.R);
        int B =  settings.mm2int(settings.pad_size / settings.xy_sample_num_each_pin * (u.R + u.XYsize));
        square.push_back(ClipperLib::IntPoint(L, T));
        square.push_back(ClipperLib::IntPoint(R, T));
        square.push_back(ClipperLib::IntPoint(R, B));
        square.push_back(ClipperLib::IntPoint(L, B));

        clipper.AddPath(square, ClipperLib::ptSubject, true);

        //computing the intersection
        ClipperLib::Paths intersection;
        clipper.Execute(ClipperLib::ctIntersection, intersection, ClipperLib::pftPositive, ClipperLib::pftPositive);
        ClipperLib::SimplifyPolygons(intersection);

        if (intersection.empty())
        {
            u.layer++;
            u.new_layer = true;
            Queue.push(u);
        }
        else
        {
            //check the collision
            double area_insec = 0;
            for (int id = 0; id < intersection.size(); id++)
                area_insec += ClipperLib::Area(intersection[id]);
            area_insec = std::abs(area_insec);

            double area_squre = std::abs(ClipperLib::Area(square));
            if (std::abs(area_insec - area_squre) < settings.ZERO_EPS || (u.XYsize == 1))
            {
                map.block(u.R, u.C, u.XYsize, u.XYsize) = Eigen::MatrixXi::Ones(u.XYsize, u.XYsize) * u.layer;
            }
            else
            {
                for (int ir = 0; ir < 2; ir++) {
                    for (int ic = 0; ic < 2; ic++) {
                        HeightMapNode v = u;
                        v.R += ir * u.XYsize / 2;
                        v.C += ic * u.XYsize / 2;
                        v.XYsize /= 2;
                        v.polys = intersection;
                        v.new_layer = false;
                        Queue.push(v);
                    }
                }
            }
        }
    }
}

void MeshLayout::get_height_map(Eigen::MatrixXd &height_map)
{
    if(!height_map_.isZero())
    {
        height_map = height_map_;
        return;
    }

    //settings.tic("read");
    std::vector<ClipperLib::Paths> stack_slices, slices;
    slicer_.get_slices(slices);

    //settings.toc();
    Eigen::MatrixXi map;
    // height map
    //settings.tic("project");
    layer_projecting(map, slices);
    //settings.toc();
    height_map.resize(map.rows(), map.cols());

    for(int ir = 0; ir < map.rows(); ir++)
    {
        for(int ic = 0; ic < map.cols(); ic++)
        {
            if(map(ir, ic) >= 0)
            {
                int num_standard_pin = map(ir, ic) * settings.layer_height / settings.pillar_standard_height;
                height_map(ir, ic) =  num_standard_pin * settings.pillar_standard_height;
            }
            else
            {
                height_map(ir, ic)  = settings.maximum_height_map;
            }
        }
    }
    height_map_ = height_map;

    return;
}

void MeshLayout::get_red_map(Eigen::MatrixXi &red_map)
{

    if(!red_map_.isZero())
    {
        red_map = red_map_;
        return;
    }

    //settings.tic("get overhang slices");
    std::vector<ClipperLib::Paths> bottom_half;
    slicer_.get_bottom_half(bottom_half);
    //settings.toc();

    //settings.tic("assemble");
    ClipperLib::Paths downward = bottom_half[0];
    ClipperLib::Clipper clipper;

    for(int layer = 1; layer < bottom_half.size(); layer++)
    {
        if(!bottom_half[layer].empty())
        {
            clipper.AddPaths(bottom_half[layer], ClipperLib::ptSubject, true);
        }
    }
    clipper.Execute(ClipperLib::ctUnion, downward, ClipperLib::pftPositive, ClipperLib::pftPositive);
    //settings.toc();

    Eigen::MatrixXi map;
    std::vector<ClipperLib::Paths> slices;
    slices.push_back(downward);
    layer_projecting(map, slices);
    red_map.resize(map.rows(), map.cols());

    for(int ir = 0; ir < map.rows(); ir++)
    {
        for(int ic = 0; ic < map.cols(); ic++)
        {
            if(map(ir, ic) >= 0)
            {
                red_map(ir, ic) = 1;
            }
            else
            {
                red_map(ir, ic) = 0;
            }
        }
    }
    //settings.toc();

    red_map_= red_map;

    return;
}

void MeshLayout::xy_layout(double &dx, double &dz, Eigen::MatrixXd &platform)
{
    assert(!slicer_.empty());

    dx = 0;
    dz = 0;

    Eigen::MatrixXd height_map;
    Eigen::MatrixXi red_map;

    //settings.tic("height map");
    get_height_map(height_map);
    //settings.toc();

    //settings.tic("red map");
    get_red_map(red_map);
    //settings.toc();

    //settings.tic("layout optimization");
    //Red Anchor
    Eigen::MatrixXi red_anchor;
    red_anchor = Eigen::MatrixXi::Zero(red_map.rows(), red_map.cols());

    for(int id = 0; id < red_map.rows(); id++)
    {
        for(int jd = 0; jd < red_map.cols(); jd++)
        {
            int L = jd > 0 ? red_anchor(id, jd - 1) : 0;
            int T = id > 0 ? red_anchor(id - 1, jd) : 0;
            int LT = id > 0 && jd > 0 ? red_anchor(id - 1, jd - 1) : 0;
            red_anchor(id, jd) = L + T - LT + red_map(id, jd);
        }
    }

    //minimum
    int n = settings.pillar_row * settings.xy_sample_num_each_pin;
    int m = settings.pillar_column * settings.xy_sample_num_each_pin;
    int kn = std::log(n) + 1;
    int km = std::log(m) + 1;

    double ****minimum;
    minimum = new double ***[kn];
    for(int jr = 0; jr < kn; jr++)
        minimum[jr] = new double **[n];

    for(int jr = 0; jr < kn; jr++)
    {
        for(int ir = 0; ir < n; ir++)
        {
            minimum[jr][ir] = new double *[km];
        }
    }

    for(int jr = 0; jr < kn; jr++)
    {
        for(int ir = 0; ir < n; ir++)
        {
            for(int jc = 0; jc < km; jc++)
                minimum[jr][ir][jc] = new double [m];
        }
    }

    for(int ir = 0; ir < n; ir++)
    {
        for (int ic = 0; ic < m; ic++)
        {
            minimum[0][ir][0][ic] = height_map(ir, ic);
        }
    }

    for(int jr = 0; jr < kn; jr ++)
    {
        for(int jc = 0; jc < km; jc++)
        {
            if(jc + jr == 0) continue;
            for(int ir = 0; ir + (1 << jr) - 1 < n; ir++)
            {
                for(int ic = 0; ic + (1 << jc) - 1 < m; ic++)
                {
                    if(jr)
                    {
                        minimum[jr][ir][jc][ic] = min(minimum[jr - 1][ir]                  [jc][ic],
                                                      minimum[jr - 1][ir + (1 << (jr - 1))][jc][ic]);
                    }
                    else
                    {
                        minimum[jr][ir][jc][ic] = min(minimum[jr][ir][jc - 1][ic],
                                                      minimum[jr][ir][jc - 1][ic + (1 << (jc - 1))]);
                    }
                }
            }
        }
    }


    //optimization
    double opt_value = 0;
    double opt_r = 0;
    double opt_c = 0;
    int num_pin = 0;
    for(int dr = 0; dr < settings.xy_sample_num_each_pin; dr++)
    {
        for(int dc = 0; dc < settings.xy_sample_num_each_pin; dc++)
        {
            double opt_tmp = 0;
            int num_pin_tmp = 0;
            Eigen::MatrixXd platform_tmp = Eigen::MatrixXd::Zero(settings.pillar_row, settings.pillar_column);

            for(int ir = 0; ir < settings.pillar_row; ir++)
            {
                for(int ic = 0; ic < settings.pillar_column; ic ++)
                {
                    int r1 = ir * settings.xy_sample_num_each_pin + dr ;
                    int r2 = (ir + 1) * settings.xy_sample_num_each_pin + dr - 1;
                    int c1 = ic * settings.xy_sample_num_each_pin + dc ;
                    int c2 = (ic + 1) * settings.xy_sample_num_each_pin + dc - 1;

                    if(c2 >= m || r2 >= n) continue;

                    //get red anchor
                    int L =  r1 > 0 ?           red_anchor(r1 - 1, c2) : 0;
                    int T =  c1 > 0 ?           red_anchor(r2,     c1 - 1) : 0;
                    int LT = r1 > 0 && c1 > 0 ? red_anchor(r1 - 1, c1 - 1) : 0;
                    int red_num =   red_anchor(r2, c2) - L - T + LT;

                    //get minimum
                    int kr = std::log2(r2 - r1 + 1);
                    int kc = std::log2(c2 - c1 + 1);
                    double minimum_1 = min(minimum[kr][r1]                [kc][c1],
                                           minimum[kr][r1]                [kc][c2 + 1 - (1 << kc)]);
                    double minimum_2 = min(minimum[kr][r2 + 1 - (1 << kr)][kc][c1],
                                           minimum[kr][r2 + 1 - (1 << kr)][kc][c2 + 1 - (1 << kc)]);
                    double minimum_height = min(minimum_1, minimum_2);

                    //std::cout << "ir " << ir << ", ic " << ic << ", min " << minimum_height << std::endl;

                    if(minimum_height < settings.maximum_height_map)
                    {
                        platform_tmp(ir, ic) = minimum_height;
                        opt_tmp += minimum_height * red_num;
                        if(red_num > 0 && minimum_height > 0) num_pin_tmp++;
                    }

                }
            }

            std::cout << dr << ",\t" << dc << ",\t" << opt_tmp << ",\t" << num_pin_tmp << std::endl;

            if(opt_tmp > opt_value || (std::abs(opt_tmp - opt_value) < settings.ZERO_EPS && num_pin_tmp < num_pin))
            {
                opt_value = opt_tmp;
                opt_r = dr * settings.pad_size / settings.xy_sample_num_each_pin;
                opt_c = dc * settings.pad_size / settings.xy_sample_num_each_pin;
                num_pin = num_pin_tmp;
                platform = platform_tmp;
            }
        }
    }

    //x is column, y is row
    dx = -opt_c;
    dz = -opt_r;

    //settings.toc();
    return;
}

#endif //EXAMPLE_MESH_LAYOUT_H
