//
// Created by 汪子琦 on 4/20/17.
//

#ifndef EXAMPLE_MESH_LAYOUT_H
#define EXAMPLE_MESH_LAYOUT_H

#define DRAW_SUPPORT 1

#include "mesh_slicer.h"
#include <cmath>
#include <queue>
#include <Eigen/Dense>

#if DRAW_SUPPORT
#include "scene_organizer.h"
#endif


typedef struct tagHeightMapNode
{
    int R;
    int C;
    int XYsize;
    int layer;
    bool new_layer;
    ClipperLib::Paths polys;
}HeightMapNode;

typedef struct tagLayoutOptOutput
{
public:

    tagLayoutOptOutput()
    {
        dx = dy = opt_value = 0;
        num_pin = 0;
        center = Eigen::Vector2d(0, 0);
        angle = 0;
        edge_sample_num = 0;
        rotate = false;
    }

public:

    bool operator < (tagLayoutOptOutput &A)
    {
        if (opt_value < A.opt_value)
            return true;

        if (std::abs(opt_value - A.opt_value) < settings.ZERO_EPS)
        {
            if (num_pin > A.num_pin) return true;
            if (num_pin == A.num_pin && edge_sample_num > A.edge_sample_num) return true;
        }

        return false;
    }

public:

    double opt_value;
    int num_pin;
    int edge_sample_num;
    Eigen::MatrixXd platform;

    double dx;
    double dy;

    bool rotate;
    double angle;
    Eigen::Vector2d center;

    Settings settings;
}LayoutOptOutput;

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

    void get_height_map(Eigen::MatrixXd &height_map);

    void get_red_map(Eigen::MatrixXi &red_map);

    void get_platform(Eigen::MatrixXd &platform, Eigen::MatrixXd &height_map, Eigen::MatrixXi &red_map);

    LayoutOptOutput xy_layout();

    LayoutOptOutput rotate_layout();

#if DRAW_SUPPORT
    LayoutOptOutput xy_opt_support(Eigen::MatrixXd &V, Eigen::MatrixXi &F);

    LayoutOptOutput non_pin_support(Eigen::MatrixXd &V, Eigen::MatrixXi &F);

    LayoutOptOutput stand_support(Eigen::MatrixXd &V, Eigen::MatrixXi &F);

    void draw_support(Eigen::MatrixXd &platform,
                      Eigen::MatrixXd &height_map,
                      Eigen::MatrixXi &red_map,
                      Eigen::MatrixXd &V,
                      Eigen::MatrixXi &F);
#endif

protected:

    void move_height_map(Eigen::MatrixXd &map, double x, double z);

    void move_red_map(Eigen::MatrixXi &map, double x, double z);

    void rotate_height_map(Eigen::MatrixXd &map, Eigen::Vector2d &center, double angle = 0);

    void rotate_red_map(Eigen::MatrixXi &map, Eigen::Vector2d &center, double angle = 0);

    void height_map_construction();

    void red_map_construction();

    LayoutOptOutput xy_layout(Eigen::MatrixXd &height_map, Eigen::MatrixXi &red_map);

private:

    //void layer_projecting(Eigen::MatrixXi &map, std::vector<ClipperLib::Paths> &slices);

public:

    MeshSlicer slicer_;

    Eigen::MatrixXi height_map_;

    Eigen::MatrixXi red_map_;

    Settings settings;
};

void MeshLayout::get_height_map(Eigen::MatrixXd &height_map)
{
    if(height_map_.isZero())
    {
        height_map_construction();
    }
    height_map.resize(height_map_.rows(), height_map_.cols());
    for(int ir = 0; ir < height_map_.rows(); ir++)
        for(int ic = 0; ic < height_map_.cols(); ic++)
            height_map(ir, ic) = slicer_.layer_pin_height(height_map_(ir, ic));
    return;
}

void MeshLayout::get_red_map(Eigen::MatrixXi &red_map)
{

    if(red_map_.isZero())
    {
        red_map_construction();
    }
    red_map = red_map_;
    return;
}

LayoutOptOutput MeshLayout::xy_layout(Eigen::MatrixXd &height_map, Eigen::MatrixXi &red_map)
{
    settings.tic("xy layout");
    assert(!height_map.isZero() && !red_map.isZero());

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
                        minimum[jr][ir][jc][ic] = std::min(minimum[jr - 1][ir]                  [jc][ic],
                                                      minimum[jr - 1][ir + (1 << (jr - 1))][jc][ic]);
                    }
                    else
                    {
                        minimum[jr][ir][jc][ic] = std::min(minimum[jr][ir][jc - 1][ic],
                                                      minimum[jr][ir][jc - 1][ic + (1 << (jc - 1))]);
                    }
                }
            }
        }
    }


    //optimization
    int Er = settings.xy_sample_num_each_pin / 10;
    int Ec = settings.xy_sample_num_each_pin / 10;

    LayoutOptOutput opt;
    for(int dr = 0; dr < settings.xy_sample_num_each_pin; dr++)
    {
        for(int dc = 0; dc < settings.xy_sample_num_each_pin; dc++)
        {
            LayoutOptOutput tmp;
            tmp.platform = Eigen::MatrixXd::Zero(settings.pillar_row, settings.pillar_column);
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


                    int edge_red_num =  red_anchor(r2 - Er, c2 - Ec)
                                       -red_anchor(r1 + Er - 1, c2 - Ec)
                                       -red_anchor(r2 - Er, c1 + Ec - 1)
                                       +red_anchor(r1 + Er - 1, c1 + Ec - 1);
                    edge_red_num = red_num - edge_red_num;

                    //get minimum
                    int kr = std::log2(r2 - r1 + 1);
                    int kc = std::log2(c2 - c1 + 1);
                    double minimum_1 = std::min(minimum[kr][r1]                [kc][c1],
                                           minimum[kr][r1]                [kc][c2 + 1 - (1 << kc)]);
                    double minimum_2 = std::min(minimum[kr][r2 + 1 - (1 << kr)][kc][c1],
                                           minimum[kr][r2 + 1 - (1 << kr)][kc][c2 + 1 - (1 << kc)]);
                    double minimum_height = std::min(minimum_1, minimum_2);

                    //std::cout << "ir " << ir << ", ic " << ic << ", min " << minimum_height << std::endl;

                    if(minimum_height < settings.maximum_height_map && red_num > 0)
                    {
                        tmp.platform(ir, ic) = minimum_height;
                        tmp.opt_value += minimum_height * red_num;
                        tmp.edge_sample_num+= edge_red_num;
                        if(red_num > 0 && minimum_height > 0) tmp.num_pin++;
                    }

                }
            }

            std::cout << dr << ",\t" << dc << ",\t" << tmp.opt_value << ",\t" << tmp.num_pin << ",\t" << tmp.edge_sample_num << std::endl;

            if(opt < tmp)
            {
                opt = tmp;
                opt.dx = settings.int2mm(- dc * settings.sample_width);
                opt.dy = settings.int2mm(- dr * settings.sample_width);
            }
        }
    }

    settings.toc();
    return opt;
}

void MeshLayout::height_map_construction()
{
    assert(!slicer_.empty());

    int nc = settings.pillar_column  * settings.xy_sample_num_each_pin;
    int nr = settings.pillar_row  * settings.xy_sample_num_each_pin;
    height_map_ =  Eigen::MatrixXi::Ones(nr, nc) * slicer_.number_layer();

    VecMatrixXi bitmaps;
    slicer_.get_bitmaps(bitmaps);

    for(int layer = 0; layer < slicer_.number_layer(); layer++)
    {
        for(int ir = 0; ir < nr; ir++)
        {
            for(int ic = 0; ic < nc; ic++)
            {
                if(bitmaps[layer](ir, ic) && height_map_(ir, ic) > layer)
                {
                    height_map_(ir, ic) = layer;
                }
            }
        }
    }

    return;
}

void MeshLayout::red_map_construction()
{
    int nc = settings.pillar_column  * settings.xy_sample_num_each_pin;
    int nr = settings.pillar_row  * settings.xy_sample_num_each_pin;
    red_map_ =  Eigen::MatrixXi::Zero(nr, nc);

    std::vector<ClipperLib::Paths> bottom_half;
    slicer_.get_bottom_half(bottom_half);

    ClipperLib::Paths downside;
    for(int layer = 1; layer < slicer_.number_layer(); layer++)
    {
        ClipperLib::Clipper clipper;
        clipper.AddPaths(downside, ClipperLib::ptSubject, true);
        clipper.AddPaths(bottom_half[layer], ClipperLib::ptClip, true);
        clipper.Execute(ClipperLib::ctUnion, downside, ClipperLib::pftPositive, ClipperLib::pftPositive);
    }

    Eigen::MatrixXi bitmap;
    ScanLineFill fill(false);
    fill.polygon_fill(downside, red_map_);
}

void MeshLayout::rotate_height_map(Eigen::MatrixXd &map, Eigen::Vector2d &center, double angle)
{
    Eigen::MatrixXd new_map;
    int nc = settings.pillar_column  * settings.xy_sample_num_each_pin;
    int nr = settings.pillar_row  * settings.xy_sample_num_each_pin;

    //init
    new_map =  Eigen::MatrixXd::Ones(nr, nc) * settings.maximum_height_map;

    //center
    Eigen::Matrix2d rot_mat;
    rot_mat <<  cos(angle), -sin(angle),
                sin(angle),  cos(angle);
    for(int ir = 0; ir < nr; ir++)
    {
        for (int ic = 0; ic < nc; ic++)
        {
            Eigen::Vector2d p(settings.pin_center_x(ic), settings.pin_center_y(ir));
            p -= center;
            p = rot_mat * p;
            p += center;

            int new_ir = std::floor(p(1) / settings.sample_width);
            int new_ic = std::floor(p(0) / settings.sample_width);
            if(0 <= new_ir && new_ir < nr && 0 <= new_ic && new_ic < nc)
                new_map(new_ir, new_ic) = map(ir, ic);
        }
    }

    map = new_map;
    return;
}

void MeshLayout::rotate_red_map(Eigen::MatrixXi &map, Eigen::Vector2d &center, double angle)
{
    Eigen::MatrixXi new_map;
    int nc = settings.pillar_column  * settings.xy_sample_num_each_pin;
    int nr = settings.pillar_row  * settings.xy_sample_num_each_pin;

    //init
    new_map =  Eigen::MatrixXi::Zero(nr, nc);

    //center
    Eigen::Matrix2d rot_mat;
    rot_mat <<  cos(angle), -sin(angle),
                sin(angle),  cos(angle);
    for(int ir = 0; ir < nr; ir++)
    {
        for (int ic = 0; ic < nc; ic++)
        {
            Eigen::Vector2d p(settings.pin_center_x(ic), settings.pin_center_y(ir));
            p -= center;
            p = rot_mat * p;
            p += center;

            int new_ir = std::floor(p(1) / settings.sample_width);
            int new_ic = std::floor(p(0) / settings.sample_width);
            if(0 <= new_ir && new_ir < nr && 0 <= new_ic && new_ic < nc)
                new_map(new_ir, new_ic) = map(ir, ic);
        }
    }

    map = new_map;
    return;
}

LayoutOptOutput MeshLayout::xy_layout()
{
    if(height_map_.isZero()) height_map_construction();
    if(red_map_.isZero()) red_map_construction();
    Eigen::MatrixXd height_map;
    get_height_map(height_map);
    return xy_layout(height_map, red_map_);
}

LayoutOptOutput MeshLayout::rotate_layout() {
    if(height_map_.isZero()) height_map_construction();
    if(red_map_.isZero()) red_map_construction();

    Eigen::Vector2d center(0, 0);
    int num_pts = 0;
    for(int ir = 0; ir < height_map_.rows(); ir++)
    {
        for(int ic = 0; ic < height_map_.cols(); ic++)
        {
            if(height_map_(ir, ic) < settings.maximum_height_map)
            {
                center += Eigen::Vector2d(settings.pin_center_x(ic), settings.pin_center_y(ir));
                num_pts++;
            }
        }
    }
    center /= num_pts;

    LayoutOptOutput opt;
    for(int id = 0; id < settings.angle_sample_num; id++)
    {
        std::cout << "angle: " << id * 360 / settings.angle_sample_num << std::endl;
        double angle = id * settings.angle_step;
        Eigen::MatrixXi red_map = red_map_;
        Eigen::MatrixXd height_map; get_height_map(height_map);
        LayoutOptOutput tmp;
        rotate_height_map(height_map, center, angle);
        rotate_red_map(red_map, center, angle);
        tmp = xy_layout(height_map, red_map);
        if(opt < tmp)
        {
            opt = tmp;
            opt.center = Eigen::Vector2d(settings.int2mm(center(0)), settings.int2mm(center(1)));
            opt.angle = angle;
        }
    }

    return opt;
}

#if DRAW_SUPPORT
LayoutOptOutput MeshLayout::xy_opt_support(Eigen::MatrixXd &V, Eigen::MatrixXi &F)
{
    LayoutOptOutput opt = xy_layout();
    Eigen::MatrixXd height_map;
    Eigen::MatrixXi red_map; get_red_map(red_map);

    if(height_map_.isZero()) height_map_construction();
    height_map.resize(height_map_.rows(), height_map_.cols());
    for(int ir = 0; ir < height_map_.rows(); ir++)
        for(int ic = 0; ic < height_map_.cols(); ic++)
            height_map(ir, ic) = height_map_(ir, ic) * settings.layer_height;

    move_height_map(height_map, opt.dx, opt.dy);
    move_red_map(red_map, opt.dx, opt.dy);

    Eigen::MatrixXd platform = opt.platform;
    draw_support(platform, height_map, red_map, V, F);
    return opt;
}

LayoutOptOutput MeshLayout::stand_support(Eigen::MatrixXd &V, Eigen::MatrixXi &F)
{
    LayoutOptOutput opt;

    Eigen::MatrixXd height_map;
    Eigen::MatrixXi red_map; get_red_map(red_map);

    if(height_map_.isZero()) height_map_construction();
    height_map.resize(height_map_.rows(), height_map_.cols());
    for(int ir = 0; ir < height_map_.rows(); ir++)
        for(int ic = 0; ic < height_map_.cols(); ic++)
            height_map(ir, ic) = height_map_(ir, ic) * settings.layer_height;

    get_platform(opt.platform, height_map, red_map);
    Eigen::MatrixXd platform = opt.platform;
    draw_support(platform, height_map, red_map, V, F);
    return opt;
}

LayoutOptOutput MeshLayout::non_pin_support(Eigen::MatrixXd &V, Eigen::MatrixXi &F)
{
    LayoutOptOutput opt;
    opt.platform = Eigen::MatrixXd::Zero(settings.pillar_row, settings.pillar_column);

    Eigen::MatrixXd height_map;
    Eigen::MatrixXi red_map; get_red_map(red_map);

    if(height_map_.isZero()) height_map_construction();
    height_map.resize(height_map_.rows(), height_map_.cols());
    for(int ir = 0; ir < height_map_.rows(); ir++)
        for(int ic = 0; ic < height_map_.cols(); ic++)
            height_map(ir, ic) = height_map_(ir, ic) * settings.layer_height;

    Eigen::MatrixXd platform = opt.platform;
    draw_support(platform, height_map, red_map, V, F);
    return opt;
}

void MeshLayout::draw_support(Eigen::MatrixXd &platform,
                              Eigen::MatrixXd &height_map,
                              Eigen::MatrixXi &red_map,
                              Eigen::MatrixXd &V,
                              Eigen::MatrixXi &F) {
    SceneOrganizer organizer;
    VecMatrixXi Fs;
    VecMatrixXd Vs;
    for(int ir = 0; ir < red_map.rows(); ir++)
    {
        for(int ic = 0; ic < red_map.cols(); ic++)
        {
            int pin_r = ir / settings.xy_sample_num_each_pin;
            int pin_c = ic / settings.xy_sample_num_each_pin;
            if(red_map(ir, ic) > 0)
            {
                Eigen::MatrixXd tV;
                Eigen::MatrixXi tF;
                Eigen::Vector3d bottom(settings.int2mm(settings.pin_center_x(ic)),
                                       platform(pin_r, pin_c),
                                       settings.int2mm(settings.pin_center_y(ir)));
                double height = height_map(ir, ic) - platform(pin_r, pin_c) - settings.layer_height;

                if(height >= settings.layer_height) {
                    organizer.draw_pillar(bottom, settings.pad_size / settings.xy_sample_num_each_pin, height, tV, tF);
                    Vs.push_back(tV);
                    Fs.push_back(tF);
                }
            }
        }
    }
    V.setZero();
    F.setZero();
    int nv = 0, nf = 0;
    for(int id = 0; id < Vs.size(); id ++)
    {
        nv += Vs[id].rows(); nf += Fs[id].rows();
    }
    V.resize(nv, 3); F.resize(nf, 3);
    nv = 0; nf = 0;
    for(int id = 0; id < Vs.size(); id ++)
    {
        //F
        for(int jd = 0; jd < Fs[id].rows(); jd++)
        {
            F.row(nf ++) = Fs[id].row(jd) + Eigen::RowVector3i(nv, nv, nv);
        }

        //V
        for(int jd = 0; jd < Vs[id].rows(); jd++)
        {
            V.row(nv ++) = Vs[id].row(jd);
        }
    }

    return;
}

#endif

void MeshLayout::move_height_map(Eigen::MatrixXd &map, double x, double z)
{
    Eigen::MatrixXd new_map;
    int nc = settings.pillar_column  * settings.xy_sample_num_each_pin;
    int nr = settings.pillar_row  * settings.xy_sample_num_each_pin;

    //init
    new_map =  Eigen::MatrixXd::Ones(nr, nc) * settings.maximum_height_map;

    for(int ir = 0; ir < nr; ir++)
    {
        for (int ic = 0; ic < nc; ic++)
        {
            Eigen::Vector2d p(settings.pin_center_x(ic), settings.pin_center_y(ir));
            p+= Eigen::Vector2d(settings.mm2int(x), settings.mm2int(z));

            int new_ir = std::floor(p(1) / settings.sample_width);
            int new_ic = std::floor(p(0) / settings.sample_width);
            if(0 <= new_ir && new_ir < nr && 0 <= new_ic && new_ic < nc)
                new_map(new_ir, new_ic) = map(ir, ic);
        }
    }

    map = new_map;
    return;
}

void MeshLayout::move_red_map(Eigen::MatrixXi &map, double x, double z)
{
    Eigen::MatrixXi new_map;
    int nc = settings.pillar_column  * settings.xy_sample_num_each_pin;
    int nr = settings.pillar_row  * settings.xy_sample_num_each_pin;

    //init
    new_map =  Eigen::MatrixXi::Zero(nr, nc);

    for(int ir = 0; ir < nr; ir++)
    {
        for (int ic = 0; ic < nc; ic++)
        {
            Eigen::Vector2d p(settings.pin_center_x(ic), settings.pin_center_y(ir));
            p+= Eigen::Vector2d(settings.mm2int(x), settings.mm2int(z));

            int new_ir = std::floor(p(1) / settings.sample_width);
            int new_ic = std::floor(p(0) / settings.sample_width);
            if(0 <= new_ir && new_ir < nr && 0 <= new_ic && new_ic < nc)
                new_map(new_ir, new_ic) = map(ir, ic);
        }
    }

    map = new_map;
    return;
}

void MeshLayout::get_platform(Eigen::MatrixXd &platform, Eigen::MatrixXd &height_map, Eigen::MatrixXi &red_map)
{
    platform = Eigen::MatrixXd::Zero(9, 11);

    for(int ir = 0; ir < settings.pillar_row; ir++)
    {
        for(int ic = 0; ic < settings.pillar_column; ic++)
        {
            int r1 = ir * settings.xy_sample_num_each_pin  ;
            int r2 = (ir + 1) * settings.xy_sample_num_each_pin - 1;
            int c1 = ic * settings.xy_sample_num_each_pin ;
            int c2 = (ic + 1) * settings.xy_sample_num_each_pin - 1;

            int num_red = 0;
            double minimum_height = settings.MAX_DOUBLE;
            for(int id = r1; id <= r2; id++)
            {
                for(int jd = c1; jd <= c2; jd++)
                {
                    if(minimum_height > height_map(id, jd)) minimum_height = height_map(id, jd);
                    num_red += red_map(id, jd);
                }
            }

            if(num_red > 0)
                platform(ir, ic) = (int)(minimum_height / settings.pillar_standard_height) * settings.pillar_standard_height;
        }
    }

    return;
}



//void MeshLayout::layer_projecting(Eigen::MatrixXi &map, std::vector<ClipperLib::Paths> &slices)
//{
//    settings.tic("\nRasterize Time: ");
//    std::queue<HeightMapNode> Queue;
//    map = Eigen::MatrixXi::Zero(settings.pillar_row  * settings.xy_sample_num_each_pin,
//                                settings.pillar_column * settings.xy_sample_num_each_pin);
//
//    for(int row_id = 0; row_id < settings.pillar_row; row_id++)
//    {
//        for(int col_id = 0; col_id < settings.pillar_column; col_id++)
//        {
//            HeightMapNode node;
//            node.R = row_id * settings.xy_sample_num_each_pin;
//            node.C = col_id * settings.xy_sample_num_each_pin;
//            node.XYsize = settings.xy_sample_num_each_pin;
//            node.layer = 0;
//            node.new_layer = true;
//            node.polys.clear();
//            Queue.push(node);
//        }
//    }
//
//    while(!Queue.empty()) {
//        HeightMapNode u = Queue.front();
//        Queue.pop();
//
//        if (u.layer >= slices.size()) {
//            map.block(u.R, u.C, u.XYsize, u.XYsize) = Eigen::MatrixXi::Ones(u.XYsize, u.XYsize) * -1;
//            continue;
//        }
//
//        ClipperLib::Clipper clipper;
//        if (u.new_layer) {
//            clipper.AddPaths(slices[u.layer], ClipperLib::ptClip, true);
//        } else {
//            clipper.AddPaths(u.polys, ClipperLib::ptClip, true);
//        }
//
//        ClipperLib::Path square;
//        int L =  settings.mm2int(settings.pad_size / settings.xy_sample_num_each_pin * u.C);
//        int R =  settings.mm2int(settings.pad_size / settings.xy_sample_num_each_pin * (u.C + u.XYsize));
//        int T =  settings.mm2int(settings.pad_size / settings.xy_sample_num_each_pin * u.R);
//        int B =  settings.mm2int(settings.pad_size / settings.xy_sample_num_each_pin * (u.R + u.XYsize));
//        square.push_back(ClipperLib::IntPoint(L, T));
//        square.push_back(ClipperLib::IntPoint(R, T));
//        square.push_back(ClipperLib::IntPoint(R, B));
//        square.push_back(ClipperLib::IntPoint(L, B));
//
//        clipper.AddPath(square, ClipperLib::ptSubject, true);
//
//        //computing the intersection
//        ClipperLib::Paths intersection;
//        clipper.Execute(ClipperLib::ctIntersection, intersection, ClipperLib::pftPositive, ClipperLib::pftPositive);
//        ClipperLib::SimplifyPolygons(intersection);
//
//        if (intersection.empty())
//        {
//            u.layer++;
//            u.new_layer = true;
//            Queue.push(u);
//        }
//        else
//        {
//            //check the collision
//            double area_insec = 0;
//            for (int id = 0; id < intersection.size(); id++)
//                area_insec += ClipperLib::Area(intersection[id]);
//            area_insec = std::abs(area_insec);
//
//            double area_squre = std::abs(ClipperLib::Area(square));
//            if (std::abs(area_insec - area_squre) < settings.ZERO_EPS || (u.XYsize == 1))
//            {
//                map.block(u.R, u.C, u.XYsize, u.XYsize) = Eigen::MatrixXi::Ones(u.XYsize, u.XYsize) * u.layer;
//            }
//            else
//            {
//                for (int ir = 0; ir < 2; ir++) {
//                    for (int ic = 0; ic < 2; ic++) {
//                        HeightMapNode v = u;
//                        v.R += ir * u.XYsize / 2;
//                        v.C += ic * u.XYsize / 2;
//                        v.XYsize /= 2;
//                        v.polys = intersection;
//                        v.new_layer = false;
//                        Queue.push(v);
//                    }
//                }
//            }
//        }
//    }
//    settings.toc();
//}


#endif //EXAMPLE_MESH_LAYOUT_H
