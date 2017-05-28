//
// Created by 汪子琦 on 5/2/17.
//

#include "mesh_layout_base.h"

#ifndef SUPPORTER_MESH_LAYOUT_OPT_H
#define SUPPORTER_MESH_LAYOUT_OPT_H
#include <cmath>
#include <string>
using std::string;
using std::endl;
using Eigen::Vector2d;

class LayoutOptResult
{
public:

    LayoutOptResult()
    {
        dx = dz = material_save = 0;
        num_pin = 0;
        center = Vector2d(0, 0);
        angle = 0;
        edge_sample_num = 0;
    }

    LayoutOptResult(const LayoutOptResult &A)
    {
        dx = A.dx;
        dz = A.dz;
        material_save = A.material_save;
        angle = A.angle;
        edge_sample_num = A.edge_sample_num;
    }

public:

    bool operator < (LayoutOptResult &A)
    {
        if (material_save < A.material_save)
            return true;

        if (std::abs(material_save - A.material_save) < 1e-3)
        {
            if (num_pin > A.num_pin) return true;
            if (num_pin == A.num_pin && edge_sample_num > A.edge_sample_num) return true;
        }

//        if(num_pin > A.num_pin)
//            return true;
//
//        if(material_save < A.material_save)
//            return true;
        return false;
    }

	void writeInfo(const string &path_name)
	{
		std::ofstream fout;
		fout.open(path_name.c_str());
		fout << material_save << endl;
		fout << num_pin << endl;
		fout << edge_sample_num << endl;
		fout << dx << endl;
		fout << dz << endl;
		fout << angle << endl;
		fout << center.x() << " " << center.y() << endl;
		fout.close();
	}

	void readInfo(const string &path_name)
	{
		std::ifstream fin;
		fin.open(path_name.c_str());
		fin >> material_save;
		fin >> num_pin;
		fin >> edge_sample_num;
		fin >> dx;
		fin >> dz;
		fin >> angle;
		double cx, cy;
		fin >> cx >> cy;
		center = Vector2d(cx, cy);
		fin.close();
	}


public:

    double material_save;
    int num_pin;
    int edge_sample_num;

    double dx;
    double dz;

    double angle;
    Vector2d center;
};

class MeshLayoutOpt : public  MeshLayoutBase
{
public:

    MeshLayoutOpt();

public:
    void clear();

    void request_layout_xz_opt(LayoutOptResult &result);

    void request_layout_rotate_opt(LayoutOptResult &result);

public:

    void get_opt_hmap(MatrixXi &hmap);

    void get_opt_hmap(MatrixXd &hmap);

    void get_opt_smap(MatrixXi &smap);

    void get_opt_platform(MatrixXd &platform);

protected:

    LayoutOptResult opt_xy_layout(MatrixXi &hmap, MatrixXi &smap);

    void get_model_center(MatrixXi &hmap, Vector2d &center);

protected:
    MatrixXi opt_hmap;

    MatrixXi opt_smap;

    MatrixXd opt_platform;
};

MeshLayoutOpt::MeshLayoutOpt()
{

}

void MeshLayoutOpt::request_layout_xz_opt(LayoutOptResult &result)
{
    opt_hmap.setZero();
    opt_smap.setZero();
    opt_platform.setZero();

    if(height_map.isZero()) height_map_construction();
    if(support_map.isZero()) support_map_construction();

    result = opt_xy_layout(height_map, support_map);

    opt_hmap = height_map;
    opt_smap = support_map;


    transform_map_xz(opt_hmap, result.dx, result.dz, slicer.number_layer());
    transform_map_xz(opt_smap, result.dx, result.dz, 0);

    get_platform(opt_hmap, opt_smap, opt_platform);

    return;
}

void MeshLayoutOpt::request_layout_rotate_opt(LayoutOptResult &result)
{
    opt_hmap.setZero();
    opt_smap.setZero();
    opt_platform.setZero();

    if(height_map.isZero()) height_map_construction();
    if(support_map.isZero()) support_map_construction();

    Eigen::Vector2d center(0, 0);
    get_model_center(height_map, center);

    LayoutOptResult opt_result;
    for(int id = 0; id < settings.angle_sample_num; id++)
    {

        double angle = id * settings.angle_step;
        std::cout << "Rotate Opt\t" << angle << std::endl;

        Eigen::MatrixXi smap = support_map;
        Eigen::MatrixXi hmap = height_map;


        rotate_map_yaxis(hmap, angle, center, slicer.number_layer());
        rotate_map_yaxis(smap, angle, center, 0);

        LayoutOptResult tmp_result = opt_xy_layout(hmap, smap);
        if(opt_result < tmp_result)
        {
            opt_result = tmp_result;
            opt_result.center = center;
            opt_result.angle = angle;
        }
    }

    result = opt_result;

    opt_hmap = height_map;
    rotate_map_yaxis(opt_hmap, result.angle, result.center, slicer.number_layer());
    transform_map_xz(opt_hmap, result.dx, result.dz, slicer.number_layer());

    opt_smap = support_map;
    rotate_map_yaxis(opt_smap, result.angle, result.center, 0);
    transform_map_xz(opt_smap, result.dx, result.dz, 0);

    get_platform(opt_hmap, opt_smap, opt_platform);
}

void MeshLayoutOpt::get_opt_hmap(MatrixXi &hmap)
{
    assert(!opt_hmap.isZero());
    hmap = opt_hmap;
    return;
}

void MeshLayoutOpt::get_opt_hmap(MatrixXd &hmap) {

}

void MeshLayoutOpt::get_opt_smap(MatrixXi &smap) {
    assert(!opt_smap.isZero());
    smap = opt_smap;
    return;
}

void MeshLayoutOpt::get_opt_platform(MatrixXd &platform) {
    assert(!opt_platform.isZero());
    platform = opt_platform;
    return;
}

LayoutOptResult MeshLayoutOpt::opt_xy_layout(MatrixXi &hmap, MatrixXi &smap) {

    assert(!hmap.isZero() && !smap.isZero());

    //settings.tic("red Anchor");
    //Red Anchor
    Eigen::MatrixXi sum_smap;
    sum_smap = Eigen::MatrixXi::Zero(smap.rows(), smap.cols());

    for(int id = 0; id < smap.rows(); id++)
    {
        for(int jd = 0; jd < smap.cols(); jd++)
        {
            int L = jd > 0 ? sum_smap(id, jd - 1) : 0;
            int T = id > 0 ? sum_smap(id - 1, jd) : 0;
            int LT = id > 0 && jd > 0 ? sum_smap(id - 1, jd - 1) : 0;
            sum_smap(id, jd) = L + T - LT + smap(id, jd);
        }
    }
    //settings.toc();

    //settings.tic("MiniHmap");
    //min_hmap
    int n = settings.pillar_row * settings.xy_sample_num_each_pin;
    int m = settings.pillar_column * settings.xy_sample_num_each_pin;
    int kn = std::log(n) + 1;
    int km = std::log(m) + 1;

    int ****min_hmap;
    min_hmap = new int ***[kn];
    for(int jr = 0; jr < kn; jr++)
        min_hmap[jr] = new int **[n];

    for(int jr = 0; jr < kn; jr++)
    {
        for(int ir = 0; ir < n; ir++)
        {
            min_hmap[jr][ir] = new int *[km];
        }
    }

    for(int jr = 0; jr < kn; jr++)
    {
        for(int ir = 0; ir < n; ir++)
        {
            for(int jc = 0; jc < km; jc++)
                min_hmap[jr][ir][jc] = new int [m];
        }
    }

    for(int ir = 0; ir < n; ir++)
    {
        for (int ic = 0; ic < m; ic++)
        {
            min_hmap[0][ir][0][ic] = hmap(ir, ic);
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
                        min_hmap[jr][ir][jc][ic] = std::min(min_hmap[jr - 1][ir]                  [jc][ic],
                                                           min_hmap[jr - 1][ir + (1 << (jr - 1))][jc][ic]);
                    }
                    else
                    {
                        min_hmap[jr][ir][jc][ic] = std::min(min_hmap[jr][ir][jc - 1][ic],
                                                           min_hmap[jr][ir][jc - 1][ic + (1 << (jc - 1))]);
                    }
                }
            }
        }
    }
    //settings.toc();

    //optimization
    int Er = settings.edge_region_num;
    int Ec = settings.edge_region_num;

    //settings.tic("Main");
    LayoutOptResult opt_result;
    for(int dr = 0; dr < settings.xy_sample_num_each_pin; dr++)
    {
        for(int dc = 0; dc < settings.xy_sample_num_each_pin; dc++)
        {
            LayoutOptResult tmp_result;
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
                    int L =  r1 > 0 ?           sum_smap(r1 - 1, c2) : 0;
                    int T =  c1 > 0 ?           sum_smap(r2,     c1 - 1) : 0;
                    int LT = r1 > 0 && c1 > 0 ? sum_smap(r1 - 1, c1 - 1) : 0;
                    int red_num =   sum_smap(r2, c2) - L - T + LT;


                    int edge_red_num =   sum_smap(r2 - Er, c2 - Ec)
                                        -sum_smap(r1 + Er - 1, c2 - Ec)
                                        -sum_smap(r2 - Er, c1 + Ec - 1)
                                        +sum_smap(r1 + Er - 1, c1 + Ec - 1);
                    edge_red_num = red_num - edge_red_num;

                    //get min_hmap
                    int kr = std::log2(r2 - r1 + 1);
                    int kc = std::log2(c2 - c1 + 1);
                    int min_hmap_1 = std::min(min_hmap[kr][r1]                [kc][c1],
                                              min_hmap[kr][r1]                [kc][c2 + 1 - (1 << kc)]);
                    int min_hmap_2 = std::min(min_hmap[kr][r2 + 1 - (1 << kr)][kc][c1],
                                              min_hmap[kr][r2 + 1 - (1 << kr)][kc][c2 + 1 - (1 << kc)]);
                    int min_hmap_h = std::min(min_hmap_1, min_hmap_2);

                    //std::cout << "ir " << ir << ", ic " << ic << ", min " << min_hmap_height << std::endl;

                    if(min_hmap_h < slicer.number_layer() && red_num > 0)
                    {
                        double h = slicer.layer_pin_height(min_hmap_h);
                        //tmp_result.platform(ir, ic) = h;
                        tmp_result.material_save += h * red_num;
                        tmp_result.edge_sample_num+= edge_red_num;
                        if(red_num > 0 && min_hmap_h > 0) tmp_result.num_pin++;
                    }

                }
            }

            //std::cout << dr << ",\t" << dc << ",\t" << tmp.opt_value << ",\t" << tmp.num_pin << ",\t" << tmp.edge_sample_num << std::endl;

            if(opt_result < tmp_result)
            {
                opt_result = tmp_result;
                opt_result.dx = settings.int2mm(- dc * settings.sample_width);
                opt_result.dz = settings.int2mm(- dr * settings.sample_width);
            }
        }
    }
    //settings.toc();


    for(int jr = 0; jr < kn; jr++)
    {
        for(int ir = 0; ir < n; ir++)
        {
            for(int jc = 0; jc < km; jc++)
               delete[] min_hmap[jr][ir][jc];
        }
    }

    for(int jr = 0; jr < kn; jr++)
    {
        for(int ir = 0; ir < n; ir++)
            delete[] min_hmap[jr][ir];
    }

    for(int jr = 0; jr < kn; jr++)
        delete[] min_hmap[jr];

    delete[] min_hmap;

    return opt_result;
}

void MeshLayoutOpt::clear()
{
    MeshLayoutBase::clear();
    opt_hmap.setZero();
    opt_smap.setZero();
    opt_platform.setZero();
}

void MeshLayoutOpt::get_model_center(MatrixXi &hmap, Vector2d &center)
{
    int num_pts = 0;
    center = Vector2d(0, 0);
    for(int ir = 0; ir < hmap.rows(); ir++)
    {
        for(int ic = 0; ic < hmap.cols(); ic++)
        {
            if(hmap(ir, ic) < slicer.number_layer())
            {
                center += Vector2d(settings.pin_center_x(ic), settings.pin_center_y(ir));
                num_pts++;
            }
        }
    }
    center /= num_pts;
}

#endif //SUPPORTER_MESH_LAYOUT_OPT_H
