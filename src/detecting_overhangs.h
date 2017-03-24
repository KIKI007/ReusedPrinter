//
// Created by 汪子琦 on 2/28/17.
//

#ifndef detect_overhangs_include_h
#define detect_overhangs_include_h

#include <igl/per_face_normals.h>
#include <Eigen/Core>

#include <cmath>
#include <iostream>
#include <random>

#include "settings.h"

class DetectingOverhangs
{
public:

    DetectingOverhangs()
    {
        nd = std::uniform_real_distribution<double>(0, 1);
    }

    ~DetectingOverhangs()
    {

    }
public:

    //C return the color of each face
    //red: overhanging faces
    //blue: non-overhanging faces
    void colorOverhangingFaces(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &C);

    //Out return the overhanging sampled points
    void sample(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &Out);

    void sample(Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::vector< Eigen::Vector3d, Eigen::aligned_allocator< Eigen::Vector3d>> &Out);
private:

    //Fi return the type of each face
    //function return the number of overhanging faces
    //1: overhanging faces
    //0: non-overhanging faces
    int findOverhangingFaces(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &Fi);

    //OF return each overhaning faces' vertices' indexs;
    //function return the number of overhanging faces
    int getOverhangingFaces(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXi &OF);

    double get_normal_distribution()
    {
        return nd(generator);
    }

private:
    //normal distribution
    std::uniform_real_distribution<double> nd;

    std::default_random_engine generator;

    Settings settings_;
};

int DetectingOverhangs::findOverhangingFaces(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &Fi)
{
    Eigen::MatrixXd faces_normals;
    Fi = Eigen::VectorXi::Zero(F.rows());
    igl::per_face_normals(V, F, faces_normals);
    int overhangs_num = 0;


    for(size_t id = 0; id < faces_normals.rows(); id++)
    {
        if(faces_normals(id, 1) * (-1) > std::cos(settings_.face_overhang_angle))
        {
            //this face is a overhanging faces(1)
            Fi(id) = 1;
            overhangs_num++;
        }
        else
        {
            //this face isn't a overhanging faces(0)
            Fi(id) = 0;
        }
    }
    return overhangs_num;
}

int DetectingOverhangs::getOverhangingFaces(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXi &OF)
{
    Eigen::VectorXi Fi;
    int n = findOverhangingFaces(V, F, Fi);
    OF = Eigen::MatrixXi::Zero(n , 3);
    size_t id = 0;
    for(size_t fd = 0; fd < F.rows(); fd++)
    {
        if(Fi(fd))
            OF.row(id++) = F.row(fd);
    }
    return n;
}

void DetectingOverhangs::colorOverhangingFaces(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &C)
{
    Eigen::VectorXi Fi;
    findOverhangingFaces(V, F, Fi);
    C = Eigen::MatrixXd::Zero(F.rows(), 3);

    for(size_t id = 0; id < Fi.rows(); id++)
    {
        C.row(id) = Fi(id) ? Eigen::RowVector3d(1, 0, 0) : Eigen::RowVector3d(0 , 0, 1);
    }
    return;
}


void DetectingOverhangs::sample(Eigen::MatrixXd &V,
                                 Eigen::MatrixXi &F,
                                 std::vector< Eigen::Vector3d, Eigen::aligned_allocator< Eigen::Vector3d>> &Out)
{
    Eigen::MatrixXd S;
    sample(V, F, S);

    Out.resize(S.rows());
    for(size_t kd = 0; kd < S.rows(); kd++)
    {
        Out[kd] = S.row(kd).transpose();
    }
    return;
}

void DetectingOverhangs::sample(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &Out)
{
    Eigen::MatrixXi OF; //whether faces is a overhanging faces
    int n = getOverhangingFaces(V, F, OF);
    Eigen::VectorXd AF(OF.rows()); //area of each overhanging faces

    double total_area = 0;
    //calculate the are of each overhanging faces;
    for(size_t fd = 0; fd < OF.rows(); fd++) {
        Eigen::RowVector3d p01 = V.row(OF(fd, 1)) - V.row(OF(fd, 0));
        Eigen::RowVector3d p02 = V.row(OF(fd, 2)) - V.row(OF(fd, 0));
        AF(fd) = p01.cross(p02).norm();
        total_area += AF(fd);
    }

    //using the faces' area to determine the number of sampling points
    int n_sample = 0; //number of sample points
    for(size_t fd = 0; fd < AF.rows(); fd++)
    {
        int times = (int) (AF(fd) / total_area * settings_.expected_sample_num);
        //times = (times == 0) ? 1 : times;
        //or using math function to map between AF to times
        n_sample += times;
    }

    std::cout << "n_sample: \t" << n_sample << std::endl;

    //sample
    Out = Eigen::MatrixXd::Zero(n_sample, 3);
    int id = 0;
    Eigen::RowVector3d A, B, C;
    for(size_t fd = 0; fd < OF.rows(); fd++) {
        int times = (int) (AF(fd) / total_area * settings_.expected_sample_num);
        //times = (times == 0) ? 1 : times;
        A = V.row(OF(fd, 0));
        B = V.row(OF(fd, 1));
        C = V.row(OF(fd, 2));
        while (times--)
        {
            double u = get_normal_distribution();
            double v = get_normal_distribution();
            Out.row(id++) = A * (1 - sqrt(u)) + B * (sqrt(u) * (1 - v))  + C * v * sqrt(u);
        }
    }

    return;
}

#endif //EXAMPLE_DETECT_OVERHANGS_H
