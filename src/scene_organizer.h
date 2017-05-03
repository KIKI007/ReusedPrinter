//
// Created by 汪子琦 on 3/2/17.
//

#ifndef SUPPORTER_SCENE_ORGANIZER_H
#define SUPPORTER_SCENE_ORGANIZER_H

#include "Eigen/Core"

class SceneOrganizer {

public:

    SceneOrganizer()
    {
        V_ = Eigen::MatrixXd(0, 0);
        F_ = Eigen::MatrixXi(0, 0);
        C_ = Eigen::MatrixXd(0, 0);
    }


public:

    void add_mesh(const Eigen::MatrixXd &V, Eigen::MatrixXi F, Eigen::RowVector3d c)
    {
        size_t n = Vrows();
        Eigen::MatrixXd C(F.rows(), 3);
        for(size_t id = 0; id < F.rows(); id++)
        {
            F.row(id) = F.row(id) + Eigen::RowVector3i(n, n, n);
            C.row(id) = c;
        }

        Eigen::MatrixXd tmpV_(V_.rows() + V.rows(), 3);
        tmpV_.block(0, 0, V_.rows(), V_.cols()) = V_;
        tmpV_.block(V_.rows(), 0, V.rows(), V.cols()) = V;
        V_ = tmpV_;

        Eigen::MatrixXi tmpF_(F_.rows() + F.rows(), 3);
        tmpF_.block(0, 0, F_.rows(), F_.cols()) = F_;
        tmpF_.block(F_.rows(), 0, F.rows(), F.cols()) = F;
        F_ = tmpF_;


        Eigen::MatrixXd tmpC_(C_.rows() + C.rows(), 3);
        tmpC_.block(0, 0, C_.rows(), C_.cols()) = C_;
        tmpC_.block(C_.rows(), 0, C.rows(), C.cols()) = C;
        C_ = tmpC_;

        return;
    }

    void get_mesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &C)
    {
        V = V_;
        F = F_;
        C = C_;
        return;
    }

public:


    void draw_pillar(Eigen::Vector3d &bottom_center, double bottom_width, double height, Eigen::MatrixXd &V, Eigen::MatrixXi &F)
    {
        Eigen::Vector3d pts[8];
        pts[0] = bottom_center + Eigen::Vector3d(-bottom_width/ 2, 0,  bottom_width / 2);
        pts[1] = bottom_center + Eigen::Vector3d( bottom_width/ 2, 0,  bottom_width / 2);
        pts[2] = bottom_center + Eigen::Vector3d( bottom_width/ 2, 0, -bottom_width / 2);
        pts[3] = bottom_center + Eigen::Vector3d(-bottom_width/ 2, 0, -bottom_width / 2);

        pts[4] = bottom_center + Eigen::Vector3d(-bottom_width/ 2, height,  bottom_width / 2);
        pts[5] = bottom_center + Eigen::Vector3d( bottom_width/ 2, height,  bottom_width / 2);
        pts[6] = bottom_center + Eigen::Vector3d( bottom_width/ 2, height, -bottom_width / 2);
        pts[7] = bottom_center + Eigen::Vector3d(-bottom_width/ 2, height, -bottom_width / 2);

        F.resize(12, 3);
        //lateral
        for(int id = 0; id < 4; id++)
        {
            F(2 * id, 0) = id;
            F(2 * id, 1) = (id + 1) % 4;
            F(2 * id, 2) = 4 + (id + 1) % 4;

            F(2 * id + 1, 0) = 4 + (id + 1) % 4;
            F(2 * id + 1, 1) = 4 +  id % 4;
            F(2 * id + 1, 2) = id;
        }

        //bottom
        F(8,  0) = 0; F(8,  1) = 3; F(8,  2) = 2;
        F(9,  0) = 2; F(9,  1) = 1; F(9,  2) = 0;
        F(10, 0) = 4; F(10, 1) = 5; F(10, 2) = 6;
        F(11, 0) = 6; F(11, 1) = 7; F(11, 2) = 4;

        V.resize(8, 3);
        for(int id = 0; id < 8; id++)
            V.row(id) = pts[id];

        return;
    }

public:

    size_t Vrows() {return V_.rows();}

private:

    Eigen::MatrixXd V_;

    Eigen::MatrixXi F_;

    Eigen::MatrixXd C_;

private:


};

class GeneratingPlatform
{
public:

    GeneratingPlatform()
    {

    }

private:

    int num_vertices_single_obj()
    {
        return 2 * settings_.circle_vertices_num + 8;
    }

    int num_faces_single_obj()
    {
        size_t num_face_cube = 2 * 2 + 4 * 2;
        size_t num_face_pillar = (settings_.circle_vertices_num - 2) + 2 * settings_.circle_vertices_num;
        return num_face_cube + num_face_pillar;
    }

public:

    void draw_platform(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &H)
    {
        int width = H.cols();
        int height = H.rows();

        V = Eigen::MatrixXd::Zero(H.size() * num_vertices_single_obj(), 3);
        F = Eigen::MatrixXi::Zero(H.size() * num_faces_single_obj(), 3);

        double px = 0;
        double pz = 0;
        int o_id = 0;
        for(int iw = 0; iw < width; iw++)
        {
            for(int ih = 0; ih < height; ih++)
            {
                px = (iw + 0.5) * settings_.pad_size;
                pz = (ih + 0.5) * settings_.pad_size;

                Eigen::MatrixXd sV;
                Eigen::MatrixXi sF;
                draw_pad_pillar(sV, sF);

                sF += Eigen::MatrixXi::Ones(sF.rows(), sF.cols()) * o_id * num_vertices_single_obj();
                for(size_t id = 0; id < sV.rows(); id++)
                {
                    sV.row(id) += Eigen::RowVector3d(px, H(ih, iw) , pz);
                }

                V.block(o_id * num_vertices_single_obj(), 0, sV.rows(), sV.cols()) = sV;
                F.block(o_id * num_faces_single_obj(), 0, sF.rows(), sF.cols()) = sF;
                o_id++;
            }
        }

        return;
    }

    void draw_pad_pillar(Eigen::MatrixXd &V, Eigen::MatrixXi &F) {

        V = Eigen::MatrixXd::Zero(num_vertices_single_obj(), 3);

        F = Eigen::MatrixXi::Zero(num_faces_single_obj(), 3);

        //cube
        double cube_size = settings_.pad_size - 0.2;
        int dx[4] = {-1, -1, 1, 1};
        int dz[4] = {-1, 1, 1, -1};
        for (size_t id = 0; id < 4; id++) {
            V.row(id) = Eigen::RowVector3d(cube_size / 2 * dx[id], 0, cube_size / 2 * dz[id]);
            V.row(id + 4) = Eigen::RowVector3d(cube_size / 2 * dx[id], -settings_.pad_thickness, cube_size / 2 * dz[id]);
        }

        size_t f_id = 0;
        //cube top
        F.row(f_id++) = Eigen::RowVector3i(0, 1, 2);
        F.row(f_id++) = Eigen::RowVector3i(2, 3, 0);

        //cube bottom
        F.row(f_id++) = Eigen::RowVector3i(4, 7, 6);
        F.row(f_id++) = Eigen::RowVector3i(6, 5, 4);

        //cube lateral
        for (size_t id = 0; id < 4; id++) {
            F.row(f_id++) = Eigen::RowVector3i(id, id + 4, (id + 1) % 4 + 4);
            F.row(f_id++) = Eigen::RowVector3i((id + 1) % 4 + 4, (id + 1) % 4, id);
        }


        //pillar
        for (size_t id = 0; id < settings_.circle_vertices_num; id++) {
            V.row(id + 8) = Eigen::RowVector3d( cos(2.0 * settings_.PI * id / settings_.circle_vertices_num) * settings_.pillar_radius,
                                               -settings_.pad_thickness,
                                               sin(2.0 * settings_.PI * id / settings_.circle_vertices_num) * settings_.pillar_radius );
            V.row(id + settings_.circle_vertices_num + 8) = Eigen::RowVector3d(cos(2.0 * settings_.PI * id / settings_.circle_vertices_num) * settings_.pillar_radius,
                                                              -settings_.pad_thickness - settings_.pillar_length,
                                                              sin(2.0 * settings_.PI * id / settings_.circle_vertices_num) * settings_.pillar_radius);
        }

        // bottom
        for(size_t id = 1; id < settings_.circle_vertices_num - 1; id++)
        {
            F.row(f_id++) = Eigen::RowVector3i(settings_.circle_vertices_num + 8, id + settings_.circle_vertices_num + 8, id + settings_.circle_vertices_num + 9);
        }

        //lateral
        for(size_t id = 0;id < settings_.circle_vertices_num; id++)
        {
            unsigned long u0 = 8 + id;
            unsigned long v0 = 8 + id + settings_.circle_vertices_num;
            unsigned long u1 = 8 + (id + 1) % settings_.circle_vertices_num;
            unsigned long v1 =  8 + (id + 1) % settings_.circle_vertices_num + settings_.circle_vertices_num;
            F.row(f_id++) = Eigen::RowVector3i(u0, u1, v1);
            F.row(f_id++) = Eigen::RowVector3i(v1, v0, u0);
        }
    }
public:

    Settings settings_;
};


#endif //SUPPORTER_SCENE_ORGANIZER_H
