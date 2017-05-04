//
// Created by 汪子琦 on 5/3/17.
//

#ifndef SUPPORTER_MESH_SLICER_SHIFT_H
#define SUPPORTER_MESH_SLICER_SHIFT_H

#include "mesh_slicer_overhang.h"
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::RowVector3i;
using Eigen::RowVector3d;
using Eigen::MatrixXd;
using Eigen::Matrix2d;
using std::vector;
using ClipperLib::IntPoint;

class MeshSlicerShift : public MeshSlicerOverhang
{
public:
    MeshSlicerShift();

public:

    void move_xz(double dx, double dz);

    void rotate_yaxis(double angle, Vector2d center);

protected:

    void move_xz(IntPoint &p, int dx, int dz);

    void move_xz(RowVector3d &p, int dx, int dz);

    void move_xz(vector<Paths> &slices, int dx, int dz);

    void rotate_yaxis(IntPoint &p, double angle, Vector2d center);

    void rotate_yaxis(RowVector3d &p, double angle, Vector2d center);

    void rotate_yaxis(vector<Paths> &slices, double angle, Vector2d center);
};

MeshSlicerShift::MeshSlicerShift() {

}

void MeshSlicerShift::move_xz(double dx, double dz) {


    for(int ir = 0; ir < V.rows(); ir++)
    {
        RowVector3d  pv = RowVector3d(V(ir, 0), V(ir, 1), V(ir, 2));
        move_xz(pv, settings.mm2int(dx), settings.mm2int(dz));
        V.row(ir) = RowVector3i(std::round(pv(0)), std::round(pv(1)), std::round(pv(2)));
    }

    move_xz(layer_slices,   settings.mm2int(dx), settings.mm2int(dz));
    move_xz(amass_layer,    settings.mm2int(dx), settings.mm2int(dz));
    move_xz(overhang_layer, settings.mm2int(dx), settings.mm2int(dz));

    return;
}

void MeshSlicerShift::rotate_yaxis(double angle, Vector2d center) {

    for(int ir = 0; ir < V.rows(); ir++)
    {
        RowVector3d  pv = RowVector3d(V(ir, 0), V(ir, 1), V(ir, 2));
        rotate_yaxis(pv, angle ,center);
        V.row(ir) = RowVector3i(std::round(pv(0)), std::round(pv(1)), std::round(pv(2)));
    }

    rotate_yaxis(layer_slices,   angle, center);
    rotate_yaxis(amass_layer,    angle, center);
    rotate_yaxis(overhang_layer, angle, center);
    return;
}

void MeshSlicerShift::move_xz(vector<Paths> &slices, int dx, int dz)
{
    for(int layer = 0; layer < slices.size(); layer++)
    {
        for(int pth = 0; pth < slices[layer].size(); pth++)
        {
            for(int id = 0; id < slices[layer][pth].size(); id++)
            {
                move_xz(slices[layer][pth][id], dx, dz);
            }
        }
    }
}

void MeshSlicerShift::rotate_yaxis(vector<Paths> &slices, double angle, Vector2d center) {
    for(int layer = 0; layer < slices.size(); layer++)
    {
        for(int pth = 0; pth < slices[layer].size(); pth++)
        {
            for(int id = 0; id < slices[layer][pth].size(); id++)
            {
                rotate_yaxis(slices[layer][pth][id], angle, center);
            }
        }
    }
}

void MeshSlicerShift::move_xz(IntPoint &p, int dx, int dz)
{
    p.X = std::round(p.X + dx);
    p.Y = std::round(p.Y + dz);
    return;
}

void MeshSlicerShift::move_xz(RowVector3d &p, int dx, int dz)
{
    p(0) += dx;
    p(2) += dz;
    return;
}

void MeshSlicerShift::rotate_yaxis(RowVector3d &p, double angle, Vector2d center)
{
    Vector2d pt(p(0), p(2));

    Matrix2d rot_mat;
    rot_mat <<  cos(angle), -sin(angle),
                sin(angle),  cos(angle);
    pt -= center;
    pt = rot_mat * pt;
    pt += center;

    p(0) = pt(0);
    p(2) = pt(1);
}

void MeshSlicerShift::rotate_yaxis(IntPoint &p, double angle, Vector2d center)
{
    Vector2d pt(p.X, p.Y);

    Matrix2d rot_mat;
    rot_mat <<  cos(angle), -sin(angle),
                sin(angle),  cos(angle);

    pt -= center;
    pt = rot_mat * pt;
    pt += center;

    p.X = std::round(pt(0));
    p.Y = std::round(pt(1));
}

#endif //SUPPORTER_MESH_SLICER_SHIFT_H
