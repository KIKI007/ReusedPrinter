//
// Created by 汪子琦 on 5/2/17.
//

#ifndef SUPPORTER_MESH_SLICER_OVERHANG_H
#define SUPPORTER_MESH_SLICER_OVERHANG_H
#include "mesh_slicer_base.h"

using ClipperLib::Paths;
using ClipperLib::ClipType;
using ClipperLib::ctIntersection;
using ClipperLib::ctDifference;

class MeshSlicerOverhang : public MeshSlicerBase
{
public:
    MeshSlicerOverhang()
    {

    }

    MeshSlicerOverhang(const MeshSlicerOverhang& slicer)
    :MeshSlicerBase(slicer)
    {
        overhang_layer = slicer.overhang_layer;
        amass_layer = slicer.amass_layer;
    }

public:/* Utility */

    //clean all data in this mesh
    void clear();

public:

    // bottom_half is those support area which touch down to the building platform
    void get_bottom_half_overhang(std::vector< ClipperLib::Paths> &bottom_half);

    // upper_half is those support area which doesn't touch down to the building platform
    void get_upper_half_overhang(std::vector< ClipperLib::Paths> &upper_half);

protected:

    // overhang[i] = difference[ layer_slices[i], layer_slices[i - 1]]
    void overhang_construction();

    // amass[i] = layer_silces[i] + amass[i - 1]
    void amass_construction();

    // construct the bottom_half or upper_half overhang PATHS
    void get_half_overhang(std::vector< ClipperLib::Paths> &half_slices, ClipType type);

protected:

    std::vector<Paths> overhang_layer;

    std::vector<Paths> amass_layer;
};

void MeshSlicerOverhang::clear() {
    MeshSlicerBase::clear();
    overhang_layer.clear();
    amass_layer.clear();
}

void MeshSlicerOverhang::get_bottom_half_overhang(std::vector<ClipperLib::Paths> &bottom_half) {
    get_half_overhang(bottom_half, ctDifference);
}

void MeshSlicerOverhang::get_upper_half_overhang(std::vector<ClipperLib::Paths> &upper_half) {
    get_half_overhang(upper_half, ctIntersection);
}

void MeshSlicerOverhang::overhang_construction() {
    assert(!V.isZero());

    if(layer_slices.empty())
        contour_construction();

    overhang_layer.clear();
    overhang_layer.resize(number_layer());

    settings.print_N();
    settings.print_TsN("OVERHANG");
    for(int layer = 1; layer < number_layer(); layer++)
    {
        memset(settings.tmp_str, 0, sizeof(settings.tmp_str));
        sprintf(settings.tmp_str, "\t layer %d, total %.3f %%...", layer,
                100.0f * (double) (layer) / (layer_slices.size() - 2));
        settings.print_Ts(settings.tmp_str);

        ClipperLib::Clipper clipper;
        ClipperLib::ClipperOffset co;
        ClipperLib::Paths upper = layer_slices[layer], lower = layer_slices[layer - 1];

        co.AddPaths(lower, ClipperLib::jtRound, ClipperLib::etClosedPolygon);
        co.Execute(lower, settings.mm2int(settings.overhang_offset));
        ClipperLib::SimplifyPolygons(lower);

        clipper.StrictlySimple(true);
        clipper.AddPaths(upper, ClipperLib::ptSubject, true);
        clipper.AddPaths(lower, ClipperLib::ptClip, true);
        clipper.Execute(ClipperLib::ctDifference, overhang_layer[layer], ClipperLib::pftPositive, ClipperLib::pftPositive);

        settings.print_TsN("done");
    }

    overhang_layer[0].clear();
}

void MeshSlicerOverhang::amass_construction() {
    if(layer_slices.empty())
        contour_construction();

    amass_layer.clear();
    amass_layer.resize(number_layer());

    amass_layer[0] = layer_slices[0];
    for(int layer = 1; layer < number_layer(); layer++)
    {
        ClipperLib::Clipper clipper;
        clipper.AddPaths(layer_slices[layer], ClipperLib::ptSubject, true);
        clipper.AddPaths(amass_layer[layer - 1], ClipperLib::ptClip, true);
        clipper.Execute(ClipperLib::ctUnion, amass_layer[layer], ClipperLib::pftPositive, ClipperLib::pftPositive);
    }

    return;
}

void MeshSlicerOverhang::get_half_overhang(std::vector<ClipperLib::Paths> &half_slices, ClipType type) {

    assert(!V.isZero());

    if(overhang_layer.empty())
        overhang_construction();

    if(amass_layer.empty())
        amass_construction();

    half_slices.clear();
    half_slices.resize(number_layer());

    settings.print_N();
    settings.print_TsN("OVERHANG");

    if(type == ctDifference)
        half_slices[0] = layer_slices[0];

    for(int layer = 1; layer < number_layer(); layer++)
    {
        memset(settings.tmp_str, 0, sizeof(settings.tmp_str));
        sprintf(settings.tmp_str, "\t layer %d, total %.3f %%...", layer,
                100.0f * (double) (layer) / (layer_slices.size() - 2));
        settings.print_Ts(settings.tmp_str);

        //half
        ClipperLib::Clipper clipper;
        clipper.AddPaths(overhang_layer[layer], ClipperLib::ptSubject, true);
        clipper.AddPaths(amass_layer[layer - 1], ClipperLib::ptClip, true);
        clipper.Execute(type, half_slices[layer], ClipperLib::pftPositive, ClipperLib::pftPositive);

        settings.print_TsN("done");
    }
    return;
}

#endif //SUPPORTER_MESH_SLICER_OVERHANG_H
