// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
// Java version by Dibyendu Majumdar
// See LICENSE-ray-optics.txt
//
// C++ port of org.redukti.rayoptics.raytr.RayGrid
#include "redukti/rayoptics/raytr/RayGrid.h"

#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/specs/FieldSpec.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/rayoptics/specs/WvlSpec.h"

#include <limits>

namespace redukti::rayoptics::raytr {

RayGrid::RayGrid(optical::OpticalModel *opt_model_, int f, std::optional<double> wl,
                 std::optional<double> foc_, std::optional<mathlib::Vector2> image_pt_2d_,
                 std::optional<mathlib::Vector2> image_delta_,
                 std::optional<std::string> output_filter,
                 std::optional<std::string> rayerr_filter, int num_rays_, bool clip_rays,
                 double value_if_none_, double oversize_)
    : RayGrid(opt_model_, *opt_model_->optical_spec->fov->fields[static_cast<std::size_t>(f)],
              wl, foc_, image_pt_2d_, image_delta_, std::move(output_filter),
              std::move(rayerr_filter), num_rays_, clip_rays, value_if_none_, oversize_) {}

RayGrid::RayGrid(optical::OpticalModel *opt_model_, specs::Field &f, std::optional<double> wl,
                 std::optional<double> foc_, std::optional<mathlib::Vector2> image_pt_2d_,
                 std::optional<mathlib::Vector2> image_delta_,
                 std::optional<std::string> output_filter,
                 std::optional<std::string> rayerr_filter, int num_rays_, bool clip_rays,
                 double value_if_none_, double oversize_) {
    this->opt_model = opt_model_;
    auto osp = opt_model_->optical_spec.get();
    this->fld = &f;
    this->wvl = !wl.has_value() ? osp->wvls->central_wvl() : *wl;

    this->foc = !foc_.has_value() ? osp->defocus()->focus_shift : *foc_;
    this->image_pt_2d = image_pt_2d_;
    this->image_delta = image_delta_;

    this->num_rays = num_rays_;
    this->value_if_none = value_if_none_;

    this->oversize = oversize_;
    this->rt_kwargs.output_filter = std::move(output_filter);
    this->rt_kwargs.rayerr_filter = std::move(rayerr_filter);
    this->rt_kwargs.check_apertures = clip_rays;

    update_data(true);
}

RayGrid::RayGrid(optical::OpticalModel *opt_model_)
    : RayGrid(opt_model_, 0, std::nullopt, std::nullopt, std::nullopt, std::nullopt,
              std::nullopt, std::nullopt, 21, true, std::numeric_limits<double>::quiet_NaN(),
              1.) {}

RayGrid &RayGrid::update_data(bool rebuild) {
    if (rebuild) {
        grid_pkg = Analyses::trace_wavefront(
            opt_model, *fld, wvl, foc,
            image_pt_2d, image_delta,
            num_rays, oversize,
            rt_kwargs);
    }

    auto opd = Analyses::focus_wavefront(opt_model, *grid_pkg,
                                         *fld, wvl, foc,
                                         image_pt_2d, image_delta,
                                         value_if_none);

    // np.rollaxis(opd, 2)
    std::vector<std::vector<std::vector<double>>> rolled(3);
    for (std::size_t k = 0; k < 3; k++) {
        rolled[k].resize(opd.size());
        for (std::size_t i = 0; i < opd.size(); i++) {
            rolled[k][i].resize(opd[i].size());
            for (std::size_t j = 0; j < opd[i].size(); j++)
                rolled[k][i][j] = opd[i][j][k];
        }
    }
    grid = std::move(rolled);

    return *this;
}

std::array<double, 2> RayGrid::calc_psf_scaling(int maxdim) const {
    return Analyses::calc_psf_scaling(opt_model, *fld, wvl, num_rays, maxdim);
}

} // namespace redukti::rayoptics::raytr
