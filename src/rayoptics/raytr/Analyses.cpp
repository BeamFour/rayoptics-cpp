// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
// Java version by Dibyendu Majumdar
// See LICENSE-ray-optics.txt
//
// C++ port of org.redukti.rayoptics.raytr.Analyses
#include "redukti/rayoptics/raytr/Analyses.h"

#include "redukti/Exceptions.h"
#include "redukti/mathlib/FFTUtil.h"
#include "redukti/mathlib/M.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/parax/ParaxTypes.h"
#include "redukti/rayoptics/specs/Field.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/rayoptics/specs/PupilSpec.h"
#include "redukti/rayoptics/raytr/Trace.h"
#include "redukti/rayoptics/raytr/WaveAbr.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace redukti::rayoptics::raytr {

using mathlib::FFTUtil;
using mathlib::Vector2;

std::vector<std::vector<GridItem>> Analyses::trace_ray_grid(
    optical::OpticalModel *opt_model, const TraceGridDef &grid_rng, specs::Field &fld,
    double wvl, double foc, bool append_if_none, const TraceOptions &trace_options_) {
    (void)foc;
    auto trace_options = trace_options_.copy();
    trace_options.apply_vignetting = false;
    auto start = grid_rng.grid_start;
    auto stop = grid_rng.grid_stop;
    auto num = grid_rng.num_rays;
    auto step = (stop.minus(start)).divide(num - 1);
    std::vector<std::vector<GridItem>> grid;
    for (int i = 0; i < num; i++) {
        std::vector<GridItem> grid_row;

        for (int j = 0; j < num; j++) {
            auto pupil = start;
            auto ray_result = Trace::trace_safe(opt_model, pupil, fld, wvl, trace_options);
            auto ray_pkg = ray_result.pkg;
            if (ray_pkg != nullptr) {
                grid_row.push_back(GridItem(pupil, ray_pkg));
            } else {  // ray outside pupil or failed
                if (append_if_none)
                    grid_row.push_back(GridItem(pupil, nullptr));
            }

            start = Vector2(start.x, start.y + step.y);
        }

        grid.push_back(std::move(grid_row));
        start = Vector2(start.x + step.x, grid_rng.grid_start.y);
    }

    return grid;
}

Analyses::OpdGrid Analyses::eval_wavefront(optical::OpticalModel *opt_model, specs::Field &fld,
                                           double wvl, double foc,
                                           std::optional<Vector2> image_pt_2d,
                                           std::optional<Vector2> image_delta, int num_rays,
                                           double value_if_none, double oversize,
                                           const TraceOptions &trace_options) {
    const auto &fod = opt_model->optical_spec->parax_data->fod;
    auto rs = Trace::setup_pupil_coords(opt_model, fld, wvl, foc, image_pt_2d, image_delta);
    auto ref_sphere = rs.ref_sphere;
    auto cr_pkg = rs.chief_ray_pkg;
    fld.chief_ray = std::const_pointer_cast<ChiefRayPkg>(cr_pkg);
    fld.ref_sphere = std::const_pointer_cast<ReferenceSphere>(ref_sphere);

    auto vig_bbox = fld.vignetting_bbox(*opt_model->optical_spec->pupil, oversize);
    TraceGridDef vig_grid_def(Vector2(vig_bbox[0][0], vig_bbox[0][1]),
                              Vector2(vig_bbox[1][0], vig_bbox[1][1]), num_rays);

    auto grid = trace_ray_grid(opt_model, vig_grid_def, fld, wvl, foc, true, trace_options);

    double convert_to_opd = 1 / opt_model->nm_to_sys_units(wvl);

    OpdGrid opd_grid(grid.size());
    for (std::size_t i = 0; i < grid.size(); i++) {
        const auto &row = grid[i];
        opd_grid[i].resize(row.size());
        for (std::size_t j = 0; j < row.size(); j++) {
            const auto &gij = row[j];
            if (gij.ray_pkg != nullptr) {
                double opdelta = WaveAbr::wave_abr_full_calc(fod, fld, wvl, foc, gij.ray_pkg,
                                                             cr_pkg, ref_sphere);
                double opd = convert_to_opd * opdelta;
                opd_grid[i][j] = {gij.pupil.x, gij.pupil.y, opd};
            } else {
                opd_grid[i][j] = {gij.pupil.x, gij.pupil.y, value_if_none};
            }
        }
    }

    return opd_grid;
}

Analyses::WavefrontGridPkg Analyses::trace_wavefront(
    optical::OpticalModel *opt_model, specs::Field &fld, double wvl, double foc,
    std::optional<Vector2> image_pt_2d, std::optional<Vector2> image_delta, int num_rays,
    double oversize, const TraceOptions &trace_options) {
    const auto &fod = opt_model->optical_spec->parax_data->fod;
    auto rs = Trace::setup_pupil_coords(opt_model, fld, wvl, foc, image_pt_2d, image_delta);
    auto ref_sphere = rs.ref_sphere;
    auto cr_pkg = rs.chief_ray_pkg;
    fld.chief_ray = std::const_pointer_cast<ChiefRayPkg>(cr_pkg);
    fld.ref_sphere = std::const_pointer_cast<ReferenceSphere>(ref_sphere);

    auto vig_bbox = fld.vignetting_bbox(*opt_model->optical_spec->pupil, oversize);
    TraceGridDef vig_grid_def(Vector2(vig_bbox[0][0], vig_bbox[0][1]),
                              Vector2(vig_bbox[1][0], vig_bbox[1][1]), num_rays);

    auto grid = trace_ray_grid(opt_model, vig_grid_def, fld, wvl, foc, true, trace_options);

    std::vector<std::vector<std::shared_ptr<const PreOpdPkg>>> upd_grid;
    for (const auto &i : grid) {
        std::vector<std::shared_ptr<const PreOpdPkg>> row;
        for (const auto &gij : i) {
            if (gij.ray_pkg != nullptr) {
                auto pre_opd_pkg = WaveAbr::wave_abr_pre_calc(fod, fld, wvl, foc, gij.ray_pkg,
                                                              cr_pkg, ref_sphere);
                row.push_back(pre_opd_pkg);
            } else {
                row.push_back(nullptr);
            }
        }
        upd_grid.push_back(std::move(row));
    }

    return WavefrontGridPkg(std::move(grid), std::move(upd_grid));
}

Analyses::OpdGrid Analyses::focus_wavefront(optical::OpticalModel *opt_model,
                                            const WavefrontGridPkg &grid_pkg, specs::Field &fld,
                                            double wvl, double foc,
                                            std::optional<Vector2> image_pt_2d,
                                            std::optional<Vector2> image_delta,
                                            double value_if_none) {
    const auto &fod = opt_model->optical_spec->parax_data->fod;
    const auto &grid = grid_pkg.grid;
    const auto &upd_grid = grid_pkg.upd_grid;
    auto rs = Trace::setup_pupil_coords(opt_model, fld, wvl, foc, image_pt_2d, image_delta);
    auto ref_sphere = rs.ref_sphere;
    auto cr_pkg = rs.chief_ray_pkg;
    double convert_to_opd = 1 / opt_model->nm_to_sys_units(wvl);

    OpdGrid refocused_grid(grid.size());
    for (std::size_t i = 0; i < grid.size(); i++) {
        const auto &ig = grid[i];
        const auto &iu = upd_grid[i];
        refocused_grid[i].resize(ig.size());
        for (std::size_t j = 0; j < ig.size(); j++) {
            const auto &gij = ig[j];
            const auto &uij = iu[j];
            if (gij.ray_pkg != nullptr) {
                double opdelta = WaveAbr::wave_abr_calc(fod, fld, wvl, foc, gij.ray_pkg, cr_pkg,
                                                        uij, ref_sphere);
                double opd = convert_to_opd * opdelta;
                refocused_grid[i][j] = {gij.pupil.x, gij.pupil.y, opd};
            } else {
                refocused_grid[i][j] = {gij.pupil.x, gij.pupil.y, value_if_none};
            }
        }
    }

    return refocused_grid;
}

// --- PSF calculation

std::array<int, 3> Analyses::psf_sampling(std::optional<int> n, std::optional<int> n_pupil,
                                          std::optional<int> n_airy) {
    // Python's round() rounds half to even, as Math.rint and std::nearbyint do
    if (!n.has_value()) {
        n = static_cast<int>(std::nearbyint((*n_pupil * *n_airy) / 2.44));
    } else if (!n_pupil.has_value()) {
        n_pupil = static_cast<int>(std::nearbyint(2.44 * *n / *n_airy));
    } else if (!n_airy.has_value()) {
        n_airy = static_cast<int>(std::nearbyint(2.44 * *n / *n_pupil));
    } else {
        // upstream: npa.index(None) raises ValueError
        throw IllegalArgumentException("one of n, n_pupil and n_airy must be null");
    }
    return {*n, *n_pupil, *n_airy};
}

std::array<double, 2> Analyses::calc_psf_scaling(optical::OpticalModel *opt_model,
                                                 const specs::Field &fld, double wvl, int ndim,
                                                 int maxdim) {
    const auto &fod = opt_model->optical_spec->parax_data->fod;
    double wl = opt_model->nm_to_sys_units(wvl);

    double fill_factor = static_cast<double>(ndim) / maxdim;
    double max_D = 2 * fod.enp_radius / fill_factor;
    double delta_x = max_D / maxdim;
    double C = wl / fod.exp_radius;

    double delta_theta = (fill_factor * C) / 2;
    double ref_sphere_radius = fld.ref_sphere->ref_sphere_radius;
    double delta_xp = delta_theta * ref_sphere_radius;

    return {delta_x, delta_xp};
}

std::vector<std::vector<double>> Analyses::calc_psf(
    const std::vector<std::vector<double>> &wavefront, int ndim, int maxdim) {
    auto mdim = static_cast<std::size_t>(maxdim);
    int maxdim_by_2 = maxdim / 2;
    std::vector<std::vector<double>> W(mdim, std::vector<double>(mdim));
    std::vector<std::vector<bool>> pupil(mdim, std::vector<bool>(mdim));
    int nd2 = ndim / 2;
    int lo = maxdim_by_2 - (nd2 - 1);
    int hi = maxdim_by_2 + (nd2 + 1);
    int rows = static_cast<int>(wavefront.size());
    // numpy raises here when the slice and the wavefront differ in shape, which
    // happens for every odd ndim
    if (lo < 0 || hi > maxdim || hi - lo != rows) {
        std::string shape = std::to_string(std::min(hi, maxdim) - std::max(lo, 0));
        throw IllegalArgumentException("could not broadcast input array from shape (" +
                                       std::to_string(rows) + "," + std::to_string(rows) +
                                       ") into shape (" + shape + "," + shape + ")");
    }
    for (int i = lo; i < hi; i++) {
        const auto &row = wavefront[static_cast<std::size_t>(i - lo)];
        if (static_cast<int>(row.size()) != hi - lo)
            throw IllegalArgumentException("wavefront must be square");
        for (int j = lo; j < hi; j++) {
            double w = row[static_cast<std::size_t>(j - lo)];
            W[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = nan_to_num(w);
            pupil[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = std::isfinite(w);
        }
    }

    // phase = np.exp(1j*2*np.pi*W), as interleaved real and imaginary parts by row
    FFTUtil::ComplexArray phase(mdim, std::vector<double>(2 * mdim));
    for (std::size_t i = 0; i < mdim; i++) {
        for (std::size_t j = 0; j < mdim; j++) {
            if (!pupil[i][j]) continue;
            double angle = 2 * mathlib::M::PI * W[i][j];
            phase[i][2 * j] = std::cos(angle);
            phase[i][2 * j + 1] = std::sin(angle);
        }
    }

    auto F = FFTUtil::fftshift(FFTUtil::fft2(FFTUtil::fftshift(phase)));
    std::vector<std::vector<double>> AP(mdim, std::vector<double>(mdim));
    double AP_max = -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < mdim; i++) {
        for (std::size_t j = 0; j < mdim; j++) {
            double a = std::hypot(F[i][2 * j], F[i][2 * j + 1]);
            AP[i][j] = a * a;
            if (!std::isnan(AP[i][j]))
                AP_max = std::max(AP_max, AP[i][j]);
        }
    }
    if (!(AP_max > 0) || !std::isfinite(AP_max))
        throw IllegalArgumentException("wavefront must produce a finite, nonzero PSF");
    for (std::size_t i = 0; i < mdim; i++)
        for (std::size_t j = 0; j < mdim; j++)
            AP[i][j] = AP[i][j] / AP_max;
    return AP;
}

double Analyses::nan_to_num(double x) {
    if (std::isnan(x)) return 0.0;
    if (x == std::numeric_limits<double>::infinity()) return std::numeric_limits<double>::max();
    if (x == -std::numeric_limits<double>::infinity()) return -std::numeric_limits<double>::max();
    return x;
}

} // namespace redukti::rayoptics::raytr
