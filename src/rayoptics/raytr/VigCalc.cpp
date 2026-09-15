// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
// Java version by Dibyendu Majumdar
// See LICENSE-ray-optics.txt
//
// C++ port of org.redukti.rayoptics.raytr.VigCalc
#include "redukti/rayoptics/raytr/VigCalc.h"

#include "redukti/Exceptions.h"
#include "redukti/mathlib/SecantSolver.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/parax/ParaxTypes.h"
#include "redukti/rayoptics/raytr/Trace.h"
#include "redukti/rayoptics/raytr/Wideangle.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/rayoptics/util/Lists.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace redukti::rayoptics::raytr {

using exceptions::TraceException;
using exceptions::TraceMissedSurfaceException;
using exceptions::TraceRayBlockedException;
using mathlib::Vector2;
using specs::ImageKey;
using specs::ValueKey;
using util::Lists::get;

std::optional<double> VigCalc::max_aperture_at_surf(
    const std::vector<std::vector<std::shared_ptr<const RayPkg>>> &rayset, int i) {
    double max_ap = -1.0e+10;
    for (const auto &f : rayset) {
        for (const auto &p : f) {
            const auto &ray = p->ray;
            if (static_cast<int>(ray.size()) > i) {
                auto idx = static_cast<std::size_t>(i);
                auto ap = std::sqrt(ray[idx].p.x * ray[idx].p.x +
                                    ray[idx].p.y * ray[idx].p.y);
                if (ap > max_ap)
                    max_ap = ap;
            } else
                return std::nullopt;
        }
    }
    return max_ap;
}

void VigCalc::set_clear_apertures(optical::OpticalModel *opt_model,
                                  const std::vector<int> *avoid_list,
                                  const std::vector<int> *include_list_in) {
    auto sm = opt_model->seq_model.get();
    auto num_surfs = sm->get_num_surfaces();
    std::vector<int> include_list;
    if (avoid_list == nullptr) {
        if (include_list_in == nullptr) {
            for (int i = 0; i < num_surfs; i++)
                include_list.push_back(i);
        } else {
            include_list = *include_list_in;
        }
    } else {
        for (int i = 0; i < num_surfs; i++) {
            if (std::find(avoid_list->begin(), avoid_list->end(), i) ==
                avoid_list->end())
                include_list.push_back(i);
        }
    }
    TraceOptions opts;
    auto rayset = Trace::trace_boundary_rays(opt_model, opts);
    auto stop_surf = sm->stop_surface;
    if (stop_surf.has_value() && std::find(include_list.begin(), include_list.end(),
                                           *stop_surf) != include_list.end()) {
        std::vector<std::vector<std::shared_ptr<const RayPkg>>> firstOnly{rayset[0]};
        auto max_ap = max_aperture_at_surf(firstOnly, *stop_surf);
        if (max_ap.has_value())
            sm->ifcs[static_cast<std::size_t>(*stop_surf)]->set_max_aperture(*max_ap);
    }
    for (auto i : include_list) {
        if (!(stop_surf.has_value() && i == *stop_surf)) {
            auto max_ap = max_aperture_at_surf(rayset, i);
            if (max_ap.has_value())
                sm->ifcs[static_cast<std::size_t>(i)]->set_max_aperture(*max_ap);
        }
    }
}

void VigCalc::set_ape(optical::OpticalModel *opm, const std::vector<int> *avoid_list,
                      const std::vector<int> *include_list) {
    set_clear_apertures(opm, avoid_list, include_list);
    // Upstream follows this with em.sync_to_seq(sm) to push the new
    // apertures into the element model. Beam43 needs no equivalent: its
    // layout reads Interface.max_aperture and surface_od() live at render
    // time rather than caching them, so there is nothing to go stale.
}

void VigCalc::set_vig(optical::OpticalModel *opm, std::optional<bool> use_bisection) {
    auto osp = opm->optical_spec.get();
    for (std::size_t fi = 0; fi < osp->fov->fields.size(); fi++) {
        auto fld_wvl_foc = osp->lookup_fld_wvl_focus(static_cast<int>(fi));
        auto fld = fld_wvl_foc.first;
        auto wvl = fld_wvl_foc.second;
        calc_vignetting_for_field(opm, *fld, wvl, use_bisection, std::nullopt);
    }
}

/**
 * Set the aperture on the stop surface to satisfy the pupil spec.
 *
 * The vignetting is recalculated after the stop aperture change.
 */
void VigCalc::set_stop_aperture(optical::OpticalModel *opm) {
    auto sm = opm->seq_model.get();
    opm->optical_spec->fov->with_index_label("axis")->clear_vignetting();
    std::vector<int> include{*sm->stop_surface};
    set_clear_apertures(opm, nullptr, &include);
    set_vig(opm, false);
}

/**
 * From existing stop size, calculate pupil spec and vignetting.
 *
 *     Use the upper Y marginal ray on-axis (field #0) and iterate until it
 *     goes through the edge of the stop surface. Use the object or image
 *     segments of this ray to update the pupil specification value
 *     e.g. EPD, NA or f/#.
 */
void VigCalc::set_pupil(optical::OpticalModel *opm, bool use_parax) {
    auto sm = opm->seq_model.get();
    if (!sm->stop_surface.has_value()) {
        std::fprintf(stderr, "Floating stop surface\n");
        return;
    }
    auto idx_stop = *sm->stop_surface;
    auto osp = opm->optical_spec.get();
    // iterate the on-axis marginal ray thru the edge of the stop.
    auto fld_foc = osp->lookup_fld_wvl_focus(0);
    auto fld_0 = fld_foc.first;
    auto cwl = fld_foc.second;
    auto stop_radius = get(sm->ifcs, idx_stop)->surface_od();
    auto start_coords = iterate_pupil_ray(opm, sm->stop_surface, 1, 1.0, stop_radius,
                                          *fld_0, cwl);
    // trace the real axial marginal ray
    TraceOptions options;
    options.output_filter = std::nullopt;
    options.rayerr_filter = std::string("full");
    options.apply_vignetting = false;
    options.check_apertures = false;
    auto ray_result = Trace::trace_safe(opm, start_coords, *fld_0, cwl, options);
    auto ray_pkg = ray_result.pkg;
    // The Java dereferences ray_pkg below without checking, so a trial geometry
    // that kills the axial ray raises NullPointerException there -- caught by
    // the optimizer's merit function and turned into a rejected step. Here that
    // would be a null dereference, so raise the ported equivalent instead.
    if (ray_pkg == nullptr)
        throw IllegalStateException(
            "the axial ray did not trace, so the pupil cannot be set");
    auto obj_img_key = osp->pupil->key.imageKey;
    auto pupil_spec = osp->pupil->key.valueKey;
    auto pupil_value_orig = osp->pupil->value;
    auto parax_data = opm->optical_spec->parax_data;
    auto &ax_ray = parax_data->ax_ray;
    auto &fod = parax_data->fod;
    if (use_parax) {
        auto scale_ratio = stop_radius / ax_ray[static_cast<std::size_t>(idx_stop)].ht;
        //logger.debug(f"{scale_ratio=:8.5f} (parax)")
        if (obj_img_key == ImageKey::Object) {
            if (pupil_spec == ValueKey::EPD) {
                osp->pupil->value = scale_ratio * (2 * fod.enp_radius);
            } else {
                auto slp0 = scale_ratio * ax_ray[0].slp;
                if (pupil_spec == ValueKey::NA) {
                    auto n0 = sm->central_rndx(0);
                    auto rs0 = get(ray_pkg->ray, 0);
                    osp->pupil->value = n0 * rs0.d.y;
                } else if (pupil_spec == ValueKey::Fnum) {
                    osp->pupil->value = 1.0 / (2.0 * slp0);
                }
            }
        } else if (obj_img_key == ImageKey::Image) {
            if (pupil_spec == ValueKey::EPD) {
                osp->pupil->value = scale_ratio * (2 * fod.exp_radius);
            } else {
                auto slpk = scale_ratio * get(ax_ray, -1).slp;
                if (pupil_spec == ValueKey::NA) {
                    auto nk = sm->central_rndx(-1);
                    auto rsm2 = get(ray_pkg->ray, -2);
                    osp->pupil->value = -nk * rsm2.d.y;
                } else if (pupil_spec == ValueKey::Fnum) {
                    osp->pupil->value = -1.0 / (2.0 * slpk);
                }
            }
        }
    } else {
        auto scale_ratio = get(ray_pkg->ray, 1).p.y / get(ax_ray, 1).ht;
        if (obj_img_key == ImageKey::Object) {
            if (pupil_spec == ValueKey::EPD) {
                osp->pupil->value *= scale_ratio;
            } else {
                auto rs0 = get(ray_pkg->ray, 0);
                auto slp0 = rs0.d.y / rs0.d.z;
                if (pupil_spec == ValueKey::NA) {
                    auto n0 = sm->central_rndx(0);
                    osp->pupil->value = n0 * rs0.d.y;
                } else if (pupil_spec == ValueKey::Fnum) {
                    osp->pupil->value = 1.0 / (2.0 * slp0);
                }
            }
        } else if (obj_img_key == ImageKey::Image) {
            auto rsm2 = get(ray_pkg->ray, -2);
            if (pupil_spec == ValueKey::EPD) {
                auto ht = rsm2.p.y;
                osp->pupil->value = 2.0 * ht;
            } else {
                auto slpk = scale_ratio * get(ax_ray, -1).slp;
                if (pupil_spec == ValueKey::NA) {
                    auto nk = sm->central_rndx(-1);
                    osp->pupil->value = -nk * rsm2.d.y;
                } else if (pupil_spec == ValueKey::Fnum) {
                    osp->pupil->value = -1.0 / (2.0 * slpk);
                }
            }
        }
    }
    // trace the real axial marginal ray with aperture clipping
    TraceOptions clipoptions;
    clipoptions.output_filter = std::nullopt;
    clipoptions.rayerr_filter = std::string("full");
    clipoptions.apply_vignetting = false;
    clipoptions.check_apertures = true;
    auto clipped_rr = Trace::trace_safe(opm, start_coords, *fld_0, cwl, clipoptions);
    auto clipped_ray_err = clipped_rr.err;
    if (clipped_ray_err != nullptr) {
        if (dynamic_cast<TraceRayBlockedException *>(clipped_ray_err.get()) != nullptr)
            std::fprintf(stderr,
                         "Axial bundle limited by surface %d not stop surface.\n",
                         clipped_ray_err->surf);
    }
    if (osp->pupil->value != pupil_value_orig) {
        opm->update_model();
    }
    // Always establish vignetting, even when the pupil value was already
    // correct. Skipping is only safe when earlier factors can stand in, and
    // set_pupil is called on freshly built models where there are none -
    // the factors would silently stay at zero. Callers cannot tell that
    // apart from a genuinely unvignetted system: contrast sampling in
    // particular would then take the full pupil as available and, since it
    // traces without aperture checking, optimize light the lens blocks.
    set_vig(opm, std::nullopt);
}

void VigCalc::calc_vignetting_for_field(optical::OpticalModel *opm, specs::Field &fld,
                                        double wvl, std::optional<bool> use_bisection_,
                                        std::optional<int> max_iter_count) {
    bool use_bisection = use_bisection_.has_value() ? *use_bisection_ : false;
    auto &pupil_starts = opm->optical_spec->pupil->pupil_rays;
    double vig_factors[4];
    for (int i = 0; i < 4; i++) {
        int xy = i / 2;
        auto &start = pupil_starts[static_cast<std::size_t>(i + 1)];
        Vector2 startv(start[0], start[1]);
        VigResult result(0.0, std::nullopt, nullptr);
        if (use_bisection) {
            result = calc_vignetted_ray_by_bisection(opm, xy, startv, fld, wvl,
                                                     max_iter_count);
        } else {
            result = calc_vignetted_ray(opm, xy, startv, fld, wvl, max_iter_count);
        }
        vig_factors[i] = result.vig;
    }
    // update the field's vignetting factors
    fld.vux = vig_factors[0];
    fld.vlx = vig_factors[1];
    fld.vuy = vig_factors[2];
    fld.vly = vig_factors[3];
}

std::optional<double> VigCalc::Fn_r_pupil_coordinate::eval(double xy_coord) {
    auto rel_p1 = Vector2::vector2_0.set(xy, xy_coord);
    std::shared_ptr<const RayPkg> ray_pkg;
    try {
        TraceOptions options;
        options.apply_vignetting = false;
        options.check_apertures = false;
        auto arr = rel_p1.as_array();
        ray_pkg = Trace::trace_base(opt_model, std::vector<double>{arr[0], arr[1]}, *fld,
                                    wvl, options);
    } catch (TraceException &ray_error) {
        ray_pkg = ray_error.ray_pkg;
        // Check if the ray error occurred at or before the indx surface.
        // if the error is at or following indx, drop thru
        if (dynamic_cast<TraceMissedSurfaceException *>(&ray_error) != nullptr) {
            // no surface intersection, so no ray data at indx
            if (ray_error.surf <= indx)
                return std::nullopt;
        } else {
            // other ray trace error exceptions
            if (ray_error.surf < indx)
                return std::nullopt;
        }
    }
    // compute the radial distance to the intersection point
    auto p = get(ray_pkg->ray, indx).p;
    auto r_ray = std::copysign(std::sqrt(p.x * p.x + p.y * p.y), r_target);
    auto delta = r_ray - r_target;
    return delta;
}

namespace {

/**
 * Same as Fn_r_pupil_coordinate but rethrows with the pupil coordinate
 * attached, which iterate_pupil_ray reads to fall back on. The Java has these
 * as two separate classes.
 */
class R_Pupil_Coordinate : public mathlib::ScalarObjectiveFunction {
public:
    optical::OpticalModel *opt_model;
    int indx;
    int xy;
    specs::Field *fld;
    double wvl;
    double r_target;

    R_Pupil_Coordinate(optical::OpticalModel *opt_model_, int indx_, int xy_,
                       specs::Field *fld_, double wvl_, double r_target_)
        : opt_model(opt_model_), indx(indx_), xy(xy_), fld(fld_), wvl(wvl_),
          r_target(r_target_) {}

    std::optional<double> eval(double xy_coord) override {
        auto rel_p1 = Vector2::vector2_0.set(xy, xy_coord);
        std::shared_ptr<const RayPkg> ray_pkg;
        try {
            TraceOptions options;
            options.apply_vignetting = false;
            options.check_apertures = false;
            auto arr = rel_p1.as_array();
            ray_pkg = Trace::trace_base(opt_model, std::vector<double>{arr[0], arr[1]},
                                        *fld, wvl, options);
        } catch (TraceException &ray_err) {
            ray_pkg = ray_err.ray_pkg;
            if (dynamic_cast<TraceMissedSurfaceException *>(&ray_err) != nullptr) {
                if (ray_err.surf <= indx) {
                    ray_err.rel_p1 = rel_p1;
                    throw;
                }
            } else if (ray_err.surf < indx) {
                ray_err.rel_p1 = rel_p1;
                throw;
            }
        }
        // compute the radial distance to the intersection point
        auto p = get(ray_pkg->ray, indx).p;
        auto r_ray = std::copysign(std::sqrt(p.x * p.x + p.y * p.y), r_target);
        auto delta = r_ray - r_target;
        //            logger.debug(f"  {xy_coord=:8.5f}   {r_ray=:8.5f}    "
        //                    f"delta={delta:9.2g}")
        //System.out.println(String.format("   xy_coord=%8.5f   r_ray=%8.5f   delta=%9.2g",xy_coord,r_ray,delta));
        return delta;
    }
};

} // namespace

/**
 * Find the limiting aperture and return the vignetting factor.
 *
 *     Args:
 *         opm: :class:`~.OpticalModel` instance
 *         xy: 0 or 1 depending on x or y axis as the pupil direction
 *         start_dir: the unit length starting pupil coordinates, e.g [1., 0.].
 *                    This establishes the radial direction of the ray iteration.
 *         fld: :class:`~.Field` point for wave aberration calculation
 *         wvl: wavelength of ray (nm)
 *         max_iter_count: fail-safe limit on aperture search
 *
 *     Returns:
 *         (**vig**, **clip_indx**, **ray_pkg**)
 *
 *         - **vig** - vignetting factor
 *         - **clip_indx** - the index of the limiting interface
 *         - **ray_pkg** - the vignetting-limited ray
 */
VigResult VigCalc::calc_vignetted_ray(optical::OpticalModel *opm, int xy,
                                      const Vector2 &start_dir, specs::Field &fld,
                                      double wvl, std::optional<int> max_iter_count_) {
    int max_iter_count = max_iter_count_.has_value() ? *max_iter_count_ : 50;
    auto rel_p1 = start_dir;
    auto sm = opm->seq_model.get();
    auto still_iterating = true;
    std::optional<int> clip_indx;
    std::optional<int> stop_indx;
    auto iter_count = 0;
    std::shared_ptr<const RayPkg> ray_pkg;
    while (still_iterating && iter_count < max_iter_count) {
        iter_count++;
        try {
            TraceOptions options;
            options.apply_vignetting = false;
            options.check_apertures = true;
            options.pt_inside_fuzz = 1e-4;
            auto arr = rel_p1.as_array();
            ray_pkg = Trace::trace_base(opm, std::vector<double>{arr[0], arr[1]}, fld,
                                        wvl, options);
            //  ray successfully traced.
            if (clip_indx.has_value()) {
                // fall through and exit
                // The Java computes r_error here and discards it; the call to
                // edge_pt_target is kept because it is the only other effect.
                (void)get(sm->ifcs, *clip_indx)->edge_pt_target(start_dir);
                //                    logger.debug(f" C {xy_str[xy]} = {rel_p1[xy]:10.6f}:   "
                //                            f"blocked at {clip_indx}, del={r_error:8.1e}, "
                //                            "exiting")
                still_iterating = false;
            } else {
                // this is the first time through
                // iterate to find the ray that goes through the edge
                // of the stop surface
                std::optional<int> indx;
                indx = stop_indx = sm->stop_surface;
                if (stop_indx.has_value()) {
                    auto r_target = get(sm->ifcs, *stop_indx)->edge_pt_target(start_dir);
                    //                        logger.debug(f" D {xy_str[xy]} = {rel_p1[xy]:10.6f}:   "
                    //                                f"passed first time, iterate to edge of stop, "
                    //                                f"ifcs[{stop_indx}]")
                    rel_p1 = iterate_pupil_ray(opm, indx, xy, rel_p1.v(xy),
                                               r_target.v(xy), fld, wvl);
                    still_iterating = true;
                    clip_indx = indx;
                } else
                    still_iterating = false;
            }
        } catch (TraceException &ray_error) {
            ray_pkg = ray_error.ray_pkg;
            std::optional<int> indx = ray_error.surf;
            if (clip_indx.has_value() && *indx == *clip_indx) {
                // As above: the Java's r_error computation here is dead, and its
                // IndexOutOfBoundsException catch guarded only that.
                (void)get(sm->ifcs, *clip_indx)->edge_pt_target(start_dir);
                //                        logger.debug(f" A {xy_str[xy]} = {rel_p1[xy]:10.6f}:   "
                //                                f"blocked at {clip_indx}, del={r_error:8.1e}, "
                //                                "exiting")
                //                        logger.debug(f" A' {xy_str[xy]} = {rel_p1[xy]:10.6f}:   "
                //                                f"blocked at {clip_indx}, "
                //                                "exiting")
                still_iterating = false;
            } else {
                auto r_target = get(sm->ifcs, *indx)->edge_pt_target(start_dir);
                // If we missed the first surface, use bisection to bracket
                // the edge. Use the result to start the newton iteration to
                // quickly find the edge.
                if (dynamic_cast<TraceMissedSurfaceException *>(&ray_error) != nullptr) {
                    Fn_r_pupil_coordinate fn(opm, *indx, xy, &fld, wvl, r_target.v(xy));
                    auto edge = Wideangle::find_edge(fn, 0.0, rel_p1.v(xy), std::nullopt);
                    rel_p1 = rel_p1.set(xy, edge.z_enp);
                }
                rel_p1 = iterate_pupil_ray(opm, indx, xy, rel_p1.v(xy), r_target.v(xy),
                                           fld, wvl);
                still_iterating = true;
                clip_indx = indx;
            }
        }
    }
    auto vig = 1.0 - (rel_p1.v(xy) / start_dir.v(xy));
    //        logger.info(f" ray: ({start_dir[0]:2.0f}, {start_dir[1]:2.0f}), "
    //                f"vig={vig:8.4f}, limited at ifcs[{clip_indx}]")
    return VigResult(vig, clip_indx, ray_pkg);
}

/**
 * Find the limiting aperture and return the vignetting factor.
 *
 *     Args:
 *         opm: :class:`~.OpticalModel` instance
 *         xy: 0 or 1 depending on x or y axis as the pupil direction
 *         start_dir: the unit length starting pupil coordinates, e.g [1., 0.].
 *                    This establishes the radial direction of the ray iteration.
 *         fld: :class:`~.Field` point for wave aberration calculation
 *         wvl: wavelength of ray (nm)
 *         max_iter_count: fail-safe limit on aperture search
 *
 *     Returns:
 *         (**vig**, **clip_indx**, **ray_pkg**)
 *
 *         - **vig** - vignetting factor
 *         - **clip_indx** - the index of the limiting interface
 *         - **ray_pkg** - the vignetting-limited ray
 */
VigResult VigCalc::calc_vignetted_ray_by_bisection(optical::OpticalModel *opm, int xy,
                                                   const Vector2 &start_dir,
                                                   specs::Field &fld, double wvl,
                                                   std::optional<int> max_iter_count_) {
    //        logger.debug(f"fld={fld.yf:5.2f}, [{start_dir[0]:5.2f}, "
    //                f"{start_dir[1]:5.2f}]")
    int max_iter_count = max_iter_count_.has_value() ? *max_iter_count_ : 10;
    auto rel_p1 = start_dir;
    std::optional<int> clip_indx;
    auto iter_count = 0;
    auto step_size = 1.0;
    std::shared_ptr<const RayPkg> ray_pkg;
    while (iter_count < max_iter_count) {
        iter_count++;
        try {
            step_size /= 2.0;
            TraceOptions options;
            options.apply_vignetting = false;
            options.check_apertures = true;
            options.pt_inside_fuzz = 1e-4;
            auto arr = rel_p1.as_array();
            ray_pkg = Trace::trace_base(opm, std::vector<double>{arr[0], arr[1]}, fld,
                                        wvl, options);
            rel_p1 = start_dir.times(step_size).plus(rel_p1);
        } catch (TraceException &ray_error) {
            ray_pkg = ray_error.ray_pkg;
            clip_indx = ray_error.surf;
            rel_p1 = start_dir.times(-step_size).plus(rel_p1);
            //                logger.debug(f"{xy_str[xy]} = {rel_p1[xy]:10.6f}: "
            //                        f"blocked at {clip_indx}")
        }
    }
    auto vig = 1.0 - (rel_p1.v(xy) / start_dir.v(xy));
    //        logger.debug(f"   {vig=:7.4f}, {clip_indx=}")
    return VigResult(vig, clip_indx, ray_pkg);
}

Vector2 VigCalc::iterate_pupil_ray(optical::OpticalModel *opt_model,
                                   std::optional<int> indx, int xy, double start_r0,
                                   double r_target, specs::Field &fld, double wvl) {
    Vector2 start_coord = Vector2::vector2_0;
    double start_r = 0;
    if (indx.has_value()) {
        R_Pupil_Coordinate objective_fn(opt_model, *indx, xy, &fld, wvl, r_target);
        try {
            start_r =
                mathlib::SecantSolver::find_root(objective_fn, start_r0, 50, 1e-6).root;
        } catch (TraceException &rt_err) {
            //                logger.debug(f"  {type(rt_err).__name__}: surf={rt_err.surf}    "
            //                        f"rel_p1={rt_err.rel_p1[xy]=:8.5f}   ")
            start_r = 0.9 * rt_err.rel_p1->v(xy);
        }
        return start_coord.set(xy, start_r);
    } else
        return start_coord.set(xy, r_target);
}

} // namespace redukti::rayoptics::raytr
