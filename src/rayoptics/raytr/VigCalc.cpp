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

void VigCalc::set_stop_aperture(optical::OpticalModel *opm) {
    auto sm = opm->seq_model.get();
    opm->optical_spec->fov->with_index_label("axis")->clear_vignetting();
    std::vector<int> include{*sm->stop_surface};
    set_clear_apertures(opm, nullptr, &include);
    set_vig(opm, false);
}

void VigCalc::set_pupil(optical::OpticalModel *opm, bool use_parax) {
    auto sm = opm->seq_model.get();
    if (!sm->stop_surface.has_value()) {
        std::fprintf(stderr, "Floating stop surface\n");
        return;
    }
    auto idx_stop = *sm->stop_surface;
    auto osp = opm->optical_spec.get();
    auto fld_foc = osp->lookup_fld_wvl_focus(0);
    auto fld_0 = fld_foc.first;
    auto cwl = fld_foc.second;
    auto stop_radius = get(sm->ifcs, idx_stop)->surface_od();
    auto start_coords = iterate_pupil_ray(opm, sm->stop_surface, 1, 1.0, stop_radius,
                                          *fld_0, cwl);
    TraceOptions options;
    options.output_filter = std::nullopt;
    options.rayerr_filter = std::string("full");
    options.apply_vignetting = false;
    options.check_apertures = false;
    auto ray_result = Trace::trace_safe(opm, start_coords, *fld_0, cwl, options);
    auto ray_pkg = ray_result.pkg;
    auto obj_img_key = osp->pupil->key.imageKey;
    auto pupil_spec = osp->pupil->key.valueKey;
    auto pupil_value_orig = osp->pupil->value;
    auto parax_data = opm->optical_spec->parax_data;
    auto &ax_ray = parax_data->ax_ray;
    auto &fod = parax_data->fod;
    if (use_parax) {
        auto scale_ratio = stop_radius / ax_ray[static_cast<std::size_t>(idx_stop)].ht;
        if (obj_img_key == ImageKey::Object) {
            if (pupil_spec == ValueKey::EPD) {
                osp->pupil->value = scale_ratio * (2 * fod.enp_radius);
            } else {
                auto slp0 = scale_ratio * ax_ray[0].slp;
                if (pupil_spec == ValueKey::NA) {
                    auto n0 = sm->central_rndx(0);
                    auto rs0 = ray_pkg->ray[0];
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
        auto scale_ratio = ray_pkg->ray[1].p.y / ax_ray[1].ht;
        if (obj_img_key == ImageKey::Object) {
            if (pupil_spec == ValueKey::EPD) {
                osp->pupil->value *= scale_ratio;
            } else {
                auto rs0 = ray_pkg->ray[0];
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
        if (dynamic_cast<TraceMissedSurfaceException *>(&ray_error) != nullptr) {
            if (ray_error.surf <= indx)
                return std::nullopt;
        } else {
            if (ray_error.surf < indx)
                return std::nullopt;
        }
    }
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
        auto p = get(ray_pkg->ray, indx).p;
        auto r_ray = std::copysign(std::sqrt(p.x * p.x + p.y * p.y), r_target);
        auto delta = r_ray - r_target;
        return delta;
    }
};

} // namespace

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
            if (clip_indx.has_value()) {
                // The Java computes r_error here and discards it; the call to
                // edge_pt_target is kept because it is the only other effect.
                (void)get(sm->ifcs, *clip_indx)->edge_pt_target(start_dir);
                still_iterating = false;
            } else {
                std::optional<int> indx;
                indx = stop_indx = sm->stop_surface;
                if (stop_indx.has_value()) {
                    auto r_target = get(sm->ifcs, *stop_indx)->edge_pt_target(start_dir);
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
                still_iterating = false;
            } else {
                auto r_target = get(sm->ifcs, *indx)->edge_pt_target(start_dir);
                if (dynamic_cast<TraceMissedSurfaceException *>(&ray_error) != nullptr &&
                    ray_error.surf == 1) {
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
    return VigResult(vig, clip_indx, ray_pkg);
}

VigResult VigCalc::calc_vignetted_ray_by_bisection(optical::OpticalModel *opm, int xy,
                                                   const Vector2 &start_dir,
                                                   specs::Field &fld, double wvl,
                                                   std::optional<int> max_iter_count_) {
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
        }
    }
    auto vig = 1.0 - (rel_p1.v(xy) / start_dir.v(xy));
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
            start_r = 0.9 * rt_err.rel_p1->v(xy);
        }
        return start_coord.set(xy, start_r);
    } else
        return start_coord.set(xy, r_target);
}

} // namespace redukti::rayoptics::raytr
