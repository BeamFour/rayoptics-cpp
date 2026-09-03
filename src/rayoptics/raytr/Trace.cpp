// C++ port of org.redukti.rayoptics.raytr.Trace
#include "redukti/rayoptics/raytr/Trace.h"

#include "redukti/Exceptions.h"
#include "redukti/Text.h"
#include "redukti/mathlib/M.h"
#include "redukti/mathlib/MinPack.h"
#include "redukti/mathlib/SecantSolver.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/parax/ParaxModel.h"
#include "redukti/rayoptics/parax/ParaxTypes.h"
#include "redukti/rayoptics/raytr/ExitPupilAiming.h"
#include "redukti/rayoptics/raytr/WaveAbr.h"
#include "redukti/rayoptics/raytr/Wideangle.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/rayoptics/util/Lists.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace redukti::rayoptics::raytr {

namespace M = mathlib::M;
using exceptions::TraceException;
using mathlib::MinPack;
using mathlib::Vector2;
using mathlib::Vector3;
using util::Lists::get;

RayDataFrame::RayDataFrame(const std::vector<RaySeg> &raySegList) {
    for (const auto &seg : raySegList) {
        inc_pt.push_back(seg.p);
        after_dir.push_back(seg.d);
        after_dst.push_back(seg.dst);
        normal.push_back(seg.nrml);
    }
}

RayResult Trace::trace_ray(optical::OpticalModel *opt_model, const Vector2 &pupil,
                           specs::Field &fld, double wvl, TraceOptions &trace_options) {
    if (!trace_options.rayerr_filter.has_value())
        trace_options.rayerr_filter = std::string("full");
    return trace_safe(opt_model, pupil, fld, wvl, trace_options);
}

RayResult Trace::trace_safe(optical::OpticalModel *opt_model, const Vector2 &pupil,
                            specs::Field &fld, double wvl,
                            const TraceOptions &trace_options) {
    RayResult result;
    std::shared_ptr<const RayPkg> ray_pkg;
    try {
        auto arr = pupil.as_array();
        ray_pkg = Trace::trace_base(opt_model, std::vector<double>{arr[0], arr[1]}, fld,
                                    wvl, trace_options);
        if (!trace_options.output_filter.has_value())
            result.pkg = ray_pkg;
        else if (*trace_options.output_filter == "last") {
            RaySeg seg = get(ray_pkg->ray, -1);
            ray_pkg = std::make_shared<const RayPkg>(std::vector<RaySeg>{seg},
                                                     ray_pkg->op_delta, ray_pkg->wvl);
            result.pkg = ray_pkg;
        } else {
            throw UnsupportedOperationException();
        }
    } catch (TraceException &rayerr) {
        // The failure is turned into data here: the partial ray comes back out
        // of the exception and the caller carries on.
        if (trace_options.rayerr_filter.has_value() &&
            *trace_options.rayerr_filter == "full") {
            ray_pkg = rayerr.ray_pkg;
            result.pkg = ray_pkg;
            result.err = std::make_shared<TraceException>(rayerr);
        } else if (trace_options.rayerr_filter.has_value() &&
                   *trace_options.rayerr_filter == "summary") {
            rayerr.ray_pkg = nullptr;
            result.err = std::make_shared<TraceException>(rayerr);
            result.pkg = nullptr;
        }
    }
    return result;
}

std::shared_ptr<const RayPkg> Trace::trace(seq::SequentialModel *seq_model,
                                           const Vector3 &pt0, const Vector3 &dir0,
                                           double wvl,
                                           const TraceOptions &trace_options) {
    RayTraceOptions options(trace_options);
    return RayTrace::trace(seq_model, pt0, dir0, wvl, options);
}

std::shared_ptr<const RayPkg> Trace::trace_base(optical::OpticalModel *opt_model,
                                                const std::vector<double> &pupil,
                                                specs::Field &fld, double wvl,
                                                const TraceOptions &trace_options) {
    std::vector<double> pupil_coords = pupil;
    if (trace_options.pupil_type == PupilType::REL_PUPIL) {
        if (trace_options.apply_vignetting)
            pupil_coords = fld.apply_vignetting(pupil);
    }
    specs::Coord coord = opt_model->optical_spec->ray_start_from_osp(
        pupil_coords, fld, trace_options.pupil_type);
    auto pt0 = coord.pt;
    auto dir0 = coord.dir;
    RayTraceOptions options(trace_options);
    if (opt_model->optical_spec->fov->is_wide_angle)
        options.intersect_obj = false;
    else {
        if (dir0.z * util::value(opt_model->seq_model->z_dir[0]) < 0)
            dir0 = dir0.negate();
    }
    auto pkg = RayTrace::trace(opt_model->seq_model.get(), pt0, dir0, wvl, options);
    return pkg->with(&fld, Vector2(pupil[0], pupil[1]),
                     Vector2(pupil_coords[0], pupil_coords[1]));
}

namespace {

/**
 * Shared state for the ray-aiming objective functions. `rr` points at the
 * caller's RayResultWithStartCoord::rr, which the Java mutates through a
 * reference, so the last traced ray is visible to the caller even on failure.
 */
class BaseObjectiveFunction {
public:
    seq::SequentialModel *seq_model;
    std::optional<int> ifcx;
    Vector3 pt0;
    double obj2enp_dist;
    double wvl;
    bool not_wa;
    RayResult *rr;

    BaseObjectiveFunction(seq::SequentialModel *seq_model_, std::optional<int> ifcx_,
                          const Vector3 &pt0_, double obj2enp_dist_, double wvl_,
                          bool not_wa_, RayResult *rr_)
        : seq_model(seq_model_), ifcx(ifcx_), pt0(pt0_), obj2enp_dist(obj2enp_dist_),
          wvl(wvl_), not_wa(not_wa_), rr(rr_) {}

    RaySeg evalSeg(double x1, double y1) {
        Vector3 pt1(x1, y1, obj2enp_dist);
        Vector3 dir0 = pt1.minus(pt0).normalize();
        if (not_wa && dir0.z * util::value(seq_model->z_dir[0]) < 0)
            dir0 = dir0.negate();
        std::shared_ptr<const RayPkg> pkg;
        try {
            pkg = RayTrace::trace(seq_model, pt0, dir0, wvl);
            rr->pkg = pkg;
            rr->err = nullptr;
        } catch (TraceException &ray_error) {
            pkg = ray_error.ray_pkg;
            rr->pkg = ray_error.ray_pkg;
            rr->err = std::make_shared<TraceException>(ray_error);
            if (ray_error.surf <= *ifcx)
                throw;
        }
        return pkg->ray[static_cast<std::size_t>(*ifcx)];
    }
};

class SecantFunction : public BaseObjectiveFunction, public mathlib::ScalarObjectiveFunction {
public:
    double y_target;

    SecantFunction(seq::SequentialModel *seq_model_, std::optional<int> ifcx_,
                   const Vector3 &pt0_, double dist, double wvl_, double y_target_,
                   bool not_wa_, RayResult *rr_)
        : BaseObjectiveFunction(seq_model_, ifcx_, pt0_, dist, wvl_, not_wa_, rr_),
          y_target(y_target_) {}

    std::optional<double> eval(double y1) override {
        RaySeg seg = evalSeg(0., y1);
        double y_ray = seg.p.y;
        return y_ray - y_target;
    }
};

class HybrdObjectiveFunction : public BaseObjectiveFunction, public mathlib::Hybrd_Function {
public:
    std::vector<double> xy_target;

    HybrdObjectiveFunction(seq::SequentialModel *seq_model_, std::optional<int> ifcx_,
                           const Vector3 &pt0_, double dist, double wvl_,
                           std::vector<double> xy_target_, bool not_wa_, RayResult *rr_)
        : BaseObjectiveFunction(seq_model_, ifcx_, pt0_, dist, wvl_, not_wa_, rr_),
          xy_target(std::move(xy_target_)) {}

    void apply(int n, std::vector<double> &x, std::vector<double> &fvec,
               std::vector<int> &iflag) override {
        (void)n;
        (void)iflag;
        RaySeg seg = evalSeg(x[0], x[1]);
        fvec[0] = seg.p.x - xy_target[0];
        fvec[1] = seg.p.y - xy_target[1];
    }
};

std::vector<Vector2> generate_hexapolar_points(const TraceRingsDef &grid_rng,
                                               double max_radius, int num_rings);
std::vector<Vector2> generate_points(const TraceRingsDef &grid_rng, int num_rings,
                                     double max_radius);
Vector2 apply_vignetting_pt(const Vector2 &pupil, const specs::Field &fld);
double vignetting_jacobian(const Vector2 &pupil, const specs::Field &fld);

/** Java's Arrays.toString(double[]). */
std::string arrayToString(const std::vector<double> &v) {
    std::string s = "[";
    for (std::size_t i = 0; i < v.size(); i++) {
        if (i > 0)
            s += ", ";
        s += doubleToString(v[i]);
    }
    return s + "]";
}

} // namespace

RayResultWithStartCoord Trace::get_1d_solution(seq::SequentialModel *seq_model,
                                               std::optional<int> ifcx,
                                               const Vector3 &pt0, double dist,
                                               double wvl, double y_target,
                                               bool not_wa) {
    RayResultWithStartCoord res;
    SecantFunction fn(seq_model, ifcx, pt0, dist, wvl, y_target, not_wa, &res.rr);
    double start_y = mathlib::SecantSolver::find_root(fn, 0., 50, 1.48e-8).root;
    res.start_coords = std::vector<double>{0, start_y};
    return res;
}

RayResultWithStartCoord Trace::get_2d_solution(seq::SequentialModel *seq_model,
                                               std::optional<int> ifcx,
                                               const Vector3 &pt0, double dist,
                                               double wvl,
                                               const std::vector<double> &xy_target,
                                               bool not_wa) {
    RayResultWithStartCoord res;
    HybrdObjectiveFunction f(seq_model, ifcx, pt0, dist, wvl, xy_target, not_wa,
                             &res.rr);
    std::vector<double> x(2, 0.0);
    std::vector<double> fvec(2, 0.0);
    int lwa = (2 * (3 * 2 + 13)) / 2;
    std::vector<double> wa(static_cast<std::size_t>(lwa), 0.0);
    std::vector<int> info(1, 0);
    double epsfcn = 1.0e-8;
    info[0] = MinPack::hybrd1(f, 2, x, fvec, 1.0e-10, wa, lwa, epsfcn);
    std::vector<int> dummy(1, 0);
    f.apply(2, x, fvec, dummy);
    double residual = std::hypot(fvec[0], fvec[1]);
    double coordinateScale = std::max(std::max(std::abs(x[0]), std::abs(x[1])),
                                      std::max(std::abs(xy_target[0]),
                                               std::abs(xy_target[1])));
    double residualTolerance = std::max(1.0e-7, 1.0e-8 * coordinateScale);
    bool converged = info[0] >= 1 && info[0] <= 4 && residual <= residualTolerance;
    if (!converged) {
        TraceException failure("2D ray aiming failed: MINPACK info=" +
                               intToString(info[0]) + ", residual=" +
                               doubleToString(residual) + ", tolerance=" +
                               doubleToString(residualTolerance) + ", start=" +
                               arrayToString(x));
        failure.surf = *ifcx;
        // The Java also sets failure.ifc here; that field is not carried.
        failure.ray_pkg = res.rr.pkg;
        throw failure;
    }
    res.start_coords = x;
    return res;
}

RayResultWithStartCoord Trace::iterate_ray(optical::OpticalModel *opt_model,
                                           std::optional<int> ifcx,
                                           const std::vector<double> &xy_target,
                                           specs::Field &fld, double wvl) {
    auto seq_model = opt_model->seq_model.get();
    auto osp = opt_model->optical_spec.get();
    auto &fod = osp->parax_data->fod;
    double obj2enp_dist = fod.obj_dist + fod.enp_dist;
    bool not_wa = !osp->fov->is_wide_angle;
    specs::Coord coord = osp->obj_coords(fld);
    auto pt0 = coord.pt;
    if (ifcx.has_value()) {
        if (pt0.x == 0.0 && xy_target[0] == 0.0) {
            auto y_target = xy_target[1];
            return get_1d_solution(seq_model, ifcx, pt0, obj2enp_dist, wvl, y_target,
                                   not_wa);
        } else {
            return get_2d_solution(seq_model, ifcx, pt0, obj2enp_dist, wvl, xy_target,
                                   not_wa);
        }
    } else {
        RayResultWithStartCoord result;
        result.start_coords = xy_target;
        return result;
    }
}

std::vector<std::shared_ptr<const RayPkg>> Trace::trace_boundary_rays_at_field(
    optical::OpticalModel *opt_model, specs::Field &fld, double wvl,
    TraceOptions &trace_options) {
    if (!trace_options.rayerr_filter.has_value())
        trace_options.rayerr_filter = std::string("full");
    auto ref_sphere_cr =
        setup_pupil_coords(opt_model, fld, wvl, 0.0, std::nullopt, std::nullopt);
    fld.chief_ray = std::const_pointer_cast<ChiefRayPkg>(ref_sphere_cr.chief_ray_pkg);
    fld.ref_sphere =
        std::const_pointer_cast<ReferenceSphere>(ref_sphere_cr.ref_sphere);
    std::vector<std::shared_ptr<const RayPkg>> rim_rays;
    auto osp = opt_model->optical_spec.get();
    for (const auto &p : osp->pupil->pupil_rays) {
        auto ray_result = trace_ray(opt_model, Vector2(p[0], p[1]), fld, wvl,
                                    trace_options);
        rim_rays.push_back(ray_result.pkg);
    }
    return rim_rays;
}

std::map<std::string, std::shared_ptr<const RayPkg>> Trace::boundary_ray_dict(
    optical::OpticalModel *opt_model,
    const std::vector<std::shared_ptr<const RayPkg>> &rim_rays) {
    std::map<std::string, std::shared_ptr<const RayPkg>> pupil_rays;
    auto &ray_labels = opt_model->optical_spec->pupil->ray_labels;
    for (std::size_t i = 0; i < rim_rays.size(); i++) {
        if (i >= ray_labels.size())
            break;
        pupil_rays[ray_labels[i]] = rim_rays[i];
    }
    return pupil_rays;
}

std::vector<std::vector<std::shared_ptr<const RayPkg>>> Trace::trace_boundary_rays(
    optical::OpticalModel *opt_model, TraceOptions &trace_options) {
    std::vector<std::vector<std::shared_ptr<const RayPkg>>> rayset;
    double wvl = opt_model->seq_model->central_wavelength();
    auto fov = opt_model->optical_spec->fov.get();
    for (std::size_t fi = 0; fi < fov->fields.size(); fi++) {
        specs::Field &fld = *fov->fields[fi];
        auto rim_rays = trace_boundary_rays_at_field(opt_model, fld, wvl, trace_options);
        fld.pupil_rays = boundary_ray_dict(opt_model, rim_rays);
        rayset.push_back(rim_rays);
    }
    return rayset;
}

std::vector<RayDataFrame> Trace::trace_ray_list_at_field(
    optical::OpticalModel *opt_model, const std::vector<std::vector<double>> &ray_list,
    specs::Field &fld, double wvl, double foc, TraceOptions &trace_options) {
    (void)foc;
    std::vector<RayDataFrame> rayset;
    for (const auto &p : ray_list) {
        auto ray_result = trace_ray(opt_model, Vector2(p[0], p[1]), fld, wvl,
                                    trace_options);
        rayset.push_back(RayDataFrame(ray_result.pkg->ray));
    }
    return rayset;
}

RayDataFrameByField Trace::trace_field(optical::OpticalModel *opt_model,
                                       specs::Field &fld, double wvl, double foc) {
    auto osp = opt_model->optical_spec.get();
    auto &pupil_rays = osp->pupil->pupil_rays;
    TraceOptions opts;
    auto rdf_list =
        trace_ray_list_at_field(opt_model, pupil_rays, fld, wvl, foc, opts);
    return RayDataFrameByField(&fld, rdf_list);
}

std::vector<RayDataFrameByField> Trace::trace_all_fields(
    optical::OpticalModel *opt_model) {
    auto osp = opt_model->optical_spec.get();
    auto t = osp->lookup_fld_wvl_focus(0);
    auto wvl = t.second;
    auto foc = t.third;
    std::vector<RayDataFrameByField> fset;
    for (auto &f : osp->fov->fields) {
        auto rset = trace_field(opt_model, *f, wvl, foc);
        fset.push_back(rset);
    }
    return fset;
}

std::shared_ptr<const ChiefRayPkg> Trace::trace_chief_ray(
    optical::OpticalModel *opt_model, specs::Field &fld, double wvl, double foc) {
    (void)foc;
    auto osp = opt_model->optical_spec.get();
    auto &fod = osp->parax_data->fod;
    TraceOptions options;
    options.rayerr_filter = std::string("full");
    auto ray_result = trace_safe(opt_model, Vector2(0., 0.), fld, wvl, options);
    auto cr = ray_result.pkg;
    auto cr_exp_seg = WaveAbr::transfer_to_exit_pupil(
        get(opt_model->seq_model->ifcs, -2),
        RayData(get(cr->ray, -2).p, get(cr->ray, -2).d), fod.exp_dist);
    return std::make_shared<const ChiefRayPkg>(cr, cr_exp_seg);
}

void Trace::apply_paraxial_vignetting(optical::OpticalModel *opt_model) {
    auto fov = opt_model->optical_spec->field_of_view();
    auto pm = opt_model->parax_model.get();
    auto mf = fov->max_field();
    auto max_field = mf.first;
    for (std::size_t j = 0; j < fov->fields.size(); j++) {
        auto &fld = *fov->fields[j];
        auto rel_fov = std::sqrt(fld.x * fld.x + fld.y * fld.y);
        if (!fov->is_relative && max_field != 0)
            rel_fov = rel_fov / max_field;
        auto vg = pm->paraxial_vignetting(rel_fov);
        auto min_vly = vg.first;
        auto min_vuy = vg.second;
        if (min_vly.second.has_value())
            fld.vly = 1.0 - min_vly.first;
        if (min_vuy.second.has_value())
            fld.vuy = 1.0 - min_vuy.first;
    }
}

std::shared_ptr<const ChiefRayPkg> Trace::get_chief_ray_pkg(
    optical::OpticalModel *opt_model, specs::Field &fld, double wvl, double foc) {
    std::shared_ptr<const ChiefRayPkg> chief_ray_pkg;
    if (fld.chief_ray == nullptr) {
        auto res = aim_chief_ray(opt_model, fld, wvl);
        if (!res.aim_pt.empty()) {
            fld.aim_info = res.aim_pt;
            fld.z_enp = std::nullopt;
        } else {
            fld.z_enp = res.z_enp;
            fld.aim_info = std::nullopt;
        }
        chief_ray_pkg = trace_chief_ray(opt_model, fld, wvl, foc);
    } else if (fld.chief_ray->chief_ray->wvl != wvl) {
        chief_ray_pkg = trace_chief_ray(opt_model, fld, wvl, foc);
    } else {
        chief_ray_pkg = fld.chief_ray;
    }
    return chief_ray_pkg;
}

RefSphereCR Trace::setup_pupil_coords(optical::OpticalModel *opt_model,
                                      specs::Field &fld, double wvl, double foc,
                                      std::optional<Vector2> image_pt,
                                      std::optional<Vector2> image_delta) {
    auto chief_ray_pkg = get_chief_ray_pkg(opt_model, fld, wvl, foc);
    auto ref_sphere = WaveAbr::calculate_reference_sphere(
        opt_model, fld, wvl, foc, *chief_ray_pkg, image_pt, image_delta);
    return RefSphereCR(ref_sphere, chief_ray_pkg);
}

AimInfo Trace::aim_chief_ray(optical::OpticalModel *opt_model, specs::Field &fld,
                             std::optional<double> wvl_) {
    auto seq_model = opt_model->seq_model.get();
    auto osp = opt_model->optical_spec.get();
    double wvl = wvl_.has_value() ? *wvl_ : seq_model->central_wavelength();
    std::optional<int> stop = seq_model->stop_surface;
    if (osp->fov->is_wide_angle) {
        auto res = Wideangle::find_real_enp(opt_model, stop, fld, wvl);
        return AimInfo(std::vector<double>{}, res.z_enp);
    } else {
        auto res = iterate_ray(opt_model, stop, std::vector<double>{0., 0.}, fld, wvl);
        return AimInfo(res.start_coords.has_value() ? *res.start_coords
                                                    : std::vector<double>{},
                       std::nullopt);
    }
}

std::vector<GridItem> Trace::trace_fan(optical::OpticalModel *opt_model,
                                       const TraceFanDef &fan_rng, specs::Field &fld,
                                       double wvl, double foc, bool append_if_none,
                                       ImageFilter *img_filter,
                                       const TraceOptions &trace_options) {
    (void)foc;
    auto start = fan_rng.start;
    auto stop = fan_rng.stop;
    auto num = fan_rng.num_rays;
    auto step = (stop.minus(start)).divide(num - 1);
    std::vector<GridItem> fan;
    for (int r = 0; r < num; r++) {
        auto pupil = start;
        auto ray_result = trace_safe(opt_model, pupil, fld, wvl, trace_options);
        if (ray_result.pkg != nullptr) {
            if (img_filter != nullptr) {
                // pkg is non-null on this path, and every filter answers a
                // value for a ray that traced, so this never throws. Java
                // would add a null here and NPE in the first consumer.
                fan.push_back(img_filter->apply(pupil, ray_result.pkg).value());
            } else {
                fan.push_back(GridItem(pupil, ray_result.pkg));
            }
        } else if (append_if_none) {
            fan.push_back(GridItem(pupil, nullptr));
        }
        start = Vector2(start.x + step.x, start.y + step.y);
    }
    return fan;
}

std::vector<GridItem> Trace::trace_grid(optical::OpticalModel *opt_model,
                                        const TraceGridDef &grid_rng, specs::Field &fld,
                                        double wvl, double foc, ImageFilter *img_filter,
                                        bool append_if_none,
                                        const TraceOptions &trace_options_in) {
    (void)foc;
    TraceOptions trace_options = trace_options_in.copy();
    trace_options.check_apertures = true;
    auto start = grid_rng.grid_start;
    auto stop = grid_rng.grid_stop;
    auto num = grid_rng.num_rays;
    auto step = (stop.minus(start)).divide(num - 1);
    std::vector<GridItem> grid;
    for (int i = 0; i < num; i++) {
        for (int j = 0; j < num; j++) {
            auto pupil = start;
            auto ray_result = trace_safe(opt_model, pupil, fld, wvl, trace_options);
            if (ray_result.pkg != nullptr) {
                if (img_filter != nullptr) {
                    grid.push_back(img_filter->apply(pupil, ray_result.pkg).value());
                } else {
                    grid.push_back(GridItem(pupil, ray_result.pkg));
                }
            } else {
                if (img_filter != nullptr) {
                    auto item = img_filter->apply(pupil, nullptr);
                    if (item.has_value())
                        grid.push_back(*item);
                    else if (append_if_none)
                        // Java adds the null itself, which NPEs in the first
                        // consumer; a failed marker is what the equivalent
                        // branch of trace_gaussian_quadrature adds.
                        grid.push_back(GridItem::failed(pupil));
                } else {
                    if (append_if_none)
                        grid.push_back(GridItem(pupil, nullptr));
                }
            }
            start = Vector2(start.x, start.y + step.y);
        }
        start = Vector2(start.x + step.x, grid_rng.grid_start.y);
    }
    return grid;
}

std::vector<GridItem> Trace::trace_rings(optical::OpticalModel *opt_model,
                                         const TraceRingsDef &grid_rng,
                                         specs::Field &fld, double wvl, double foc,
                                         ImageFilter *img_filter, bool append_if_none,
                                         const TraceOptions &trace_options_in) {
    (void)foc;
    TraceOptions trace_options = trace_options_in.copy();
    trace_options.check_apertures = true;
    trace_options.pupil_type = PupilType::REL_PUPIL;
    trace_options.apply_vignetting = true;
    std::vector<GridItem> grid;
    int num_rings = grid_rng.num_rings;
    double max_radius = grid_rng.max_radius;
    std::vector<Vector2> points;
    if (grid_rng.hexapolar) {
        points = generate_hexapolar_points(grid_rng, max_radius, num_rings);
    } else {
        points = generate_points(grid_rng, num_rings, max_radius);
    }
    for (std::size_t i = 0; i < points.size(); i++) {
        auto pupil = points[i];
        auto ray_result = trace_safe(opt_model, pupil, fld, wvl, trace_options);
        if (ray_result.pkg != nullptr) {
            if (img_filter != nullptr) {
                grid.push_back(img_filter->apply(pupil, ray_result.pkg).value());
            } else {
                grid.push_back(GridItem(pupil, ray_result.pkg));
            }
        } else {
            if (img_filter != nullptr) {
                auto item = img_filter->apply(pupil, nullptr);
                if (item.has_value())
                    grid.push_back(*item);
                else if (append_if_none)
                    // See the note in trace_grid.
                    grid.push_back(GridItem::failed(pupil));
            } else {
                if (append_if_none)
                    grid.push_back(GridItem(pupil, nullptr));
            }
        }
    }
    return grid;
}

void Trace::list_ray(std::string &sb, const RayPkg &ray_pkg,
                     const std::optional<math::Tfm3d> &tfrms, std::optional<int> start_) {
    int start = start_.has_value() ? *start_ : 0;
    auto &ray = ray_pkg.ray;
    sb += "            X            Y            Z           L            M            "
          "N               Len\n";
    char buf[256];
    for (std::size_t i = static_cast<std::size_t>(start); i < ray.size(); i++) {
        const auto &r = ray[i];
        Vector3 p = r.p;
        Vector3 d = r.d;
        if (tfrms.has_value()) {
            auto rot = *tfrms->rt;
            auto trns = tfrms->t;
            p = rot.multiply(r.p).plus(trns);
            d = rot.multiply(r.d);
        }
        // "%3d %12.5f %12.5f %12.5g %12.6f %12.6f %12.6f %12.5g\n" -- the two
        // %g fields use Java's semantics, not C's.
        std::snprintf(buf, sizeof(buf), "%3d %12.5f %12.5f ", static_cast<int>(i), p.x,
                      p.y);
        sb += buf;
        sb += formatG(p.z, 12, 5);
        std::snprintf(buf, sizeof(buf), " %12.6f %12.6f %12.6f ", d.x, d.y, d.z);
        sb += buf;
        sb += formatG(r.dst, 12, 5);
        sb += "\n";
    }
}

// ---------------------------------------------------------------------------
// Pupil sampling
// ---------------------------------------------------------------------------

namespace {

std::vector<Vector2> generate_hexapolar_points(const TraceRingsDef &grid_rng,
                                               double max_radius, int num_rings) {
    std::vector<Vector2> points;
    double cx = grid_rng.cx, cy = grid_rng.cy;
    points.push_back(Vector2(cx, cy));
    double step = max_radius / num_rings;
    for (double r = max_radius; r > 1e-8; r -= step) {
        double astep = (step / r) * (M::PI / 3);
        for (double a = 0; a < 2 * M::PI - 1e-8; a += astep)
            points.push_back(Vector2(cx + std::sin(a) * r, cy + std::cos(a) * r));
    }
    return points;
}

std::vector<Vector2> generate_points(const TraceRingsDef &grid_rng, int num_rings,
                                     double max_radius) {
    std::vector<Vector2> points;
    double cx = grid_rng.cx, cy = grid_rng.cy;
    points.push_back(Vector2(cx, cy));
    int num_points_in_ring_one = grid_rng.num_points_in_ring_one;
    double angle_deg_ring_one = 360.0 / num_points_in_ring_one;
    for (int ring = 1; ring <= num_rings; ring++) {
        double daz = angle_deg_ring_one / ring;
        double offset = (ring % 2 == 0) ? 0.0 : 0.5 * daz;
        double r = ring * max_radius / num_rings;
        int numPoints = num_points_in_ring_one * ring;
        for (int jaz = 0; jaz < numPoints; jaz++) {
            double angle_deg = offset + jaz * daz;
            double angle_rad = M::toRadians(angle_deg);
            double x = cx + r * std::cos(angle_rad);
            double y = cy + r * std::sin(angle_rad);
            points.push_back(Vector2(x, y));
        }
    }
    return points;
}

/** Ported for completeness; nothing in the codebase calls it. */
[[maybe_unused]] std::vector<Vector2> generate_gaussian(const TraceRingsDef &grid_rng,
                                                        int ncircles,
                                                        double max_radius) {
    std::vector<Vector2> points;
    double cx = grid_rng.cx, cy = grid_rng.cy;
    points.push_back(Vector2(cx, cy));
    double sigma = max_radius / std::sqrt(2.0 * std::log(1 + ncircles));
    for (int icirc = 1; icirc <= ncircles; icirc++) {
        double daz = 60.0 / icirc;
        double offset = (icirc % 2 == 0) ? 0.0 : 0.5 * daz;
        double p = M::square(icirc) / (ncircles + M::square(ncircles));
        double r = sigma * std::sqrt(2.0 * std::log(1 / (1 - p)));
        for (int jaz = 0; jaz < 6 * icirc; jaz++) {
            double angle_deg = offset + jaz * daz;
            double angle_rad = M::toRadians(angle_deg);
            double x = cx + r * std::cos(angle_rad);
            double y = cy + r * std::sin(angle_rad);
            points.push_back(Vector2(x, y));
        }
    }
    return points;
}

std::vector<std::vector<double>> gauss_legendre_nodes_and_weights(int order) {
    std::vector<double> nodes(static_cast<std::size_t>(order), 0.0);
    std::vector<double> weights(static_cast<std::size_t>(order), 0.0);
    int rootsToFind = (order + 1) / 2;
    for (int i = 0; i < rootsToFind; i++) {
        double x = std::cos(M::PI * (i + 0.75) / (order + 0.5));
        double derivative;
        double delta;
        do {
            double p0 = 1.0;
            double p1 = x;
            for (int degree = 2; degree <= order; degree++) {
                double p2 =
                    ((2.0 * degree - 1.0) * x * p1 - (degree - 1.0) * p0) / degree;
                p0 = p1;
                p1 = p2;
            }
            double polynomial = p1;
            double previousPolynomial = p0;
            derivative = order * (x * polynomial - previousPolynomial) / (x * x - 1.0);
            delta = polynomial / derivative;
            x -= delta;
        } while (std::abs(delta) > 1.0e-15);
        double weight = 2.0 / ((1.0 - x * x) * derivative * derivative);
        nodes[static_cast<std::size_t>(i)] = -x;
        nodes[static_cast<std::size_t>(order - 1 - i)] = x;
        weights[static_cast<std::size_t>(i)] = weight;
        weights[static_cast<std::size_t>(order - 1 - i)] = weight;
    }
    return {nodes, weights};
}

Vector2 apply_vignetting_pt(const Vector2 &pupil, const specs::Field &fld) {
    auto arr = pupil.as_array();
    auto vignetted = fld.apply_vignetting(std::vector<double>{arr[0], arr[1]});
    return Vector2(vignetted[0], vignetted[1]);
}

double vignetting_jacobian(const Vector2 &pupil, const specs::Field &fld) {
    return fld.vignetting_scale_x(pupil.x) * fld.vignetting_scale_y(pupil.y);
}

const double MIN_CONTRAST_CONTRACTION = 0.1;

bool valid_contrast_point(const Vector2 &pupil, const Vector2 &sagittalShift,
                          const Vector2 &tangentialShift, const specs::Field &fld) {
    return Trace::inside_vignetted_pupil(pupil, fld) &&
           Trace::inside_vignetted_pupil(pupil.plus(sagittalShift), fld) &&
           Trace::inside_vignetted_pupil(pupil.plus(tangentialShift), fld);
}

bool valid_contrast_pattern(const std::vector<Vector2> &offsets,
                            const Vector2 &overlapCenter, double contraction,
                            const Vector2 &sagittalShift, const Vector2 &tangentialShift,
                            const specs::Field &fld) {
    for (const auto &offset : offsets) {
        auto pupil = overlapCenter.plus(offset.times(contraction));
        if (!valid_contrast_point(pupil, sagittalShift, tangentialShift, fld)) {
            return false;
        }
    }
    return true;
}

} // namespace

bool Trace::inside_vignetted_pupil(const Vector2 &pupil, const specs::Field &fld) {
    double xScale = fld.vignetting_scale_x(pupil.x);
    double yScale = fld.vignetting_scale_y(pupil.y);
    if (!(xScale > 0.0) || !(yScale > 0.0)) {
        return false;
    }
    double x = pupil.x / xScale;
    double y = pupil.y / yScale;
    return x * x + y * y <= 1.0 + 1.0e-14;
}

std::vector<GaussianQuadraturePoint> Trace::generate_gaussian_quadrature(
    const TraceRingsDef &grid_rng, int num_rings, std::optional<int> num_spokes) {
    if (num_rings < 1 || (num_spokes.has_value() && *num_spokes < 3)) {
        throw IllegalArgumentException(
            "The number of rings must be at least 1 and spokes must be at least 3");
    }
    if (!std::isfinite(grid_rng.min_radius) || !std::isfinite(grid_rng.max_radius) ||
        grid_rng.min_radius < 0.0 || grid_rng.max_radius <= grid_rng.min_radius) {
        throw IllegalArgumentException(
            "Pupil radii must be finite and satisfy 0 <= min_radius < max_radius");
    }
    int spokes = !num_spokes.has_value() ? 4 * (num_rings + 1) : *num_spokes;
    auto nodesAndWeights = gauss_legendre_nodes_and_weights(num_rings);
    std::vector<GaussianQuadraturePoint> points;
    points.reserve(static_cast<std::size_t>(num_rings * spokes));
    for (int angle = 1; angle <= spokes; angle++) {
        double theta = 2.0 * M::PI * angle / spokes;
        double cosTheta = std::cos(theta);
        double sinTheta = std::sin(theta);
        for (int ring = 0; ring < num_rings; ring++) {
            double radialFraction =
                0.5 + 0.5 * nodesAndWeights[0][static_cast<std::size_t>(ring)];
            double innerRadiusSquared = grid_rng.min_radius * grid_rng.min_radius;
            double outerRadiusSquared = grid_rng.max_radius * grid_rng.max_radius;
            double radius = std::sqrt(innerRadiusSquared +
                                      radialFraction *
                                          (outerRadiusSquared - innerRadiusSquared));
            Vector2 pupil(grid_rng.cx + radius * cosTheta,
                          grid_rng.cy + radius * sinTheta);
            double weight =
                0.5 * nodesAndWeights[1][static_cast<std::size_t>(ring)] / spokes;
            points.push_back(GaussianQuadraturePoint(pupil, weight));
        }
    }
    return points;
}

std::vector<GaussianQuadraturePoint> Trace::generate_contrast_quadrature(
    const TraceRingsDef &grid_rng, std::optional<int> num_spokes,
    const Vector2 &sagittal_shift, const Vector2 &tangential_shift, specs::Field &fld) {
    auto nominal = generate_gaussian_quadrature(grid_rng, grid_rng.num_rings, num_spokes);
    Vector2 nominalCenter(grid_rng.cx, grid_rng.cy);
    auto physicalCenter = apply_vignetting_pt(nominalCenter, fld);
    auto overlapCenter =
        physicalCenter.minus(sagittal_shift.plus(tangential_shift).times(0.5));
    if (!valid_contrast_point(overlapCenter, sagittal_shift, tangential_shift, fld)) {
        throw IllegalArgumentException(
            "The requested contrast shear has no common vignetted pupil overlap");
    }
    std::vector<Vector2> offsets;
    offsets.reserve(nominal.size());
    for (const auto &point : nominal)
        offsets.push_back(apply_vignetting_pt(point.pupil, fld).minus(physicalCenter));
    double contraction = 1.0;
    if (!valid_contrast_pattern(offsets, overlapCenter, contraction, sagittal_shift,
                                tangential_shift, fld)) {
        double low = 0.0;
        double high = 1.0;
        for (int iteration = 0; iteration < 60; iteration++) {
            double trial = 0.5 * (low + high);
            if (valid_contrast_pattern(offsets, overlapCenter, trial, sagittal_shift,
                                       tangential_shift, fld)) {
                low = trial;
            } else {
                high = trial;
            }
        }
        contraction = low;
    }
    if (contraction < MIN_CONTRAST_CONTRACTION) {
        throw IllegalArgumentException(
            "The requested contrast shear leaves too little vignetted pupil overlap to "
            "sample: the pattern would contract to " +
            doubleToString(contraction) + " of its nominal extent");
    }
    std::vector<GaussianQuadraturePoint> points;
    points.reserve(nominal.size());
    double weightSum = 0.0;
    for (std::size_t i = 0; i < nominal.size(); i++) {
        const auto &point = nominal[i];
        auto pupil = overlapCenter.plus(offsets[i].times(contraction));
        double weight = point.weight * vignetting_jacobian(point.pupil, fld);
        points.push_back(GaussianQuadraturePoint(pupil, weight));
        weightSum += weight;
    }
    if (!(weightSum > 0.0)) {
        throw IllegalArgumentException("The vignetted pupil has zero area");
    }
    const double normalization = weightSum;
    std::vector<GaussianQuadraturePoint> normalized;
    normalized.reserve(points.size());
    for (const auto &point : points)
        normalized.push_back(
            GaussianQuadraturePoint(point.pupil, point.weight / normalization));
    return normalized;
}

std::vector<GridItem> Trace::trace_gaussian_quadrature(
    optical::OpticalModel *opt_model, const TraceRingsDef &grid_rng,
    std::optional<int> num_spokes, specs::Field &fld, double wvl, double foc,
    ImageFilter *img_filter, bool append_if_none, const TraceOptions &trace_options_in) {
    (void)foc;
    TraceOptions trace_options = trace_options_in.copy();
    trace_options.pupil_type = PupilType::REL_PUPIL;
    trace_options.apply_vignetting = true;
    std::vector<GridItem> grid;
    auto points = generate_gaussian_quadrature(grid_rng, grid_rng.num_rings, num_spokes);
    for (const auto &point : points) {
        auto pupil = point.pupil;
        auto ray_result = trace_safe(opt_model, pupil, fld, wvl, trace_options);
        if (ray_result.pkg != nullptr) {
            GridItem item = img_filter != nullptr
                                ? img_filter->apply(pupil, ray_result.pkg).value()
                                : GridItem(pupil, ray_result.pkg);
            grid.push_back(item.withWeight(point.weight));
        } else if (img_filter != nullptr) {
            auto item = img_filter->apply(pupil, nullptr);
            if (item.has_value())
                grid.push_back(item->withWeight(point.weight));
            else if (append_if_none)
                grid.push_back(GridItem::failed(pupil).withWeight(point.weight));
        } else if (append_if_none) {
            grid.push_back(GridItem::failed(pupil).withWeight(point.weight));
        }
    }
    return grid;
}

std::vector<ContrastRayTriplet> Trace::trace_contrast(
    optical::OpticalModel *opt_model, const TraceRingsDef &grid_rng,
    std::optional<int> num_spokes, const Vector2 &sagittal_shift,
    const Vector2 &tangential_shift, specs::Field &fld, double wvl,
    const TraceOptions &trace_options) {
    return trace_contrast(opt_model, grid_rng, num_spokes, sagittal_shift,
                          tangential_shift, std::nullopt, std::nullopt, fld, wvl,
                          trace_options, false);
}

std::vector<ContrastRayTriplet> Trace::trace_contrast(
    optical::OpticalModel *opt_model, const TraceRingsDef &grid_rng,
    std::optional<int> num_spokes, const Vector2 &sagittal_shift,
    const Vector2 &tangential_shift, std::optional<Vector2> sagittal_exit_shift,
    std::optional<Vector2> tangential_exit_shift, specs::Field &fld, double wvl,
    const TraceOptions &trace_options_in, bool aim_exit_pupil) {
    if (aim_exit_pupil &&
        (!sagittal_exit_shift.has_value() || !tangential_exit_shift.has_value())) {
        throw IllegalArgumentException(
            "Physical exit-pupil shifts are required when aiming is enabled");
    }
    TraceOptions trace_options = trace_options_in.copy();
    trace_options.pupil_type = PupilType::REL_PUPIL;
    trace_options.apply_vignetting = false;
    trace_options.rayerr_filter = std::string("summary");
    std::vector<ContrastRayTriplet> samples;
    auto points = generate_contrast_quadrature(grid_rng, num_spokes, sagittal_shift,
                                               tangential_shift, fld);
    for (const auto &point : points) {
        auto pupil = point.pupil;
        auto sagittalPupil = pupil.plus(sagittal_shift);
        auto tangentialPupil = pupil.plus(tangential_shift);
        auto reference = trace_safe(opt_model, pupil, fld, wvl, trace_options);
        RayResult sagittal;
        RayResult tangential;
        if (aim_exit_pupil && reference.pkg != nullptr) {
            auto referenceCoordinate = ExitPupilAiming::sphere_coord(
                reference.pkg, fld.chief_ray, fld.ref_sphere);
            if (!referenceCoordinate.has_value()) {
                auto error = std::make_shared<ExitPupilAiming::ExitPupilAimException>(
                    "Reference ray has no finite exit-pupil coordinate");
                error->surf = -1;
                sagittal = RayResult(nullptr, error);
                tangential = RayResult(nullptr, error);
            } else {
                Vector2 sagittalTarget(referenceCoordinate->x + sagittal_exit_shift->x,
                                       referenceCoordinate->y + sagittal_exit_shift->y);
                Vector2 tangentialTarget(
                    referenceCoordinate->x + tangential_exit_shift->x,
                    referenceCoordinate->y + tangential_exit_shift->y);
                sagittal = ExitPupilAiming::aim(opt_model, sagittalPupil, sagittalTarget,
                                                fld, wvl, trace_options)
                               .ray;
                tangential = ExitPupilAiming::aim(opt_model, tangentialPupil,
                                                  tangentialTarget, fld, wvl,
                                                  trace_options)
                                 .ray;
            }
        } else {
            sagittal = trace_safe(opt_model, sagittalPupil, fld, wvl, trace_options);
            tangential = trace_safe(opt_model, tangentialPupil, fld, wvl, trace_options);
        }
        samples.push_back(ContrastRayTriplet(pupil, reference.pkg, sagittal.pkg,
                                             tangential.pkg, reference.err,
                                             sagittal.err, tangential.err,
                                             point.weight));
    }
    return samples;
}

// ---------------------------------------------------------------------------
// "raw" variants, iterating over an explicit path rather than a model
// ---------------------------------------------------------------------------

namespace {

class BaseObjectiveFunctionRaw {
public:
    const std::vector<seq::PathSeg> *pthlist;
    std::optional<int> ifcx;
    Vector3 pt0;
    double dist;
    double wvl;
    bool not_wa;
    RayResult *rr;

    BaseObjectiveFunctionRaw(const std::vector<seq::PathSeg> *pthlist_,
                             std::optional<int> ifcx_, const Vector3 &pt0_, double dist_,
                             double wvl_, bool not_wa_, RayResult *rr_)
        : pthlist(pthlist_), ifcx(ifcx_), pt0(pt0_), dist(dist_), wvl(wvl_),
          not_wa(not_wa_), rr(rr_) {}

    RaySeg evalSeg(double x1, double y1) {
        Vector3 pt1(x1, y1, dist);
        Vector3 dir0 = pt1.minus(pt0).normalize();
        if (not_wa && dir0.z * util::value(*(*pthlist)[0].Zdir) < 0)
            dir0 = dir0.negate();
        std::shared_ptr<const RayPkg> pkg;
        try {
            RayTraceOptions options;
            options.check_apertures = false;
            options.intersect_obj = true;
            options.filter_out_phantoms = false;
            pkg = RayTrace::trace_raw(*pthlist, pt0, dir0, wvl, options);
            rr->pkg = pkg;
            rr->err = nullptr;
        } catch (TraceException &ray_error) {
            pkg = ray_error.ray_pkg;
            rr->pkg = ray_error.ray_pkg;
            rr->err = std::make_shared<TraceException>(ray_error);
            if (ray_error.surf <= *ifcx)
                throw;
        }
        return pkg->ray[static_cast<std::size_t>(*ifcx)];
    }
};

class SecantFunctionRaw : public BaseObjectiveFunctionRaw,
                          public mathlib::ScalarObjectiveFunction {
public:
    double y_target;

    SecantFunctionRaw(const std::vector<seq::PathSeg> *pthlist_, std::optional<int> ifcx_,
                      const Vector3 &pt0_, double dist_, double wvl_, double y_target_,
                      bool not_wa_, RayResult *rr_)
        : BaseObjectiveFunctionRaw(pthlist_, ifcx_, pt0_, dist_, wvl_, not_wa_, rr_),
          y_target(y_target_) {}

    std::optional<double> eval(double y1) override {
        RaySeg seg = evalSeg(0., y1);
        double y_ray = seg.p.y;
        return y_ray - y_target;
    }
};

class HybrdObjectiveFunctionRaw : public BaseObjectiveFunctionRaw,
                                  public mathlib::Hybrd_Function {
public:
    std::vector<double> xy_target;

    HybrdObjectiveFunctionRaw(const std::vector<seq::PathSeg> *pthlist_,
                              std::optional<int> ifcx_, const Vector3 &pt0_, double dist_,
                              double wvl_, std::vector<double> xy_target_, bool not_wa_,
                              RayResult *rr_)
        : BaseObjectiveFunctionRaw(pthlist_, ifcx_, pt0_, dist_, wvl_, not_wa_, rr_),
          xy_target(std::move(xy_target_)) {}

    void apply(int n, std::vector<double> &x, std::vector<double> &fvec,
               std::vector<int> &iflag) override {
        (void)n;
        (void)iflag;
        RaySeg seg = evalSeg(x[0], x[1]);
        fvec[0] = seg.p.x - xy_target[0];
        fvec[1] = seg.p.y - xy_target[1];
    }
};

} // namespace

RayResultWithStartCoord Trace::get_1d_solution_raw(
    const std::vector<seq::PathSeg> &pthlist, std::optional<int> ifcx,
    const Vector3 &pt0, double dist, double wvl, double y_target, bool not_wa) {
    RayResultWithStartCoord res;
    SecantFunctionRaw fn(&pthlist, ifcx, pt0, dist, wvl, y_target, not_wa, &res.rr);
    double start_y = mathlib::SecantSolver::find_root(fn, 0., 50, 1.48e-8).root;
    res.start_coords = std::vector<double>{0, start_y};
    return res;
}

RayResultWithStartCoord Trace::get_2d_solution_raw(
    const std::vector<seq::PathSeg> &pthlist, std::optional<int> ifcx,
    const Vector3 &pt0, double dist, double wvl, const std::vector<double> &xy_target,
    bool not_wa) {
    RayResultWithStartCoord res;
    HybrdObjectiveFunctionRaw f(&pthlist, ifcx, pt0, dist, wvl, xy_target, not_wa,
                                &res.rr);
    std::vector<double> x(2, 0.0);
    std::vector<double> fvec(2, 0.0);
    int lwa = (2 * (3 * 2 + 13)) / 2;
    std::vector<double> wa(static_cast<std::size_t>(lwa), 0.0);
    std::vector<int> info(1, 0);
    double epsfcn = 1.0e-8;
    info[0] = MinPack::hybrd1(f, 2, x, fvec, 1.0e-10, wa, lwa, epsfcn);
    std::vector<int> dummy(1, 0);
    f.apply(2, x, fvec, dummy);
    double residual = std::hypot(fvec[0], fvec[1]);
    double coordinateScale = std::max(std::max(std::abs(x[0]), std::abs(x[1])),
                                      std::max(std::abs(xy_target[0]),
                                               std::abs(xy_target[1])));
    double residualTolerance = std::max(1.0e-7, 1.0e-8 * coordinateScale);
    bool converged = info[0] >= 1 && info[0] <= 4 && residual <= residualTolerance;
    if (!converged) {
        TraceException failure("2D ray aiming failed: MINPACK info=" +
                               intToString(info[0]) + ", residual=" +
                               doubleToString(residual) + ", tolerance=" +
                               doubleToString(residualTolerance) + ", start=" +
                               arrayToString(x));
        failure.surf = *ifcx;
        failure.ray_pkg = res.rr.pkg;
        throw failure;
    }
    res.start_coords = x;
    return res;
}

RayResultWithStartCoord Trace::iterate_ray_raw(
    const std::vector<seq::PathSeg> &pthlist, std::optional<int> ifcx,
    const std::vector<double> &xy_target, const Vector3 &pt0, const Vector3 &d0,
    double obj2pup_dist, double eprad, double wvl, bool not_wa) {
    (void)d0;
    (void)eprad;
    if (ifcx.has_value()) {
        if (pt0.x == 0.0 && xy_target[0] == 0.0) {
            auto y_target = xy_target[1];
            try {
                return get_1d_solution_raw(pthlist, ifcx, pt0, obj2pup_dist, wvl,
                                           y_target, not_wa);
            } catch (TraceException &ray_err) {
                RayResultWithStartCoord result;
                result.start_coords = std::vector<double>{0, 0};
                result.rr = RayResult(ray_err.ray_pkg,
                                      std::make_shared<TraceException>(ray_err));
                return result;
            }
        } else {
            try {
                return get_2d_solution_raw(pthlist, ifcx, pt0, obj2pup_dist, wvl,
                                           xy_target, not_wa);
            } catch (TraceException &ray_err) {
                RayResultWithStartCoord result;
                result.start_coords = std::vector<double>{0, 0};
                result.rr = RayResult(ray_err.ray_pkg,
                                      std::make_shared<TraceException>(ray_err));
                return result;
            }
        }
    } else {
        RayResultWithStartCoord result;
        result.start_coords = xy_target;
        return result;
    }
}

} // namespace redukti::rayoptics::raytr
