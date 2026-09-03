// C++ port of org.redukti.rayoptics.raytr.RayTrace
#include "redukti/rayoptics/raytr/RayTrace.h"

#include "redukti/rayoptics/exceptions/TraceException.h"
#include "redukti/rayoptics/seq/SequentialModel.h"

#include <cmath>
#include <optional>

namespace redukti::rayoptics::raytr {

using exceptions::TraceMissedSurfaceException;
using exceptions::TraceRayBlockedException;
using exceptions::TraceTIRException;
using mathlib::Matrix3;
using mathlib::Vector3;
using seq::InteractMode;
using seq::PathSeg;
using util::ZDir;

mathlib::Vector3 RayTrace::bend(const Vector3 &d_in, const Vector3 &normal, double n_in,
                                double n_out) {
    double normal_len = normal.length();
    double cosI = d_in.dot(normal) / normal_len;
    double sinI_sqr = 1.0 - cosI * cosI;
    double sqrrt_in = n_out * n_out - n_in * n_in * sinI_sqr;
    if (sqrrt_in <= 0)
        throw TraceTIRException();
    double sqrrt = std::sqrt(sqrrt_in);
    double n_cosIp = cosI > 0 ? sqrrt : -sqrrt;
    double alpha = n_cosIp - n_in * cosI;
    Vector3 d_out = (d_in.times(n_in).plus(normal.times(alpha))).divide(n_out);
    return d_out;
}

mathlib::Vector3 RayTrace::reflect(const Vector3 &d_in, const Vector3 &normal) {
    double normal_len = normal.length();
    double cosI = d_in.dot(normal) / normal_len;
    Vector3 d_out = d_in.minus(normal.times(2.0 * cosI));
    return d_out;
}

double RayTrace::eic_distance_from_axis(const Vector3 &p, const Vector3 &d, ZDir z_dir) {
    // eq 3.20/3.21
    double e = ((p.dot(d) + util::value(z_dir) * p.z) / (1.0 + util::value(z_dir) * d.z));
    return e;
}

std::shared_ptr<const RayPkg> RayTrace::trace(seq::SequentialModel *seq_model,
                                              const Vector3 &pt0, const Vector3 &dir0,
                                              double wvl, RayTraceOptions &options) {
    // Bind by reference: path() now returns a reference into the model's path
    // cache, and copying it here would give back the per-trace vector copy the
    // cache exists to remove.
    const auto &path = seq_model->path(wvl, std::nullopt, std::nullopt, 1);
    if (!options.first_surf.has_value())
        options.first_surf = 1;
    if (!options.last_surf.has_value())
        options.last_surf = seq_model->get_num_surfaces() - 2;
    return trace_raw(path, pt0, dir0, wvl, options);
}

std::shared_ptr<const RayPkg> RayTrace::trace(seq::SequentialModel *seq_model,
                                              const Vector3 &pt0, const Vector3 &dir0,
                                              double wvl) {
    RayTraceOptions options;
    options.first_surf = 1;
    options.last_surf = seq_model->get_num_surfaces() - 2;
    return trace(seq_model, pt0, dir0, wvl, options);
}

namespace {

bool in_gap_range(int first_surf, std::optional<int> last_surf, int gap_indx,
                  bool include_last_surf) {
    if (last_surf.has_value() && first_surf == *last_surf)
        return false;
    if (gap_indx < first_surf)
        return false;
    if (!last_surf.has_value())
        return true;
    else
        return include_last_surf ? gap_indx <= *last_surf : gap_indx < *last_surf;
}

bool in_surface_range(int first_surf, std::optional<int> last_surf, int s) {
    if (s < first_surf)
        return false;
    if (!last_surf.has_value())
        return true;
    else if (s > *last_surf)
        return false;
    else
        return true;
}

} // namespace

std::shared_ptr<const RayPkg> RayTrace::trace_raw(const std::vector<PathSeg> &path,
                                                  const Vector3 &pt0,
                                                  const Vector3 &dir0, double wvl,
                                                  const RayTraceOptions &options) {
    int first_surf = options.first_surf.has_value() ? *options.first_surf : 0;
    std::optional<int> last_surf = options.last_surf;
    std::vector<RaySeg> ray;
    std::size_t it = 0;
    // trace object surface
    const PathSeg *obj = &path[it++];
    const PathSeg *before = obj;
    InteractMode b4_interact_mode = InteractMode::DUMMY;
    Vector3 before_pt = Vector3::ZERO;
    Vector3 before_normal = Vector3::ZERO;
    if (options.intersect_obj) {
        auto srf_obj = obj->ifc;
        b4_interact_mode = srf_obj->interact_mode;
        auto intersection = srf_obj->intersect(pt0, dir0, options.eps, *obj->Zdir);
        before_pt = intersection.intersection_point;
        before_normal = srf_obj->normal(before_pt);
    } else {
        before_pt = pt0;
        before_normal = Vector3::vector3_001;
    }
    Vector3 before_dir = dir0;
    math::Tfm3d tfrm_from_before = *before->Tfrm;
    // Nullable, exactly as in the Java: path() zips 13 interfaces against only
    // 12 z_dir entries, so the final segment's Zdir is null. Java never unboxes
    // it there (it is only assigned to z_dir_before, and the loop then ends),
    // so this must stay an optional and only be dereferenced where Java unboxes.
    std::optional<ZDir> z_dir_before = before->Zdir;
    double op_delta = 0.0;
    double opl = 0.0;
    int surf = 0;
    Vector3 inc_pt = Vector3::ZERO;
    Vector3 after_dir = Vector3::ZERO;
    Vector3 normal = Vector3::ZERO;
    // loop of remaining surfaces in path
    while (true) {
        double pp_dst = 0.0;
        try {
            // Java relies on Iterator.next() throwing NoSuchElementException to
            // end the loop; the bound is checked directly here and the same
            // terminating work is done below.
            if (it >= path.size()) {
                ray.push_back(RaySeg(inc_pt, after_dir, 0.0, normal));
                op_delta += opl;
                break;
            }
            const PathSeg &after = path[it++];
            surf += 1;
            Matrix3 rt = *tfrm_from_before.rt;
            Vector3 t = tfrm_from_before.t;
            Vector3 b4_pt = rt.multiply(before_pt.minus(t));
            Vector3 b4_dir = rt.multiply(before_dir);
            pp_dst = -b4_pt.dot(b4_dir);
            Vector3 pp_pt_before = b4_pt.plus(b4_dir.times(pp_dst));
            auto ifc = after.ifc;
            InteractMode interact_mode = ifc->interact_mode;
            std::optional<ZDir> z_dir_after = after.Zdir;
            // intersect ray with profile
            auto intersection =
                ifc->intersect(pp_pt_before, b4_dir, options.eps, *z_dir_before);
            double pp_dst_intrsct = intersection.distance;
            inc_pt = intersection.intersection_point;
            double dst_b4 = pp_dst + pp_dst_intrsct;
            if (b4_interact_mode == InteractMode::PHANTOM && options.filter_out_phantoms) {
                // if a phantom interface, don't add intersection point
                // but do add the path length.
                auto pos = ray.size() - 1;
                auto prev_seg = ray[pos];
                ray[pos] = RaySeg(prev_seg, dst_b4);
            } else {
                // add *previous* intersection point, direction, etc., to ray
                ray.push_back(RaySeg(before_pt, before_dir, dst_b4, before_normal));
            }
            if (in_gap_range(first_surf, last_surf, surf - 1, false))
                opl += *before->Indx * dst_b4;
            normal = ifc->normal(inc_pt);
            if (options.check_apertures &&
                in_surface_range(first_surf, last_surf, surf) &&
                interact_mode != InteractMode::PHANTOM) {
                if (!ifc->point_inside(inc_pt.x, inc_pt.y, options.pt_inside_fuzz))
                    throw TraceRayBlockedException(inc_pt);
            }
            // TODO phase element: if present, use it to calculate after_dir
            // refract or reflect ray at interface
            if (interact_mode == InteractMode::REFLECT)
                after_dir = reflect(b4_dir, normal);
            else if (interact_mode == InteractMode::TRANSMIT)
                after_dir = bend(b4_dir, normal, *before->Indx, *after.Indx);
            else if (interact_mode == InteractMode::DUMMY)
                after_dir = b4_dir;
            else if (interact_mode == InteractMode::PHANTOM)
                after_dir = b4_dir;
            else // no action, input becomes output
                after_dir = b4_dir;

            before_pt = inc_pt;
            before_normal = normal;
            before_dir = after_dir;
            z_dir_before = z_dir_after;
            b4_interact_mode = interact_mode;
            before = &after;
            tfrm_from_before = *before->Tfrm;

        } catch (TraceMissedSurfaceException &ray_miss) {
            // Annotate the in-flight exception with where the ray got to and
            // rethrow it: callers pull ray_pkg back out and carry on. The `ifc`
            // field the Java also sets here is dropped -- nothing reads it.
            ray.push_back(RaySeg(before_pt, before_dir, pp_dst, before_normal));
            ray_miss.surf = surf;
            ray_miss.prev_tfrm = *before->Tfrm;
            ray_miss.ray_pkg = std::make_shared<const RayPkg>(ray, opl, wvl);
            throw;
        } catch (TraceTIRException &ray_tir) {
            ray.push_back(RaySeg(inc_pt, before_dir, 0.0, normal));
            ray_tir.surf = surf;
            ray_tir.int_pt = inc_pt;
            ray_tir.ray_pkg = std::make_shared<const RayPkg>(ray, opl, wvl);
            throw;
        } catch (TraceRayBlockedException &ray_blocked) {
            ray.push_back(RaySeg(inc_pt, before_dir, 0.0, normal));
            ray_blocked.surf = surf;
            ray_blocked.ray_pkg = std::make_shared<const RayPkg>(ray, opl, wvl);
            throw;
        }
    }
    return std::make_shared<const RayPkg>(ray, op_delta, wvl);
}

} // namespace redukti::rayoptics::raytr
