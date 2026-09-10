// C++ port of org.redukti.rayoptics.raytr.Wideangle
#include "redukti/rayoptics/raytr/Wideangle.h"

#include "redukti/Exceptions.h"
#include "redukti/Text.h"
#include "redukti/mathlib/BrentSolver.h"
#include "redukti/mathlib/M.h"
#include "redukti/mathlib/Matrix3.h"
#include "redukti/mathlib/SecantSolver.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/parax/ParaxTypes.h"
#include "redukti/rayoptics/raytr/RayTrace.h"
#include "redukti/rayoptics/raytr/Trace.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/rayoptics/util/Lists.h"

#include <cmath>
#include <cstdio>

namespace redukti::rayoptics::raytr {

namespace M = mathlib::M;
using exceptions::TraceException;
using exceptions::TraceMissedSurfaceException;
using mathlib::Matrix3;
using mathlib::Vector2;
using mathlib::Vector3;
using util::Lists::get;

std::string Wideangle::ZEnpStopHt::toString() const {
    return "{z_enp=" + doubleToString(z_enp) +
           ", ht_at_stop=" + doubleToString(ht_at_stop) + "}";
}

/**
 * Trace a ray thru the center of the entrance pupil at z_enp.
 *
 *     Args:
 *
 *         z_enp:      The z distance from the 1st interface of the
 *                     entrance pupil for the field angle dir0
 *         seq_model:  The sequential model
 *         stop_idx:   index of the aperture stop interface
 *         dir0:       direction cosine vector in object space
 *         obj_dist:   object distance to first interface
 *         wvl:        wavelength of raytrace (nm)
 *
 * Returns the intercept coord on the stop surface on success
 * Along with the ray tracing results
 */
RayResultWithStopCoord Wideangle::enp_z_coordinate(double z_enp,
                                                   seq::SequentialModel *seq_model,
                                                   int stop_idx, const Vector3 &dir0,
                                                   double obj_dist, double wvl) {
    auto obj2enp_dist = (obj_dist + z_enp);
    Vector3 pt1(0., 0., obj2enp_dist);
    auto rot_mat = Matrix3::rot_v1_into_v2(Vector3::vector3_001, dir0);
    auto pt0 = rot_mat.multiply(pt1.negate()).plus(pt1);
    std::shared_ptr<const RayPkg> ray_pkg;
    RayResult rr;
    Vector3 final_coord = Vector3::ZERO;
    try {
        RayTraceOptions options;
        options.intersect_obj = false;
        ray_pkg = RayTrace::trace(seq_model, pt0, dir0, wvl, options);
        rr = RayResult(ray_pkg, nullptr);
        final_coord = get(ray_pkg->ray, stop_idx).p;
    } catch (TraceException &ray_error) {
        //print(f'  ray_error: "{type(ray_error).__name__}", '
        //              f'{ray_error.surf=}')
        //logger.debug(f'   ray_error: "{type(ray_error).__name__}", '
        //f'{ray_error.surf=}')
        ray_pkg = ray_error.ray_pkg;
        rr = RayResult(ray_pkg, ray_error.clone());
        // FIXME should below be null?
        final_coord = Vector3::ZERO;
    }
    return RayResultWithStopCoord(final_coord, rr, stop_idx);
}

std::optional<double> Wideangle::Enp_z_coordinate_wrapper::eval(double z_enp) {
    auto coord_rr = enp_z_coordinate(z_enp, seq_model, stop_idx, dir0, obj_dist, wvl);
    auto final_coord = coord_rr.stop_coord;
    auto rr = coord_rr.rr;
    if (rr.err == nullptr) {
        auto ht_at_stop = final_coord.y;
        return ht_at_stop;
    } else
        return std::nullopt;
}

std::optional<double> Wideangle::Eval_Z_Enp_Function::eval(double z_enp) {
    auto coord_rr = enp_z_coordinate(z_enp, seq_model, stop_idx, dir0, obj_dist, wvl);
    auto final_coord = coord_rr.stop_coord;
    rr = coord_rr.rr;
    rr_set = true;
    return final_coord.y - y_target;
}

RayResultWithZEnp Wideangle::find_real_enp(optical::OpticalModel *opm,
                                           std::optional<int> stop_idx,
                                           specs::Field &fld, double wvl,
                                           const std::string &selector) {
    if (selector == "rev1")
        return find_real_enp_rev1(opm, stop_idx, fld, wvl, std::nullopt);
    else
        return find_real_enp_orig(opm, stop_idx, fld, wvl);
}

RayResultWithZEnp Wideangle::find_real_enp(optical::OpticalModel *opm,
                                           std::optional<int> stop_idx,
                                           specs::Field &fld, double wvl) {
    return find_real_enp(opm, stop_idx, fld, wvl, "rev1");
}

Wideangle::ZEnpStopHt Wideangle::find_edge(mathlib::ScalarObjectiveFunction &f, double a,
                                           double b, std::optional<int> max_iter_) {
    // use binary search to find the edge of the fct's range.
    int max_iter = max_iter_.has_value() ? *max_iter_ : 3;
    auto fa = f.eval(a);
    auto fb = f.eval(b);
    for (int i = 0; i < max_iter; i++) {
        auto c = a + (b - a) / 2;
        auto fc = f.eval(c);
        // A null result means the ray failed at c, so c becomes the far end of
        // the bracket; otherwise it becomes the near end.
        if (!fc.has_value()) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }
    // Java unboxes fa/fb here and would NPE on a null; .value() throws in the
    // same situation.
    if (!fb.has_value())
        return ZEnpStopHt(a, fa.value());
    else
        return ZEnpStopHt(b, fb.value());
}

RayResultWithZEnp Wideangle::find_real_enp_rev1(optical::OpticalModel *opm,
                                                std::optional<int> stop_idx_,
                                                specs::Field &fld, double wvl,
                                                std::optional<bool> check_direction_) {
    bool check_direction = check_direction_.has_value() ? *check_direction_ : true;
    auto sm = opm->seq_model.get();
    auto osp = opm->optical_spec.get();
    auto &fod = osp->parax_data->fod;
    int stop_idx = stop_idx_.has_value() ? *stop_idx_ : 1;
    auto coord = osp->obj_coords(fld);
    auto dir0 = coord.dir;
    // If there is aim_info, try it and return if good.
    if (fld.z_enp.has_value()) {
        auto z_enp_f = *fld.z_enp;
        auto coord_rr =
            enp_z_coordinate(z_enp_f, sm, stop_idx, dir0, fod.obj_dist, wvl);
        auto final_coord = coord_rr.stop_coord;
        auto rr_f = coord_rr.rr;
        auto tol = 1.48e-08;
        if (std::abs(final_coord.y) < tol)
            return RayResultWithZEnp(z_enp_f, rr_f);
    }
    // filter on-axis chief ray. z_enp is the paraxial result.
    auto z_enp_0 = fod.enp_dist;
    if (dir0.z == 1.) {
        // axial chief ray
        auto coord_rr = enp_z_coordinate(z_enp_0, sm, stop_idx, dir0, fod.obj_dist, wvl);
        auto rr_f = coord_rr.rr;
        //logger.info(f"  axial chief {z_enp_0=:8.4f}  {rr.err is None}")
        return RayResultWithZEnp(z_enp_0, rr_f);
    }
    std::optional<ZEnpStopHt> start_z;
    std::optional<ZEnpStopHt> prev_z;
    std::optional<ZEnpStopHt> end_z;
    auto del_z = -z_enp_0 / 16.0;  // step z
    auto z_enp = z_enp_0;
    bool keep_going = true;
    std::string direction = "first";
    int first_surf_misses = 0;
    // protect against infinite loops
    int trial = 0;
    // if the trace succeeds 5 times in a row, go on to the next phase
    int successes = 0;
    // last ray result seen, reported to the caller when no pupil is found
    RayResult rr;
    while (keep_going && trial < 64 && first_surf_misses < 2) {
        auto coord_rr = enp_z_coordinate(z_enp, sm, stop_idx, dir0, fod.obj_dist, wvl);
        rr = coord_rr.rr;
        auto final_coord = coord_rr.stop_coord;
        if (rr.err == nullptr) {
            auto ht_at_stop = final_coord.y;
            //            logger.debug(f"  ray passed at z_enp={z_enp:10.5f},  "
            //                         f"{ht_at_stop=:7.3f}")
            successes++;
            if (!start_z.has_value())
                start_z = ZEnpStopHt(z_enp, ht_at_stop);
            prev_z = end_z;
            end_z = ZEnpStopHt(z_enp, ht_at_stop);
            if (successes > 1) {
                // check for a zero crossing, if so, we're done
                if (prev_z->ht_at_stop * end_z->ht_at_stop < 0)
                    keep_going = false;
            }
            if (successes == 2 && check_direction) {
                // check that we're searching in the right direction
                if (std::abs(start_z->ht_at_stop) < std::abs(end_z->ht_at_stop)) {
                    if (direction == "first") {
                        // first time through, reverse direction and start on
                        // the other side of z_enp_0.
                        //logger.debug("  --> reverse search direction")
                        del_z = -del_z;
                        z_enp = z_enp_0;
                        direction = "reverse";
                        auto tmp = end_z;
                        end_z = start_z;
                        start_z = tmp;
                    }
                }
            }
        } else {
            // logger.debug(f"  ray failed at z_enp={z_enp:10.5f}, "
            //                         f"{type(rr.err).__name__} at surf {rr.err.surf}")
            if (dynamic_cast<TraceMissedSurfaceException *>(rr.err.get()) != nullptr) {
                // if the first surface was missed, then exit
                if (rr.err->surf == 1) {
                    //logger.debug(f"Num 1st surf misses {first_surf_misses}: "
                    //                                 +msg1)
                    del_z = -del_z;
                    z_enp = z_enp_0;
                    first_surf_misses++;
                }
            }
            // if the first surface was missed, then exit
            if (start_z.has_value()) {
                if (direction == "first") {
                    // logger.debug("  --> reverse search direction")
                    del_z = -del_z;
                    z_enp = z_enp_0;
                    direction = "reverse";
                    auto tmp = end_z;
                    end_z = start_z;
                    start_z = tmp;
                } else
                    keep_going = false;
            }
        }
        z_enp += del_z;
        if (M::isZero(z_enp)) {
            // if we sample directly at z_enp = 0, i.e. the vertex of the first
            // surface, the surface intersection method won't find the right
            // root.Perturb the sample point slightly away from zero to avoid this problem
            z_enp = del_z / 10.0;
        }
        trial += 1;
    }
    if (!start_z.has_value() || !end_z.has_value())
        return RayResultWithZEnp(std::nullopt, rr);
    auto z_enp_a = start_z->z_enp;
    auto ht_at_stop_a = start_z->ht_at_stop;
    auto z_enp_b = end_z->z_enp;
    auto ht_at_stop_b = end_z->ht_at_stop;
    double a = 0.0, b = 0.0;
    // If start and end are equal, then only one ray was successful.
    // Sample z_enp evenly 1 del_z to either side.
    if (z_enp_a == z_enp_b) {
        auto start_new = z_enp_a - del_z;
        auto end_new = z_enp_b + del_z;
        start_z = std::nullopt;
        end_z = std::nullopt;
        for (auto x : linspace(start_new, end_new, 8)) {
            z_enp = x;
            auto coord_rr =
                enp_z_coordinate(z_enp, sm, stop_idx, dir0, fod.obj_dist, wvl);
            auto final_coord = coord_rr.stop_coord;
            rr = coord_rr.rr;
            if (rr.err == nullptr) {
                auto ht_at_stop = final_coord.y;
                if (!start_z.has_value())
                    start_z = ZEnpStopHt(z_enp, ht_at_stop);
                end_z = ZEnpStopHt(z_enp, ht_at_stop);
            }
        }
        if (!start_z.has_value() || !end_z.has_value())
            return RayResultWithZEnp(std::nullopt, rr);
        a = start_z->z_enp;
        b = end_z->z_enp;
    // test for crossing between the end points
    } else if (ht_at_stop_a * ht_at_stop_b < 0) {
        // yes, there was a crossing somewhere in the interval
        //        # set the bracket using start_z
        a = z_enp_a;
        b = z_enp_b;
        if (prev_z.has_value()) {
            // if there's a previous sample, see if the crossing can be
            //            # more tightly bracketed.
            auto z_enp_c = prev_z->z_enp;
            auto ht_at_stop_c = prev_z->ht_at_stop;
            if (ht_at_stop_c * ht_at_stop_b < 0) {
                // set the smallest bracket using prev_z
                start_z = prev_z;
                a = z_enp_c;
                b = z_enp_b;
            }
        }
    } else {
        // we haven't found a zero crossing yet
        //        # refine the interval by finding the effective "edge" of the beam
        Enp_z_coordinate_wrapper z_enp_coordinate_wrapper(sm, stop_idx, dir0,
                                                          fod.obj_dist, wvl);
        auto edge_b = find_edge(z_enp_coordinate_wrapper, z_enp_b, z_enp_b + del_z, 6);
        auto z_enp_edge_b = edge_b.z_enp;
        auto ht_at_stop_edg_b = edge_b.ht_at_stop;
        //logger.debug(f"  edge_b found at at z_enp={z_enp_edge_b:10.5f},  "
        //                     f"{ht_at_stop_edg_b=:7.3f}")
        if (ht_at_stop_edg_b * ht_at_stop_b < 0) {
            start_z = ZEnpStopHt(z_enp_b, ht_at_stop_b);
            end_z = ZEnpStopHt(z_enp_edge_b, ht_at_stop_edg_b);
            a = z_enp_b;
            b = z_enp_edge_b;
        } else {
            // find the other effective "edge" of the beam
            Enp_z_coordinate_wrapper wrapper2(sm, stop_idx, dir0, fod.obj_dist, wvl);
            auto edge_a = find_edge(wrapper2, z_enp_a, z_enp_a - del_z, 6);
            auto z_enp_edge_a = edge_a.z_enp;
            auto ht_at_stop_edg_a = edge_a.ht_at_stop;
            // logger.debug(f"  edge_a found at at z_enp={z_enp_edge_a:10.5f},  "
            //                         f"{ht_at_stop_edg_a=:7.3f}")
            if (ht_at_stop_edg_a * ht_at_stop_a < 0) {
                // found an interval containing a crossover point
                start_z = ZEnpStopHt(z_enp_a, ht_at_stop_a);
                end_z = ZEnpStopHt(z_enp_edge_a, ht_at_stop_edg_a);
                a = z_enp_a;
                b = z_enp_edge_a;
            } else {
                //  there is no ray that passes thru the center of the stop
                //                # surface.
                std::fprintf(stderr, "chief ray trace failed at field %3.1f\n",
                             fld.yv());
                auto z_enp_cntr = z_enp_edge_a + (z_enp_edge_b - z_enp_edge_a) / 2;
                auto coord_rr =
                    enp_z_coordinate(z_enp_cntr, sm, stop_idx, dir0, fod.obj_dist, wvl);
                rr = coord_rr.rr;
                // logger.debug(f"  fld: {fld.yv:3.1f}:   {z_enp_edge_a=:8.4f}  "
                //                    f"{z_enp_edge_b=:8.4f}  {z_enp_cntr=:8.4f}  "
                //                    f"{ht_at_stop=:10.2e}")
                return RayResultWithZEnp(z_enp_b, rr);
            }
        }
    }
    // compute the straightline crossing pt given the interval
    double z_estimate;
    if (M::is_fuzzy_zero(end_z->ht_at_stop - start_z->ht_at_stop)) {
        z_estimate = start_z->z_enp;
    } else {
        z_estimate = start_z->z_enp - ((end_z->z_enp - start_z->z_enp) /
                                       (end_z->ht_at_stop - start_z->ht_at_stop)) *
                                          start_z->ht_at_stop;
    }
    //    logger.debug(f"  trials: {trial},   {successes=}")
    //    logger.debug(f"  z_enp: start_z={a:10.5f} z_estimate={z_estimate:10.5f}  "
    //                 f"end_z={b:10.5f}")
    //    logger.debug(f"  ht_at_stop: start_z={start_z[1]:10.5f} "
    //                 f"end_z={end_z[1]:10.5f}")
    // Vector3 has no default constructor, so seed the pair explicitly.
    util::Pair<Vector3, RayResult> result(Vector3::ZERO, RayResult());
    try {
        result = find_z_enp_on_interval(opm, stop_idx, a, b, z_estimate, fld, wvl);
    } catch (const IndexOutOfBoundsException &) {
        return RayResultWithZEnp(std::nullopt, rr);
    } catch (const IllegalArgumentException &) {
        // upstream: except (IndexError, ValueError): return None, rr
        return RayResultWithZEnp(std::nullopt, rr);
    }
    auto start_coord = result.first;
    rr = result.second;
    z_enp = start_coord.z;
    return RayResultWithZEnp(z_enp, rr);
}

util::Pair<Vector3, RayResult> Wideangle::find_z_enp_on_interval(
    optical::OpticalModel *opt_model, std::optional<int> stop_idx, double start_z,
    double end_z, double z_estimate, specs::Field &fld, double wvl) {
    RayResult rr;
    auto sm = opt_model->seq_model.get();
    auto osp = opt_model->optical_spec.get();
    auto &fod = osp->parax_data->fod;
    auto z_enp = z_estimate;
    auto obj_dist = fod.obj_dist;
    auto coord = osp->obj_coords(fld);
    auto dir0 = coord.dir;
    double y_target = 0.;
    Vector3 start_coords = Vector3::ZERO;
    bool converged = false;
    if (stop_idx.has_value()) {
        // do 1D iteration if field and target points are zero in x
        Eval_Z_Enp_Function fn(sm, *stop_idx, dir0, obj_dist, wvl, y_target);
        try {
            auto result = mathlib::SecantSolver::find_root(fn, z_enp, 50, 1.48e-8);
            z_enp = result.root;
            converged = result.converged;
            rr = fn.rr;
            auto ht_at_stop = get(rr.pkg->ray, *stop_idx).p.y;
            if (std::abs(ht_at_stop - y_target) < 1e-6)
                converged = true;
            start_coords = Vector3(0., 0., z_enp);
        } catch (const TraceException &) {
            // the objective records the ray result on every evaluation;
            // hold on to the last one so the caller still has ray data.
            rr = fn.rr;
            converged = false;
            start_coords = Vector3(0., 0., z_enp);
        }
        if (!converged) {
            //                logger.debug(f'  {results.method} converged: '
            //                             f'{results.converged},  # fct evals='
            //                             f'{results.function_calls}  msg: "{results.flag}" '
            //                             f'{z_enp=:9.4f}')
            try {
                auto result = mathlib::BrentSolver::find_root(start_z, end_z, fn);
                if (result.converged) {
                    z_enp = result.root;
                    converged = true;
                    start_coords = Vector3(0., 0., z_enp);
                    // else keep the estimate from the secant iteration, as
                    // upstream does when brentq fails to converge. A bad
                    // bracket makes BrentSolver report root=0, which must not
                    // be taken as an answer.
                }
            } catch (const TraceException &) {
            // ditto - retain the secant estimate rather than failing
            }
            if (fn.rr_set)
                rr = fn.rr;
        }
    } else
        start_coords = Vector3(0., 0., fod.enp_dist);
    //    logger.debug(f'  {results.method} converged: {results.converged},  '
    //                 f'# fct evals={results.function_calls}  msg: "{results.flag}"')
    return util::Pair<Vector3, RayResult>(start_coords, rr);
}

std::vector<double> Wideangle::linspace(double start, double end, int num) {
    std::vector<double> result(static_cast<std::size_t>(num), 0.0);
    if (num == 1) {
        result[0] = start;
        return result;
    }
    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; i++) {
        result[static_cast<std::size_t>(i)] = start + step * i;
    }
    return result;
}

/**
 * Locate the z center of the real pupil for `fld`, wrt 1st ifc
 *
 *     This function implements a 2 step process to finding the chief ray
 *     for `fld` and `wvl` for wide angle systems. `fld` should be of type
 *     ('object', 'angle'), even for finite object distances.
 *
 *     The first phase searches for the window of pupil locations by sampling the
 *     z coordinate from the paraxial pupil location towards the first interface
 *     vertex. Failed rays are discarded until a range of z coordinates is found
 *     where rays trace successfully. If the search forward is unsuccessful (i.e.
 *     winds up missing the 1st surface), the search is restarted moving away from
 *     the first interface. If only a single successful trace is in hand, a
 *     second, more finely subdivided search is conducted about the successful
 *     point.
 *
 *     The outcome is a range, start_z -> end_z, that is divided in 3 and a ray
 *     iteration (using :func:`~.raytr.wideangle.find_z_enp`) to find the center
 *     of the stop surface is done. Sometimes the start point doesn't produce a
 *     solution; use of the mid-point as a start is a reliable second try.
 */
RayResultWithZEnp Wideangle::find_real_enp_orig(optical::OpticalModel *opm,
                                                std::optional<int> stop_idx_,
                                                specs::Field &fld, double wvl) {
    auto sm = opm->seq_model.get();
    auto osp = opm->optical_spec.get();
    auto &fod = osp->parax_data->fod;
    int stop_idx = stop_idx_.has_value() ? *stop_idx_ : 1;
    RayResult rr;
    auto coord = osp->obj_coords(fld);
    auto dir0 = coord.dir;
    if (fld.z_enp.has_value()) {
        auto z_enp_f = *fld.z_enp;
        auto coord_rr =
            enp_z_coordinate(z_enp_f, sm, stop_idx, dir0, fod.obj_dist, wvl);
        double tol = 1.48e-8;
        if (std::abs(coord_rr.stop_coord.y) < tol)
            return RayResultWithZEnp(z_enp_f, coord_rr.rr);
    }
    auto z_enp_0 = fod.enp_dist;
    if (dir0.z == 1.0) {  // axial chief ray
        auto coord_rr = enp_z_coordinate(z_enp_0, sm, stop_idx, dir0, fod.obj_dist, wvl);
        return RayResultWithZEnp(z_enp_0, coord_rr.rr);
    }
    std::optional<double> start_z;
    std::optional<double> end_z;
    auto del_z = -z_enp_0 / 16.0;
    auto z_enp = z_enp_0;
    bool keep_going = true;
    int first_surf_misses = 0;
    // protect against infinite loops
    int trial = 0;
    // if the trace succeeds 5 times in a row, go on to the next phase
    int successes = 0;
    while (keep_going && successes < 4 && trial < 64 && first_surf_misses < 2) {
        auto coord_rr = enp_z_coordinate(z_enp, sm, stop_idx, dir0, fod.obj_dist, wvl);
        rr = coord_rr.rr;
        if (rr.err == nullptr) {
            successes++;
            if (!start_z.has_value())
                start_z = z_enp;
            end_z = z_enp;
        } else {
            if (dynamic_cast<TraceMissedSurfaceException *>(rr.err.get()) != nullptr) {
                if (rr.err->surf == 1) {
                    del_z = -del_z;
                    z_enp = z_enp_0;
                    first_surf_misses++;
                }
                // if the first surface was missed, then exit
                if (start_z.has_value())
                    keep_going = false;
            }
        }
        z_enp += del_z;
        trial += 1;
    }
    // If start and end are equal, then only one ray was successful.
    // Sample z_enp evenly 1 del_z to either side.
    if (start_z.has_value() && end_z.has_value() && *start_z == *end_z) {
        auto start_new = *start_z - del_z;
        auto end_new = *end_z + del_z;
        start_z = std::nullopt;
        end_z = std::nullopt;
        for (auto x : linspace(start_new, end_new, 8)) {
            z_enp = x;
            auto coord_rr =
                enp_z_coordinate(z_enp, sm, stop_idx, dir0, fod.obj_dist, wvl);
            rr = coord_rr.rr;
            if (rr.err == nullptr) {
                if (!start_z.has_value())
                    start_z = z_enp;
                end_z = z_enp;
            }
        }
    }
    // Java unboxes start_z/end_z here; a null would NPE, and .value() throws in
    // the same situation.
    std::vector<double> starting_pts = {start_z.value(),
                                        (start_z.value() + end_z.value()) / 2.0,
                                        // Now that candidate z_enps have been identified that trace without
                                        // ray failures, iterate to find the ray thru the stop center
                                        end_z.value()};
    for (auto init_z : starting_pts) {
        auto result = find_z_enp(opm, stop_idx, init_z, fld, wvl);
        rr = result.rr;
        z_enp = result.z_enp.has_value() ? *result.z_enp : z_enp;
        if (rr.err == nullptr)
            break;
    }
    return RayResultWithZEnp(z_enp, rr);
}

/**
 * iterates a ray to [0, 0] on interface stop_ifc, returning aim info
 *
 *     Args:
 *
 *         opt_model:  input OpticalModel
 *         stop_idx:   index of the aperture stop interface
 *         z_enp_0:    estimate of pupil location. this estimate must support
 *                     a raytrace up to stop_ifc
 *         fld:        field point
 *         wvl:        wavelength of raytrace (nm)
 *
 *     Returns z distance from 1st interface to the entrance pupil.
 *
 *     If stop_ifc is None, i.e. a floating stop surface, returns paraxial
 *     entrance pupil.
 *
 *     If the iteration fails, a TraceError will be raised
 */
RayResultWithZEnp Wideangle::find_z_enp(optical::OpticalModel *opt_model,
                                        std::optional<int> stop_idx, double z_enp_0,
                                        specs::Field &fld, double wvl) {
    RayResult rr;
    auto seq_model = opt_model->seq_model.get();
    auto osp = opt_model->optical_spec.get();
    auto &fod = osp->parax_data->fod;
    auto z_enp = z_enp_0;
    auto obj_dist = fod.obj_dist;
    auto coord = osp->obj_coords(fld);
    auto dir0 = coord.dir;
    double y_target = 0.;
    if (stop_idx.has_value()) {
        // do 1D iteration if field and target points are zero in x
        Eval_Z_Enp_Function func(seq_model, *stop_idx, dir0, obj_dist, wvl, y_target);
        try {
            auto result = mathlib::SecantSolver::find_root(func, z_enp, 50, 1.48e-8);
            z_enp = result.root;
        } catch (const TraceException &) {
            z_enp = z_enp_0;
            // Preserve the caller's estimate rather than replacing it with zero.
        }
        rr = func.rr;
    } else
        z_enp = fod.enp_dist;
    return RayResultWithZEnp(z_enp, rr);
}

/**
 * Trace reverse ray from image point to get object space inputs.
 *
 *     This function traces the chief ray for `fld` and `wvl` through the center of the stop surface, starting from the specified real image height.
 *
 *     This is the implementation of :meth:`~.raytr.opticalspec.FieldSpec.obj_coords` for ('image', 'real height')
 */
RayDataWithZ_Enp Wideangle::eval_real_image_ht(optical::OpticalModel *opt_model,
                                               specs::Field &fld, double wvl) {
    auto sm = opt_model->seq_model.get();
    auto osp = opt_model->optical_spec.get();
    auto fov = osp->fov.get();
    auto &fod = osp->parax_data->fod;
    auto not_wa = !fov->is_wide_angle;
    auto stop_idx = !sm->stop_surface.has_value() ? 1 : *sm->stop_surface;
    auto ifcx = static_cast<int>(sm->ifcs.size()) - stop_idx - 1;
    auto rpath = sm->reverse_path(wvl, static_cast<int>(sm->ifcs.size()), std::nullopt,
                                  -1);
    auto eprad = fod.exp_radius;
    auto obj2pup_dist = fod.exp_dist - fod.img_dist;
    Vector3 p_exp(0., 0., obj2pup_dist);
    auto xy_target = Vector2::vector2_0;
    Vector3 p_i(fld.x, fld.y, 0);
    if (fov->is_relative)
        p_i = p_i.times(fov->value);
    auto d_i = p_exp.minus(p_i).normalize();
    auto arr = xy_target.as_array();
    auto result = Trace::iterate_ray_raw(rpath, ifcx, std::vector<double>{arr[0], arr[1]},
                                         p_i, d_i, obj2pup_dist, eprad, wvl, not_wa);
    auto rrev_cr = result.rr;
    auto p_k = get(rrev_cr.pkg->ray, -2).p;
    auto p_k01 = std::sqrt(p_k.x * p_k.x + p_k.y * p_k.y);
    auto d_k = get(rrev_cr.pkg->ray, -2).d;
    auto d_o = d_k.negate();
    auto d_k01 = std::sqrt(d_k.x * d_k.x + d_k.y * d_k.y);
    double z_enp;
    if (d_k01 == 0.)
        z_enp = fod.enp_dist;
    else
        z_enp = p_k.z + p_k01 * d_o.z / d_k01;
    auto p_o = get(rrev_cr.pkg->ray, -1).p;
    if (osp->conjugate_type(specs::ImageKey::Object) == specs::ConjugateType::INFINITE) {
        auto obj2enp_dist = fod.obj_dist + z_enp;
        Vector3 enp_pt(0, 0, obj2enp_dist);
        p_o = enp_pt.plus(d_k.times(obj2enp_dist));
    }
    return RayDataWithZ_Enp(RayData(p_o, d_o), z_enp);
}

} // namespace redukti::rayoptics::raytr
