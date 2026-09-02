// C++ port of org.redukti.rayoptics.raytr.{WaveAbr,FinitePupilWaveAberrationResult}
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_RAYTR_WAVEABR_H
#define REDUKTI_RAYOPTICS_RAYTR_WAVEABR_H

#include "redukti/rayoptics/parax/ParaxTypes.h"
#include "redukti/rayoptics/raytr/RayTypes.h"
#include "redukti/rayoptics/util/Tuples.h"

#include <memory>
#include <optional>

namespace redukti::rayoptics::optical {
class OpticalModel;
}

namespace redukti::rayoptics::raytr {

/** Java's `record FinitePupilWaveAberrationResult(...)`. */
class FinitePupilWaveAberrationResult {
public:
    std::shared_ptr<const RayPkg> ray_pkg;
    std::shared_ptr<const ChiefRayPkg> chief_ray_pkg;
    std::shared_ptr<const ReferenceSphere> ref_sphere;
    double e1;
    double ekp;
    double ep;
    /** Null when the discriminant is negative. */
    std::optional<mathlib::Vector3> ray_exit_pupil_coord;
    double ray_op;
    double cr_op;

    FinitePupilWaveAberrationResult(
        std::shared_ptr<const RayPkg> ray_pkg_,
        std::shared_ptr<const ChiefRayPkg> chief_ray_pkg_,
        std::shared_ptr<const ReferenceSphere> ref_sphere_, double e1_, double ekp_,
        double ep_, std::optional<mathlib::Vector3> ray_exit_pupil_coord_,
        double ray_op_, double cr_op_)
        : ray_pkg(std::move(ray_pkg_)), chief_ray_pkg(std::move(chief_ray_pkg_)),
          ref_sphere(std::move(ref_sphere_)), e1(e1_), ekp(ekp_), ep(ep_),
          ray_exit_pupil_coord(ray_exit_pupil_coord_), ray_op(ray_op_), cr_op(cr_op_) {}
};

class WaveAbr {
public:
    static std::shared_ptr<const ReferenceSphere> calculate_reference_sphere(
        optical::OpticalModel *opt_model, specs::Field &fld, double wvl, double foc,
        const ChiefRayPkg &chief_ray_pkg, std::optional<mathlib::Vector2> image_pt_2d,
        std::optional<mathlib::Vector2> image_delta);

    static std::shared_ptr<const ChiefRayExitPupilSegment> transfer_to_exit_pupil(
        std::shared_ptr<seq::Interface> ifc, const RayData &ray_seg,
        double exp_dst_parax);

    static double eic_distance(const RayData &r, const RayData &r0);

    static double ray_dist_to_perp_from_origin(const RayData &r);

    /** Closest points on two skew rays, each with its distance along the ray. */
    static util::Pair<util::Pair<mathlib::Vector3, double>,
                      util::Pair<mathlib::Vector3, double>>
    dist_to_shortest_join(const RayData &r1, const RayData &r2);

    static double wave_abr_full_calc(const parax::FirstOrderData &fod, specs::Field &fld,
                                     double wvl, double foc,
                                     const std::shared_ptr<const RayPkg> &ray_pkg,
                                     const std::shared_ptr<const ChiefRayPkg> &chief_ray_pkg,
                                     const std::shared_ptr<const ReferenceSphere> &ref_sphere);

    static FinitePupilWaveAberrationResult wave_abr_calc_finite_pupil(
        const std::shared_ptr<const RayPkg> &ray_pkg,
        const std::shared_ptr<const ChiefRayPkg> &chief_ray_pkg,
        const std::shared_ptr<const ReferenceSphere> &ref_sphere);

private:
    static double wave_abr_full_calc_finite_pup(
        const parax::FirstOrderData &fod, specs::Field &fld, double wvl, double foc,
        const std::shared_ptr<const RayPkg> &ray_pkg,
        const std::shared_ptr<const ChiefRayPkg> &chief_ray_pkg,
        const std::shared_ptr<const ReferenceSphere> &ref_sphere);

    static double wave_abr_full_calc_inf_ref(
        const parax::FirstOrderData &fod, specs::Field &fld, double wvl, double foc,
        const std::shared_ptr<const RayPkg> &ray_pkg,
        const std::shared_ptr<const ChiefRayPkg> &chief_ray_pkg,
        const std::shared_ptr<const ReferenceSphere> &ref_sphere);
};

} // namespace redukti::rayoptics::raytr

#endif // REDUKTI_RAYOPTICS_RAYTR_WAVEABR_H
