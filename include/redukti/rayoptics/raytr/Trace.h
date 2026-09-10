// C++ port of org.redukti.rayoptics.raytr.Trace
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_RAYTR_TRACE_H
#define REDUKTI_RAYOPTICS_RAYTR_TRACE_H

#include "redukti/rayoptics/math/Tfm3d.h"
#include "redukti/rayoptics/raytr/RayTrace.h"
#include "redukti/rayoptics/raytr/RayTypes.h"
#include "redukti/rayoptics/seq/SurfaceData.h"

#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace redukti::rayoptics::optical {
class OpticalModel;
}
namespace redukti::rayoptics::seq {
class SequentialModel;
}

namespace redukti::rayoptics::raytr {

/** Java's `record GaussianQuadraturePoint(Vector2 pupil, double weight)`. */
class GaussianQuadraturePoint {
public:
    mathlib::Vector2 pupil;
    double weight;

    GaussianQuadraturePoint(const mathlib::Vector2 &pupil_, double weight_)
        : pupil(pupil_), weight(weight_) {}
};

// r = pd.DataFrame(ray, columns=['inc_pt', 'after_dir',
//                                   'after_dst', 'normal'])
/** Java's `record RayDataFrame(...)` neighbours, kept with their owner. */
class RayDataFrame {
public:
    std::vector<mathlib::Vector3> inc_pt;
    std::vector<mathlib::Vector3> after_dir;
    std::vector<double> after_dst;
    std::vector<mathlib::Vector3> normal;

    explicit RayDataFrame(const std::vector<RaySeg> &raySegList);
};

class RayDataFrameByField {
public:
    specs::Field *fld;
    std::vector<RayDataFrame> frames;

    RayDataFrameByField(specs::Field *fld_, std::vector<RayDataFrame> frames_)
        : fld(fld_), frames(std::move(frames_)) {}
};

class Trace {
public:
    /**
     * Trace a single ray via pupil, field and wavelength specs.
     *
     *     This function traces a single ray at a given wavelength, pupil and field specification.
     *
     *     Ray failures (miss surface, TIR) and aperture clipping are handled via RayError exceptions. If a failure occurs, a second item is returned (if  *rayerr_filter* is set to 'summary' or 'full') that contains information about the failure. Apertures are tested using the :meth:`~.seq.interface.Interface.point_inside` API when *check_apertures* is True.
     *
     *     The pupil coordinates by default are normalized to the vignetted pupil extent. Alternatively, the pupil coordinates can be taken as actual coordinates on the pupil plane (and similarly for ray direction) using the **pupil_type** keyword argument.
     *
     *     The amount of output that is returned can range from the entire ray (default) to the image segment only or even the return from a user-supplied filtering function.
     *
     *     Args:
     *         opt_model: :class:`~.OpticalModel` instance
     *         pupil: 2d vector of relative pupil coordinates
     *         fld: :class:`~.Field` point for wave aberration calculation
     *         wvl: wavelength of ray (nm)
     *
     *         check_apertures: if True, do point_inside() test on inc_pt
     *         apply_vignetting: if True, apply the `fld` vignetting factors to **pupil**
     *
     *         pupil_type: ::
     *
     *             - 'rel pupil': relative pupil coordinates
     *             - 'aim pt': aim point on pupil plane
     *             - 'aim dir': aim direction in object space
     *
     *         use_named_tuples: if True, returns data as RayPkg and RaySeg.
     *
     *         output_filter: ::
     *
     *             - if None, append entire ray
     *             - if 'last', append the last ray segment only
     *             - else treat as callable and append the return value
     *
     *         rayerr_filter: ::
     *
     *             - if None, on ray error append nothing
     *             - if 'summary', append the exception without ray data
     *             - if 'full', append the exception with ray data up to error
     *             - else append nothing
     *
     *         eps: accuracy tolerance for surface intersection calculation
     *
     *     Returns:
     *         tuple: ray_pkg, trace_error | None
     */
    static RayResult trace_ray(optical::OpticalModel *opt_model,
                               const mathlib::Vector2 &pupil, specs::Field &fld,
                               double wvl, TraceOptions &trace_options);

    /**
     * Wrapper for trace_base that handles exceptions.
     *
     *     Args:
     *         opt_model: :class:`~.OpticalModel` instance
     *         pupil: 2d vector of relative pupil coordinates
     *         fld: :class:`~.Field` point for wave aberration calculation
     *         wvl: wavelength of ray (nm)
     *         output_filter: ::
     *
     *             - if None, append entire ray
     *             - if 'last', append the last ray segment only
     *             - else treat as callable and append the return value
     *
     *         rayerr_filter: ::
     *
     *             - if None, on ray error append nothing
     *             - if 'summary', append the exception without ray data
     *             - if 'full', append the exception with ray data up to error
     *             - else append nothing
     *
     *     Returns:
     *         ray_result: see discussion of filters, above.
     *
     */
    /**
     * Trace one pupil coordinate, returning the failure as data rather than
     * letting the exception escape.
     */
    static RayResult trace_safe(optical::OpticalModel *opt_model,
                                const mathlib::Vector2 &pupil, specs::Field &fld,
                                double wvl, const TraceOptions &trace_options);

    /**
     * returns (ray, ray_opl, wvl)
     *
     * Args:
     * seq_model: the :class:`~.SequentialModel` to be traced
     * pt0: starting coordinate at object interface
     * dir0: starting direction cosines following object interface
     * wvl: ray trace wavelength in nm
     * **kwargs: keyword arguments
     *
     * Returns:
     * (**ray**, **op_delta**, **wvl**)
     *
     * - **ray** is a list for each interface in **path_pkg** of these
     * elements: [pt, after_dir, after_dst, normal]
     *
     * - pt: the intersection point of the ray
     * - after_dir: the ray direction cosine following the interface
     * - after_dst: after_dst: the geometric distance to the next
     * interface
     * - normal: the surface normal at the intersection point
     *
     * - **op_delta** - optical path wrt equally inclined chords to the
     * optical axis
     * - **wvl** - wavelength (in nm) that the ray was traced in
     */
    static std::shared_ptr<const RayPkg> trace(seq::SequentialModel *seq_model,
                                               const mathlib::Vector3 &pt0,
                                               const mathlib::Vector3 &dir0, double wvl,
                                               const TraceOptions &trace_options);

    /**
     * Trace ray specified by relative aperture and field point.
     *
     * `pupil_type` controls how `pupil` data is interpreted when calculating the starting ray coordinates.
     *
     * Args:
     * opt_model: instance of :class:`~.OpticalModel` to trace
     * pupil: aperture coordinates of ray
     * fld: instance of :class:`~.Field`
     * wvl: ray trace wavelength in nm
     * apply_vignetting: if True, apply the `fld` vignetting factors to **pupil**
     * pupil_type: ::
     *
     * - 'rel pupil': relative pupil coordinates
     * - 'aim pt': aim point on pupil plane
     * - 'aim dir': aim direction in object space
     *
     * **kwargs: keyword arguments
     *
     * Returns:
     * (**ray**, **op_delta**, **wvl**)
     *
     * - **ray** is a list for each interface in **path_pkg** of these
     * elements: [pt, after_dir, after_dst, normal]
     *
     * - pt: the intersection point of the ray
     * - after_dir: the ray direction cosine following the interface
     * - after_dst: after_dst: the geometric distance to the next
     * interface
     * - normal: the surface normal at the intersection point
     *
     * - **op_delta** - optical path wrt equally inclined chords to the
     * optical axis
     * - **wvl** - wavelength (in nm) that the ray was traced in
     *
     * @param opt_model instance of :class:`~.OpticalModel` to trace
     * @param pupil     relative pupil coordinates of ray
     * @param fld       instance of :class:`~.Field`
     * @param wvl       ray trace wavelength in nm
     */
    static std::shared_ptr<const RayPkg> trace_base(optical::OpticalModel *opt_model,
                                                    const std::vector<double> &pupil,
                                                    specs::Field &fld, double wvl,
                                                    const TraceOptions &trace_options);

    static RayResultWithStartCoord get_1d_solution(seq::SequentialModel *seq_model,
                                                   std::optional<int> ifcx,
                                                   const mathlib::Vector3 &pt0,
                                                   double dist, double wvl,
                                                   double y_target, bool not_wa);

    static RayResultWithStartCoord get_2d_solution(seq::SequentialModel *seq_model,
                                                   std::optional<int> ifcx,
                                                   const mathlib::Vector3 &pt0,
                                                   double dist, double wvl,
                                                   const std::vector<double> &xy_target,
                                                   bool not_wa);

    /**
     * iterates a ray to xy_target on interface ifcx, returns aim points on
     * the paraxial entrance pupil plane
     *
     * If idcx is None, i.e. a floating stop surface, returns xy_target.
     *
     * If the iteration fails, a TraceError will be raised
     *
     */
    static RayResultWithStartCoord iterate_ray(optical::OpticalModel *opt_model,
                                               std::optional<int> ifcx,
                                               const std::vector<double> &xy_target,
                                               specs::Field &fld, double wvl);

    /**
     * returns a list of RayPkgs for the boundary rays for field fld
     */
    static std::vector<std::shared_ptr<const RayPkg>> trace_boundary_rays_at_field(
        optical::OpticalModel *opt_model, specs::Field &fld, double wvl,
        TraceOptions &trace_options);

    static std::map<std::string, std::shared_ptr<const RayPkg>> boundary_ray_dict(
        optical::OpticalModel *opt_model,
        const std::vector<std::shared_ptr<const RayPkg>> &rim_rays);

    static std::vector<std::vector<std::shared_ptr<const RayPkg>>> trace_boundary_rays(
        optical::OpticalModel *opt_model, TraceOptions &trace_options);

    /* returns a list of ray |DataFrame| for the ray_list at field fld */
    static std::vector<RayDataFrame> trace_ray_list_at_field(
        optical::OpticalModel *opt_model, const std::vector<std::vector<double>> &ray_list,
        specs::Field &fld, double wvl, double foc, TraceOptions &trace_options);

    static RayDataFrameByField trace_field(optical::OpticalModel *opt_model,
                                           specs::Field &fld, double wvl, double foc);

    static std::vector<RayDataFrameByField> trace_all_fields(
        optical::OpticalModel *opt_model);

    /**
     * Trace a chief ray at fld and wvl.
     *
     *     Returns:
     *         tuple: **chief_ray**, **cr_exp_seg**
     *
     *             - **chief_ray**: RayPkg of chief ray
     *             - **cr_exp_seg**: exp_pt, exp_dir, exp_dst, ifc, b4_pt, b4_dir
     */
    static std::shared_ptr<const ChiefRayPkg> trace_chief_ray(
        optical::OpticalModel *opt_model, specs::Field &fld, double wvl, double foc);

    static void apply_paraxial_vignetting(optical::OpticalModel *opt_model);

    /**
     * Get the chief ray package at **fld**, computing it if necessary.
     *
     *     Args:
     *         opt_model: :class:`~.OpticalModel` instance
     *         fld: :class:`~.Field` point for wave aberration calculation
     *         wvl: wavelength of ray (nm)
     *         foc: defocus amount
     *
     *     Returns:
     *         tuple: **chief_ray**, **cr_exp_seg**
     *
     *             - **chief_ray**: chief_ray, chief_ray_op, wvl
     *             - **cr_exp_seg**: chief ray exit pupil segment (pt, dir, dist)
     *
     *                 - pt: chief ray intersection with exit pupil plane
     *                 - dir: direction cosine of the chief ray in exit pupil space
     *                 - dist: distance from interface to the exit pupil point
     *
     */
    static std::shared_ptr<const ChiefRayPkg> get_chief_ray_pkg(
        optical::OpticalModel *opt_model, specs::Field &fld, double wvl, double foc);

    /**
     * Trace chief ray and setup reference sphere for `fld`.
     *
     *     Returns:
     *         tuple: **ref_sphere**, **chief_ray_pkg**
     *
     *             - **ref_sphere**: image_pt, ref_dir, ref_sphere_radius, lcl_tfrm_last
     *             - **chief_ray_pkg**: chief_ray, cr_exp_seg
     */
    static RefSphereCR setup_pupil_coords(optical::OpticalModel *opt_model,
                                          specs::Field &fld, double wvl, double foc,
                                          std::optional<mathlib::Vector2> image_pt,
                                          std::optional<mathlib::Vector2> image_delta);

    /** Iterate the chief ray for a field to the stop surface. */
    static AimInfo aim_chief_ray(optical::OpticalModel *opt_model, specs::Field &fld,
                                 std::optional<double> wvl);

    static std::vector<GridItem> trace_fan(optical::OpticalModel *opt_model,
                                           const TraceFanDef &fan_rng, specs::Field &fld,
                                           double wvl, double foc, bool append_if_none,
                                           ImageFilter *img_filter,
                                           const TraceOptions &trace_options);

    static std::vector<GridItem> trace_grid(optical::OpticalModel *opt_model,
                                            const TraceGridDef &grid_rng,
                                            specs::Field &fld, double wvl, double foc,
                                            ImageFilter *img_filter, bool append_if_none,
                                            const TraceOptions &trace_options);

    static std::vector<GridItem> trace_rings(optical::OpticalModel *opt_model,
                                             const TraceRingsDef &grid_rng,
                                             specs::Field &fld, double wvl, double foc,
                                             ImageFilter *img_filter,
                                             bool append_if_none,
                                             const TraceOptions &trace_options);

    static std::vector<GridItem> trace_gaussian_quadrature(
        optical::OpticalModel *opt_model, const TraceRingsDef &grid_rng,
        std::optional<int> num_spokes, specs::Field &fld, double wvl, double foc,
        ImageFilter *img_filter, bool append_if_none,
        const TraceOptions &trace_options);

    /**
     * Trace the reference and the two displaced rays required by contrast
     * optimization. Pupil displacements are expressed in relative pupil
     * coordinates of the physically vignetted pupil. Every quadrature sample
     * is returned so clients can keep a stable residual layout; a ray is null
     * when the trace fails.
     */
    static std::vector<ContrastRayTriplet> trace_contrast(
        optical::OpticalModel *opt_model, const TraceRingsDef &grid_rng,
        std::optional<int> num_spokes, const mathlib::Vector2 &sagittal_shift,
        const mathlib::Vector2 &tangential_shift, specs::Field &fld, double wvl,
        const TraceOptions &trace_options);

    /**
     * Trace contrast triplets, optionally inverse-aiming each displaced partner at its
     * requested coordinate separation on the exit-pupil reference sphere.
     */
    static std::vector<ContrastRayTriplet> trace_contrast(
        optical::OpticalModel *opt_model, const TraceRingsDef &grid_rng,
        std::optional<int> num_spokes, const mathlib::Vector2 &sagittal_shift,
        const mathlib::Vector2 &tangential_shift,
        std::optional<mathlib::Vector2> sagittal_exit_shift,
        std::optional<mathlib::Vector2> tangential_exit_shift, specs::Field &fld,
        double wvl, const TraceOptions &trace_options, bool aim_exit_pupil);

    /**
     * Generate fixed-count quadrature samples in the physical vignetted pupil.
     * The complete pattern is contracted about the overlap centre until every
     * sample and both of its displaced partners are valid. A common contraction
     * preserves relative quadrature weights and avoids frequency-dependent ray
     * loss at the pupil boundary.
     */
    static std::vector<GaussianQuadraturePoint> generate_contrast_quadrature(
        const TraceRingsDef &grid_rng, std::optional<int> num_spokes,
        const mathlib::Vector2 &sagittal_shift,
        const mathlib::Vector2 &tangential_shift, specs::Field &fld);

    static bool inside_vignetted_pupil(const mathlib::Vector2 &pupil,
                                       const specs::Field &fld);

    /**
     * Generates a Gaussian quadrature pattern over a circular or concentric
     * annular pupil.
     *
     * The radial coordinates are Gauss-Legendre nodes transformed from
     * [-1, 1] to squared pupil radius
     * [min_radius^2, max_radius^2]. Each radial node
     * is repeated at uniformly spaced angles. The returned weights are
     * normalized to sum to one, so they integrate a pupil average rather than
     * the area (pi) of the unit disk.
     *
     * Based on B. J. Bauman, H. Xiao, "Gaussian Quadrature for Optical Design
     * with Non-circular Pupils and Fields, and Broad Wavelength Ranges".
     *
     * Also see https://optics.ansys.com/hc/en-us/articles/42661826659347-How-to-use-vignetting-factors
     */
    static std::vector<GaussianQuadraturePoint> generate_gaussian_quadrature(
        const TraceRingsDef &grid_rng, int num_rings, std::optional<int> num_spokes);

    // ---- "raw" variants, which iterate over an explicit path -----------------

    static RayResultWithStartCoord get_1d_solution_raw(
        const std::vector<seq::PathSeg> &pthlist, std::optional<int> ifcx,
        const mathlib::Vector3 &pt0, double dist, double wvl, double y_target,
        bool not_wa);

    static RayResultWithStartCoord get_2d_solution_raw(
        const std::vector<seq::PathSeg> &pthlist, std::optional<int> ifcx,
        const mathlib::Vector3 &pt0, double dist, double wvl,
        const std::vector<double> &xy_target, bool not_wa);

    /**
     * iterates a ray to xy_target on interface ifcx, returns aim points on
     *     the paraxial entrance pupil plane
     *
     *     If idcx is None, i.e. a floating stop surface, returns xy_target.
     *
     *     If the iteration fails, a TraceError will be raised
     */
    static RayResultWithStartCoord iterate_ray_raw(
        const std::vector<seq::PathSeg> &pthlist, std::optional<int> ifcx,
        const std::vector<double> &xy_target, const mathlib::Vector3 &pt0,
        const mathlib::Vector3 &d0, double obj2pup_dist, double eprad, double wvl,
        bool not_wa);

    static void list_ray(std::string &sb, const RayPkg &ray_pkg,
                         const std::optional<math::Tfm3d> &tfrms,
                         std::optional<int> start);
};

} // namespace redukti::rayoptics::raytr

#endif // REDUKTI_RAYOPTICS_RAYTR_TRACE_H
