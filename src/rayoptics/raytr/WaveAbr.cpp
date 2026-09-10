// C++ port of org.redukti.rayoptics.raytr.WaveAbr
#include "redukti/rayoptics/raytr/WaveAbr.h"

#include "redukti/mathlib/M.h"
#include "redukti/rayoptics/elem/transform/Transform.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/util/Lists.h"

#include <cmath>
#include <limits>

namespace redukti::rayoptics::raytr {

namespace M = mathlib::M;
using elem::transform::Transform;
using mathlib::Vector2;
using mathlib::Vector3;
using util::Lists::get;

/**
 * Compute the reference sphere for a defocussed image point at **fld**.
 *
 *         The local transform from the final interface to the image interface is
 *         included to facilitate infinite reference sphere calculations.
 *
 *     Args:
 *         opt_model: :class:`~.OpticalModel` instance
 *         fld: :class:`~.Field` point for wave aberration calculation
 *         wvl: wavelength of ray (nm)
 *         foc: defocus amount
 *         chief_ray_pkg: input tuple of chief_ray, cr_exp_seg
 *         image_pt_2d: x, y image point in (defocussed) image plane, if None, use
 *                      the chief ray coordinate.
 *         image_delta: x, y displacements from image_pt_2d in (defocussed)
 *                      image plane, if not None.
 *
 *     Returns:
 *         ref_sphere: tuple of image_pt, ref_dir, ref_sphere_radius, lcl_tfrm_last
 */
std::shared_ptr<const ReferenceSphere> WaveAbr::calculate_reference_sphere(
    optical::OpticalModel *opt_model, specs::Field &fld, double wvl, double foc,
    const ChiefRayPkg &chief_ray_pkg, std::optional<Vector2> image_pt_2d,
    std::optional<Vector2> image_delta) {
    (void)fld;
    (void)wvl;
    auto cr = chief_ray_pkg.chief_ray;
    auto cr_exp_seg = chief_ray_pkg.cr_exp_seg;
    Vector3 image_pt = Vector3::ZERO;
    if (!image_pt_2d.has_value()) {
        // get distance along cr corresponding to a z shift of the defocus
        auto dist = foc / get(cr->ray, -1).d.z;
        image_pt = get(cr->ray, -1).p.plus(get(cr->ray, -1).d.times(dist));
    } else {
        image_pt = Vector3(image_pt_2d->x, image_pt_2d->y, foc);
    }
    if (image_delta.has_value())
        image_pt = Vector3(image_pt.x + image_delta->x, image_pt.y + image_delta->y,
                           image_pt.z);
    // get the image point wrt the final surface
    auto seq_model = opt_model->seq_model.get();
    auto lcl_tfrm_last = get(seq_model->lcl_tfrms, -2);
    auto image_thi = get(seq_model->gaps, -1)->thi;
    Vector3 img_pt(image_pt.x, image_pt.y, image_pt.z + image_thi);
    // R' radius of reference sphere for O'
    auto ref_sphere_vec = img_pt.minus(cr_exp_seg->exp_pt);  //p=exp_pt
    auto ref_sphere_radius = ref_sphere_vec.length();
    auto ref_dir = ref_sphere_vec.normalize();
    return std::make_shared<const ReferenceSphere>(image_pt, ref_dir, ref_sphere_radius,
                                                   lcl_tfrm_last);
}

/**
 * Given the exiting interface and chief ray data, return exit pupil ray coords.
 *
 *     Args:
 *         interface: the exiting :class:'~.Interface' for the path sequence, the last surface before image plane
 *         ray_seg: ray segment exiting from **interface**  - this is the ray from last surface to image plane
 *         exp_dst_parax: z distance to the paraxial exit pupil
 *
 *     Returns:
 *         (**exp_pt**, **exp_dir**, **exp_dst**)
 *
 *         - **exp_pt** - ray intersection with exit pupil plane
 *         - **exp_dir** - direction cosine of the ray in exit pupil space
 *         - **exp_dst** - distance from interface to exit pupil pt
 *         - **interface** - exiting :class:'~.Interface' for the path sequence
 *         - **b4_pt** - ray intersection pt wrt image gap coordinates
 *         - **b4_dir** - ray direction cosine wrt image gap coordinates
 */
std::shared_ptr<const ChiefRayExitPupilSegment> WaveAbr::transfer_to_exit_pupil(
    std::shared_ptr<seq::Interface> ifc, const RayData &ray_seg,
    double exp_dst_parax) {
    RayData b4_ray = Transform::transform_after_surface(*ifc, ray_seg);
    Vector3 b4_pt = b4_ray.pt;
    Vector3 b4_dir = b4_ray.dir;
    // h = b4_pt[0]**2 + b4_pt[1]**2
    // u = b4_dir[0]**2 + b4_dir[1]**2
    // handle field points in the YZ plane
    double h = b4_pt.y;  // y=1
    double u = b4_dir.y;  // y=1
    double exp_dst;
    if (std::abs(u) < 1e-14) {
        exp_dst = exp_dst_parax;
    } else {
        // exp_dst = -np.sign(b4_dir[2])*sqrt(h/u)
        exp_dst = -h / u;
    }
    Vector3 exp_pt = b4_pt.plus(b4_dir.times(exp_dst));
    Vector3 exp_dir = b4_dir;
    return std::make_shared<const ChiefRayExitPupilSegment>(exp_pt, exp_dir, exp_dst,
                                                            ifc, b4_pt, b4_dir);
}

/**
 * calculate equally inclined chord distance between 2 rays
 *
 * Args:
 * r: (p, d), where p is a point on the ray r and d is the direction
 * cosine of r
 * r0: (p0, d0), where p0 is a point on the ray r0 and d0 is the direction
 * cosine of r0
 *
 * Returns:
 * float: distance along r from equally inclined chord point to p
 */
double WaveAbr::eic_distance(const RayData &r, const RayData &r0) {
    // eq 3.9 Hopkins paper
    double e = (r.dir.plus(r0.dir).dot(r.pt.minus(r0.pt))) / (1. + r.dir.dot(r0.dir));
    return e;
}

/**
 * compute distance along ray to perpendicular to the origin.
 *
 *     Args:
 *         p, d: a ray, defined by point p and unit direction d
 *
 *     Returns:
 *         t: distance from p to perpendicular to the origin
 */
double WaveAbr::ray_dist_to_perp_from_origin(const RayData &r) {
    auto p = r.pt;
    auto d = r.dir;
    return d.dot(p.negate());
}

util::Pair<util::Pair<Vector3, double>, util::Pair<Vector3, double>>
/**
 * compute distance to pts at the closest join between 2 rays.
 *
 *     Args:
 *         r1: ray 1, defined by point p1 and unit direction d1
 *         r2: ray 2, defined by point p2 and unit direction d2
 *
 *     Returns: (p1_min, t1), (p2_min, t2)
 *         p1_min: point on ray 1 at the closest join
 *         t1: distance from p1 to p1_min
 *         p2_min: point on ray 2 at the closest join
 *         t2: distance from p2 to p2_min
 */
WaveAbr::dist_to_shortest_join(const RayData &r1, const RayData &r2) {
    auto p1 = r1.pt;
    auto d1 = r1.dir;
    auto p2 = r2.pt;
    auto d2 = r2.dir;
    auto del_p = p2.minus(p1);
    auto n = d1.cross(d2);
    auto nn = n.dot(n);
    if (nn == 0.0) {
        auto t2 = p1.minus(p2).dot(d1) * d1.dot(d2);
        auto p2_min = p2.plus(d2.times(t2));
        return util::Pair<util::Pair<Vector3, double>, util::Pair<Vector3, double>>(
            util::Pair<Vector3, double>(p1, 0.0),
            util::Pair<Vector3, double>(p2_min, t2));
    } else {
        auto t1 = d2.cross(n).dot(del_p) / nn;
        auto t2 = d1.cross(n).dot(del_p) / nn;
        auto p1_min = p1.plus(d1.times(t1));
        auto p2_min = p2.plus(d2.times(t2));
        return util::Pair<util::Pair<Vector3, double>, util::Pair<Vector3, double>>(
            util::Pair<Vector3, double>(p1_min, t1),
            util::Pair<Vector3, double>(p2_min, t2));
    }
}

/**
 * Given a ray, a chief ray and an image pt, evaluate the OPD.
 *
 *     The main references for the calculations are in the H. H. Hopkins paper
 *     `Calculation of the Aberrations and Image Assessment for a General Optical
 *     System <https://doi.org/10.1080/713820605>`_
 *
 *     Args:
 *         fod: :class:`~.FirstOrderData` for object and image space refractive
 *              indices
 *         fld: :class:`~.Field` point for wave aberration calculation
 *         wvl: wavelength of ray (nm)
 *         foc: defocus amount
 *         ray_pkg: input tuple of ray, ray_op, wvl
 *         chief_ray_pkg: input tuple of chief_ray, cr_exp_seg
 *         ref_sphere: input tuple of image_pt, ref_dir, ref_sphere_radius
 *
 *     Returns:
 *         opd: OPD of ray wrt chief ray at **fld**
 */
double WaveAbr::wave_abr_full_calc(
    const parax::FirstOrderData &fod, specs::Field &fld, double wvl, double foc,
    const std::shared_ptr<const RayPkg> &ray_pkg,
    const std::shared_ptr<const ChiefRayPkg> &chief_ray_pkg,
    const std::shared_ptr<const ReferenceSphere> &ref_sphere) {
    if (M::is_kinda_big(ref_sphere->ref_sphere_radius))
        return wave_abr_full_calc_inf_ref(fod, fld, wvl, foc, ray_pkg, chief_ray_pkg,
                                          ref_sphere);
    else
        return wave_abr_full_calc_finite_pup(fod, fld, wvl, foc, ray_pkg, chief_ray_pkg,
                                             ref_sphere);
}

/**
 * Given a ray, a chief ray and an image pt, compute the data required to calculate the
 * optical path difference. A by-product of the calculation is the coordinate of the test ray’s
 * reference-sphere intersection relative to the chief-ray exit-pupil point.
 *
 *     The main references for the calculations are in the H. H. Hopkins paper
 *     `Calculation of the Aberrations and Image Assessment for a General Optical
 *     System <https://doi.org/10.1080/713820605>`_
 *
 *     Args:
 *         ray_pkg: input tuple of ray, ray_op, wvl
 *         chief_ray_pkg: input tuple of chief_ray, cr_exp_seg
 *         ref_sphere: input tuple of image_pt, ref_dir, ref_sphere_radius
 *
 *     Returns:
 *         various components required to calculate OPD
 *
 * Note that this is refactored from upstream wave_abr_full_calc_finite_pup because it is
 * useful to be able to reuse the calculation to obtain the aperture ray's
 * intersection with the reference sphere.
 */
FinitePupilWaveAberrationResult WaveAbr::wave_abr_calc_finite_pupil(
    const std::shared_ptr<const RayPkg> &ray_pkg,
    const std::shared_ptr<const ChiefRayPkg> &chief_ray_pkg,
    const std::shared_ptr<const ReferenceSphere> &ref_sphere) {
    Vector3 ref_dir = ref_sphere->ref_dir;
    double ref_sphere_radius = ref_sphere->ref_sphere_radius;
    auto cr = chief_ray_pkg->chief_ray;
    auto cr_exp_seg = chief_ray_pkg->cr_exp_seg;
    const std::vector<RaySeg> &cr_ray = cr->ray;
    double cr_op = cr->op_delta;
    Vector3 cr_exp_pt = cr_exp_seg->exp_pt;
    double cr_exp_dist = cr_exp_seg->exp_dst;
    auto ifc = cr_exp_seg->ifc;
    const std::vector<RaySeg> &ray = ray_pkg->ray;
    double ray_op = ray_pkg->op_delta;
    const int k = -2;  // last interface in sequence
    // eq 3.12
    double e1 = eic_distance(RayData(ray[1].p, ray[0].d),
                             RayData(cr_ray[1].p, cr_ray[0].d));
    // eq 3.13
    double ekp = eic_distance(RayData(get(ray, k).p, get(ray, k).d),
                              RayData(get(cr_ray, k).p, get(cr_ray, k).d));
    RayData tafter = Transform::transform_after_surface(
        *ifc, RayData(get(ray, k).p, get(ray, k).d));
    // cr_exp_pt = Ē′ in HH paper
    // eic_exp_pt = B̃′ in HH paper
    Vector3 b4_pt = tafter.pt;  // Test-ray point at the last optical surface, expressed in image-gap coordinates
    Vector3 b4_dir = tafter.dir;  // Test-ray direction after the last optical surface, expressed in image-gap coordinates
    // -dst = cr_exp_dist - ekp is the signed distance from b4_pt to B̃′ along b4_dir.
    double dst = ekp - cr_exp_dist;
    Vector3 eic_exp_pt = b4_pt.minus(b4_dir.times(dst));  // B̃′: EIC point on the test ray near the exit pupil
    Vector3 p_coord = eic_exp_pt.minus(cr_exp_pt);  // Vector Ē′B̃′ = B̃′ - Ē′
    // eq 4.4
    double F = ref_dir.dot(b4_dir) - b4_dir.dot(p_coord) / ref_sphere_radius;
    // eq 4.5
    double J = p_coord.dot(p_coord) / ref_sphere_radius - 2.0 * ref_dir.dot(p_coord);
    double sign_soln = ref_dir.z * get(cr->ray, -1).d.z < 0 ? -1 : 1;
    // denominator in eq 4.6
    double ep;
    double discriminant = F * F - J / ref_sphere_radius;
    std::optional<Vector3> ray_exit_pupil_coord;
    if (discriminant < 0) {
        ep = std::numeric_limits<double>::quiet_NaN();
    } else {
        double denom = F + sign_soln * std::sqrt(discriminant);
        // Eq 4.6: signed distance e′ from B̃′ along the test ray to its reference-sphere intersection B′.
        ep = denom == 0 ? 0.0 : J / denom;
        ray_exit_pupil_coord = p_coord.plus(b4_dir.times(ep));
    }
    return FinitePupilWaveAberrationResult(ray_pkg, chief_ray_pkg, ref_sphere, e1, ekp,
                                           ep, ray_exit_pupil_coord, ray_op, cr_op);
}

/**
 * Given a ray, a chief ray and an image pt, evaluate the OPD.
 *
 *     The main references for the calculations are in the H. H. Hopkins paper
 *     `Calculation of the Aberrations and Image Assessment for a General Optical
 *     System <https://doi.org/10.1080/713820605>`_
 *
 *     Args:
 *         fod: :class:`~.FirstOrderData` for object and image space refractive
 *              indices
 *         fld: :class:`~.Field` point for wave aberration calculation
 *         wvl: wavelength of ray (nm)
 *         foc: defocus amount
 *         ray_pkg: input tuple of ray, ray_op, wvl
 *         chief_ray_pkg: input tuple of chief_ray, cr_exp_seg
 *         ref_sphere: input tuple of image_pt, ref_dir, ref_sphere_radius
 *
 *     Returns:
 *         opd: OPD of ray wrt chief ray at **fld**
 */
double WaveAbr::wave_abr_full_calc_finite_pup(
    const parax::FirstOrderData &fod, specs::Field &fld, double wvl, double foc,
    const std::shared_ptr<const RayPkg> &ray_pkg,
    const std::shared_ptr<const ChiefRayPkg> &chief_ray_pkg,
    const std::shared_ptr<const ReferenceSphere> &ref_sphere) {
    (void)fld;
    (void)wvl;
    (void)foc;
    FinitePupilWaveAberrationResult result =
        wave_abr_calc_finite_pupil(ray_pkg, chief_ray_pkg, ref_sphere);
    double n_obj = std::abs(fod.n_obj);
    double n_img = std::abs(fod.n_img);
    // OPD = -n_obj e1 + (cr_op - ray_op) + n_img(ekp - ep)
    double opd = -n_obj * result.e1 - result.ray_op + n_img * result.ekp +
                 result.cr_op - n_img * result.ep;
    return opd;
}

/**
 * Given a ray, a chief ray and an image pt, evaluate the OPD.
 *
 *     This set of functions calculates the wavefront aberration using an
 *     infinite reference sphere.
 *     The main references for the calculations are in the paper
 *     `Dependence of the wave-front aberration on the radius of the reference sphere <https://doi.org/10.1364/JOSAA.19.001187>`_ by Antonı́n Mikš.
 *
 *     Args:
 *         fod: :class:`~.FirstOrderData` for object and image space refractive
 *              indices
 *         fld: :class:`~.Field` point for wave aberration calculation
 *         wvl: wavelength of ray (nm)
 *         foc: defocus amount
 *         ray_pkg: input tuple of ray, ray_op, wvl
 *         chief_ray_pkg: input tuple of chief_ray, cr_exp_seg
 *         ref_sphere: input tuple of image_pt, ref_dir, ref_sphere_radius, lcl_tfrm_last
 *
 *     Returns:
 *         opd: OPD of ray wrt chief ray at **fld**
 */
double WaveAbr::wave_abr_full_calc_inf_ref(
    const parax::FirstOrderData &fod, specs::Field &fld, double wvl, double foc,
    const std::shared_ptr<const RayPkg> &ray_pkg,
    const std::shared_ptr<const ChiefRayPkg> &chief_ray_pkg,
    const std::shared_ptr<const ReferenceSphere> &ref_sphere) {
    (void)fld;
    (void)foc;
    auto image_pt = ref_sphere->image_pt;
    auto ref_dir = ref_sphere->ref_dir;
    auto ref_sphere_radius = ref_sphere->ref_sphere_radius;
    auto lcl_tfrm_last = ref_sphere->lcl_tfrm_last;
    (void)ref_dir;
    (void)ref_sphere_radius;
    auto cr = chief_ray_pkg->chief_ray;
    auto &cr_ray = cr->ray;
    auto cr_op = cr->op_delta;
    wvl = cr->wvl;
    auto &ray = ray_pkg->ray;
    auto ray_op = ray_pkg->op_delta;
    wvl = ray_pkg->wvl;
    (void)wvl;
    int k = -2;  // last interface in sequence
    auto n_obj = std::abs(fod.n_obj);
    auto n_img = std::abs(fod.n_img);
    // eq 3.12
    auto e1 = eic_distance(RayData(ray[1].p, ray[0].d),
                           RayData(cr_ray[1].p, cr_ray[0].d));
    // eq 3.13
    auto ekp = eic_distance(RayData(get(ray, k).p, get(ray, k).d),
                            RayData(get(cr_ray, k).p, get(cr_ray, k).d));
    (void)ekp;
    Vector3 p_b4 = Vector3::ZERO;
    Vector3 d_b4 = Vector3::ZERO;
    Vector3 p_cr_b4 = Vector3::ZERO;
    Vector3 d_cr_b4 = Vector3::ZERO;
    // lcl_tfrm_last is a value here, not a reference, so the Java's null check
    // becomes a check on its rotation being present.
    if (lcl_tfrm_last.rt.has_value()) {
        auto rt = *lcl_tfrm_last.rt;
        auto t = lcl_tfrm_last.t;
        p_b4 = rt.multiply(get(ray, k).p.minus(t));
        d_b4 = rt.multiply(get(ray, k).d);
        p_cr_b4 = rt.multiply(get(cr_ray, k).p.minus(t));
        d_cr_b4 = rt.multiply(get(cr_ray, k).d);
    } else {
        p_b4 = get(ray, k).p;
        d_b4 = get(ray, k).d;
        p_cr_b4 = get(cr_ray, k).p;
        d_cr_b4 = get(cr_ray, k).d;
    }
    auto op_b4 = ray_dist_to_perp_from_origin(RayData(p_b4, d_b4));
    auto op_cr_b4 = ray_dist_to_perp_from_origin(RayData(p_cr_b4, d_cr_b4));
    auto P1_P2 = dist_to_shortest_join(
        RayData(get(cr_ray, -1).p, get(cr_ray, -1).d),
        RayData(get(ray, -1).p, get(ray, -1).d));
    auto P1 = P1_P2.first;
    auto P2 = P1_P2.second;
    auto rF0 = (P1.first.plus(P2.first)).divide(2.0);
    auto V_B = ray_op + op_b4;
    auto V_BE = cr_op + op_cr_b4;
    auto W0 = V_B - V_BE + n_img * d_b4.minus(d_cr_b4).dot(rF0);
    auto ta = get(ray, -1).p.minus(image_pt);
    auto numer = d_cr_b4.minus(d_b4.times(d_b4.dot(d_cr_b4))).dot(ta);
    auto denom = 1.0 + d_b4.dot(d_cr_b4);
    auto W_inf = W0 + n_img * numer / denom;
    auto opd = -n_obj * e1 - W_inf;
    return opd;
}

} // namespace redukti::rayoptics::raytr
