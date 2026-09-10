// C++ port of the parax value types:
//   Etendue, ParaxComponent, ParaxPathComp, ParaxData, FirstOrderData,
//   PrincipalPointsInfo, ThirdOrderData.
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_PARAX_PARAXTYPES_H
#define REDUKTI_RAYOPTICS_PARAX_PARAXTYPES_H

#include "redukti/mathlib/M.h"
#include "redukti/rayoptics/seq/Medium.h"

#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

namespace redukti::rayoptics::parax {

/** Conversions between numerical aperture, slope and field angle. */
namespace Etendue {

/**
 * convert numerical aperture to slope
 */
inline double na2slp(double na, double n) { return n * std::tan(std::asin(na / n)); }
inline double na2slp(double na) { return na2slp(na, 1.0); }

/**
 * convert a ray slope to numerical aperture
 */
inline double slp2na(double slp, double n) { return n * std::sin(std::atan(slp / n)); }
inline double slp2na(double slp) { return slp2na(slp, 1.0); }

/**
 * convert paraxial numerical aperture to slope
 */
inline double na2slp_parax(double na, double n) { return na / n; }

/**
 * convert a ray slope to paraxial numerical aperture
 */
inline double slp2na_parax(double slp, double n) { return n * slp; }

/**
 * convert an angle in degrees to a slope
 */
inline double ang2slp(double ang) { return std::tan(mathlib::M::toRadians(ang)); }

/**
 * convert a slope to an angle in degrees
 */
inline double slp2ang(double slp) { return mathlib::M::toDegrees(std::atan(slp)); }

} // namespace Etendue

/** Height, slope and angle of incidence at one surface. */
class ParaxComponent {
public:
    double ht;
    double slp;
    double aoi;

    ParaxComponent(double ht_, double slp_, double aoi_) : ht(ht_), slp(slp_), aoi(aoi_) {}
};

class ParaxPathComp {
public:
    double pwr;
    double tau;
    double indx;
    seq::InteractMode rmd;

    ParaxPathComp(double power, double tau_, double indx_, seq::InteractMode imode)
        : pwr(power), tau(tau_), indx(indx_), rmd(imode) {}
};

/**
 * Container class for first order optical properties
 *
 *     All quantities are based on paraxial ray tracing. The last interface is
 *     the image-1 interface.
 *
 *     Attributes:
 *         opt_inv: optical invariant
 *         efl: effective focal length
 *         pp1: distance of front principle plane from 1st interface
 *         ppk: distance of rear principle plane from last interface
 *         ffl: front focal length
 *         bfl: back focal length
 *         fno: focal ratio at working conjugates, f/#
 *         red: reduction ratio
 *         n_obj: refractive index at central wavelength in object space
 *         n_img: refractive index at central wavelength in image space
 *         obj_dist: object distance
 *         img_dist: paraxial image distance
 *         obj_ang: maximum object angle (degrees)
 *         img_ht: image height
 *         enp_dist: entrance pupil distance from 1st interface
 *         enp_radius: entrance pupil radius
 *         exp_dist: exit pupil distance from last interface
 *         exp_radius: exit pupil radius
 *         obj_na: numerical aperture in object space
 *         img_na: numerical aperture in image space
 */
/** First order properties of the optical model. */
class FirstOrderData {
public:
    /**
     * optical invariant
     */
    double opt_inv = 0.0;
    /**
     * optical power of system
     */
    double power = 0.0;
    /**
     * effective focal length
     */
    double efl = 0.0;
    /**
     * object space focal length, f
     */
    double fl_obj = 0.0;
    /**
     * image space focal length, f'
     */
    double fl_img = 0.0;
    /**
     * distance from the 1st interface to the front principle plane
     */
    double pp1 = 0.0;
    /**
     * distance from the last interface to the rear principle plane
     */
    double ppk = 0.0;
    /**
     * distance from the front principle plane to the rear principle
     * plane
     */
    double pp_sep = 0.0;
    /**
     * front focal length, distance from the 1st interface to the front
     * focal point
     */
    double ffl = 0.0;
    /**
     * back focal length, distance from the last interface to the back
     * focal point
     */
    double bfl = 0.0;
    /**
     * focal ratio at working conjugates, f/#
     */
    double fno = 0.0;
    /**
     * transverse magnification
     */
    double m = 0.0;
    /**
     * reduction ratio, -1/m
     */
    double red = 0.0;
    /**
     * refractive index at central wavelength in object space
     */
    double n_obj = 0.0;
    /**
     * refractive index at central wavelength in image space
     */
    double n_img = 0.0;
    /**
     * object distance
     */
    double obj_dist = 0.0;
    /**
     * paraxial image distance
     */
    double img_dist = 0.0;
    /**
     * maximum object angle (degrees)
     */
    double obj_ang = 0.0;
    /**
     * image height
     */
    double img_ht = 0.0;
    /**
     * entrance pupil distance from 1st interface
     */
    double enp_dist = 0.0;
    /**
     * entrance pupil radius
     */
    double enp_radius = 0.0;
    /**
     * exit pupil distance from last interface
     */
    double exp_dist = 0.0;
    /**
     * exit pupil radius
     */
    double exp_radius = 0.0;
    /**
     * numerical aperture in object space
     */
    double obj_na = 0.0;
    /**
     * numerical aperture in image space
     */
    double img_na = 0.0;

    /** Formatted exactly as the Java, using its %g semantics. */
    void toString(std::string &sb) const;
    std::string toString() const;
};

/**
 * tuple grouping together paraxial rays and first order properties
 *
 *     Attributes:
 *         ax_ray: axial marginal ray data, y, u, i
 *         pr_ray: chief ray data, y, u, i
 *         fod: instance of :class:`~.FirstOrderData`
 */
class ParaxData {
public:
    std::vector<ParaxComponent> ax_ray;
    std::vector<ParaxComponent> pr_ray;
    FirstOrderData fod;

    ParaxData(std::vector<ParaxComponent> ax_ray_, std::vector<ParaxComponent> pr_ray_,
              const FirstOrderData &fod_)
        : ax_ray(std::move(ax_ray_)), pr_ray(std::move(pr_ray_)), fod(fod_) {}
};

class PrincipalPointsInfo {
public:
    std::vector<ParaxComponent> p_ray;
    std::vector<ParaxComponent> q_ray;
    double power;
    double fl_obj;
    double fl_img;
    double efl;
    double pp1;
    double ppk;
    double ffl;
    double bfl;
    double pp_sep;

    PrincipalPointsInfo(std::vector<ParaxComponent> p_ray_,
                        std::vector<ParaxComponent> q_ray_, double power_, double efl_,
                        double fl_obj_, double fl_img_, double pp1_, double ppk_,
                        double pp_sep_, double ffl_, double bfl_)
        : p_ray(std::move(p_ray_)), q_ray(std::move(q_ray_)), power(power_),
          fl_obj(fl_obj_), fl_img(fl_img_), efl(efl_), pp1(pp1_), ppk(ppk_), ffl(ffl_),
          bfl(bfl_), pp_sep(pp_sep_) {}
};

/** Seidel aberration contributions for one surface. */
class ThirdOrderData {
public:
    int c;
    double SI, SII, SIII, SIV, SV;
    double SI_star = 0.0, SII_star = 0.0, SIII_star = 0.0, SIV_star = 0.0,
           SV_star = 0.0;

    ThirdOrderData(int c_, double SI_, double SII_, double SIII_, double SIV_, double SV_)
        : c(c_), SI(SI_), SII(SII_), SIII(SIII_), SIV(SIV_), SV(SV_) {}

    /** Surface index this contribution belongs to. */
    int surface() const { return c; }

    /** Surface contribution as {S-I, S-II, S-III, S-IV, S-V}. */
    std::array<double, 5> seidel() const { return {SI, SII, SIII, SIV, SV}; }

    /**
     * Aspheric contribution as {S-I, S-II, S-III, S-IV, S-V}, all zero unless
     * the surface has a non-zero 4th order aspheric term.
     */
    std::array<double, 5> aspheric() const {
        return {SI_star, SII_star, SIII_star, SIV_star, SV_star};
    }

    /** True when this surface carries a non-zero aspheric contribution. */
    bool has_aspheric() const {
        for (double v : aspheric()) {
            if (v != 0.0)
                return true;
        }
        return false;
    }
};

} // namespace redukti::rayoptics::parax

#endif // REDUKTI_RAYOPTICS_PARAX_PARAXTYPES_H
