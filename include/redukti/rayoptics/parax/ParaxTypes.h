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

inline double na2slp(double na, double n) { return n * std::tan(std::asin(na / n)); }
inline double na2slp(double na) { return na2slp(na, 1.0); }

inline double slp2na(double slp, double n) { return n * std::sin(std::atan(slp / n)); }
inline double slp2na(double slp) { return slp2na(slp, 1.0); }

inline double na2slp_parax(double na, double n) { return na / n; }
inline double slp2na_parax(double slp, double n) { return n * slp; }

inline double ang2slp(double ang) { return std::tan(mathlib::M::toRadians(ang)); }
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

/** First order properties of the optical model. */
class FirstOrderData {
public:
    double opt_inv = 0.0;
    double power = 0.0;
    double efl = 0.0;
    double fl_obj = 0.0;
    double fl_img = 0.0;
    double pp1 = 0.0;
    double ppk = 0.0;
    double pp_sep = 0.0;
    double ffl = 0.0;
    double bfl = 0.0;
    double fno = 0.0;
    double m = 0.0;
    double red = 0.0;
    double n_obj = 0.0;
    double n_img = 0.0;
    double obj_dist = 0.0;
    double img_dist = 0.0;
    double obj_ang = 0.0;
    double img_ht = 0.0;
    double enp_dist = 0.0;
    double enp_radius = 0.0;
    double exp_dist = 0.0;
    double exp_radius = 0.0;
    double obj_na = 0.0;
    double img_na = 0.0;

    /** Formatted exactly as the Java, using its %g semantics. */
    void toString(std::string &sb) const;
    std::string toString() const;
};

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
