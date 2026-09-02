// C++ port of org.redukti.rayoptics.elem.profiles.RadialPolynomial
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_ELEM_PROFILES_RADIALPOLYNOMIAL_H
#define REDUKTI_RAYOPTICS_ELEM_PROFILES_RADIALPOLYNOMIAL_H

#include "redukti/rayoptics/elem/profiles/SurfaceProfile.h"
#include "redukti/rayoptics/exceptions/TraceException.h"

#include <optional>

namespace redukti::rayoptics::elem::profiles {

class RadialPolynomial : public SurfaceProfile {
public:
    double ec = 1.0;
    std::vector<double> coefs;
    int max_nonzero_coef = 0;

    /** All five parameters are nullable in the Java. */
    RadialPolynomial(std::optional<double> c, std::optional<double> cc_,
                     std::optional<double> r_, std::optional<double> ec_,
                     const std::vector<double> &coefs_) {
        double ecv = ec_.has_value() ? *ec_ : 1.0;
        double cval = c.has_value() ? *c : 0.0;
        if (r_.has_value()) {
            r(*r_);
        } else {
            cv = cval;
        }
        if (cc_.has_value())
            cc(*cc_);
        else
            this->ec = ecv;
        if (!coefs_.empty()) {
            this->coefs = coefs_;
        } else {
            this->coefs.assign(10, 0.0);
        }
        max_nonzero_coef = 0;
        update();
    }

    RadialPolynomial()
        : RadialPolynomial(std::nullopt, std::nullopt, std::nullopt, std::nullopt, {}) {}

    RadialPolynomial *setCoefs(const std::vector<double> &_coefs) {
        this->coefs = _coefs;
        update();
        return this;
    }

    RadialPolynomial *r(double radius) {
        SurfaceProfile::r(radius);
        return this;
    }

    RadialPolynomial *cc(double cc_) {
        this->ec = cc_ + 1.0;
        return this;
    }

    double cc() const { return ec - 1.0; }

    void calc_max_nonzero_coef() {
        max_nonzero_coef = -1;
        for (std::size_t i = 0; i < coefs.size(); i++) {
            if (coefs[i] != 0.0)
                max_nonzero_coef = static_cast<int>(i);
        }
        max_nonzero_coef++;
    }

    void apply_scale_factor(double scale_factor) override {
        cv /= scale_factor;
        auto sf = 1.0 / scale_factor;
        for (std::size_t i = 0; i < coefs.size(); i++) {
            coefs[i] = std::pow(sf, static_cast<double>(i)) * coefs[i];
        }
    }

    SurfaceProfile *update() override {
        calc_max_nonzero_coef();
        return this;
    }

    double sag(double x, double y) const override {
        double r2 = x * x + y * y;
        double r = std::sqrt(r2);
        // sphere + conic contribution
        double t = 1. - ec * cv * cv * r2;
        if (t < 0)
            throw exceptions::TraceMissedSurfaceException();
        double z = cv * r2 / (1. + std::sqrt(t));
        // polynomial asphere contribution
        double z_asp = 0.0;
        double r_pow = r;
        for (int i = 0; i < max_nonzero_coef; i++) {
            z_asp += coefs[static_cast<std::size_t>(i)] * r_pow;
            r_pow *= r;
        }
        double z_tot = z + z_asp;
        return z_tot;
    }

    double f(const mathlib::Vector3 &p) const override { return p.z - sag(p.x, p.y); }

    mathlib::Vector3 df(const mathlib::Vector3 &p) const override {
        // sphere + conic contribution
        double r2 = p.x * p.x + p.y * p.y;
        double r = std::sqrt(r2);
        double t = 1. - ec * cv * cv * r2;
        if (t < 0)
            throw exceptions::TraceMissedSurfaceException();
        double e = cv / std::sqrt(t);
        // polynomial asphere contribution
        double e_asp = 0.0;
        double r_pow;
        if (r == 0.0)
            r_pow = 1.0;
        else
            // Initialize to 1/r because we multiply by r's components p[0] and
            // p[1] at the final normalization step.
            r_pow = 1 / r;
        double c_coef = 1.0;
        for (int i = 0; i < max_nonzero_coef; i++) {
            e_asp += c_coef * coefs[static_cast<std::size_t>(i)] * r_pow;
            c_coef += 1.0;
            r_pow *= r;
        }
        double e_tot = e + e_asp;
        return mathlib::Vector3(-e_tot * p.x, -e_tot * p.y, 1.0);
    }

    std::vector<mathlib::Vector2> profile(const std::vector<double> &sd, int dir,
                                          int steps) const override {
        (void)sd;
        (void)dir;
        (void)steps;
        return {}; // returns null in the Java
    }

    std::string toString() const override;
};

} // namespace redukti::rayoptics::elem::profiles

#endif // REDUKTI_RAYOPTICS_ELEM_PROFILES_RADIALPOLYNOMIAL_H
