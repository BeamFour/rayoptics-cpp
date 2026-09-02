// C++ port of org.redukti.rayoptics.elem.profiles.EvenPolynomial
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_ELEM_PROFILES_EVENPOLYNOMIAL_H
#define REDUKTI_RAYOPTICS_ELEM_PROFILES_EVENPOLYNOMIAL_H

#include "redukti/rayoptics/elem/profiles/SurfaceProfile.h"
#include "redukti/rayoptics/exceptions/TraceException.h"

#include <optional>

namespace redukti::rayoptics::elem::profiles {

class EvenPolynomial : public SurfaceProfile {
public:
    double cc = 0.0;
    std::vector<double> coefs;
    int max_nonzero_coef = 0;

    /** `r` and `ec` are nullable Doubles in the Java. */
    EvenPolynomial(double c, double cc_, std::optional<double> r_,
                   std::optional<double> ec_, const std::vector<double> &coefs_) {
        if (r_.has_value()) {
            r(*r_);
        } else {
            cv = c;
        }
        if (ec_.has_value())
            ec(*ec_);
        else
            this->cc = cc_;
        if (!coefs_.empty()) {
            this->coefs = coefs_;
        } else {
            this->coefs.assign(10, 0.0);
        }
        max_nonzero_coef = 0;
        update();
    }

    EvenPolynomial() { this->coefs.assign(10, 0.0); }

    EvenPolynomial *c(double _c) {
        cv = _c;
        return this;
    }

    EvenPolynomial *r(double radius) {
        SurfaceProfile::r(radius);
        return this;
    }

    EvenPolynomial *setCoefs(const std::vector<double> &_coefs) {
        this->coefs = _coefs;
        update();
        return this;
    }

    double ec() const { return cc + 1.0; }
    void ec(double ec_) { cc = ec_ - 1.0; }

    EvenPolynomial *setCc(double _cc) {
        this->cc = _cc;
        return this;
    }

    SurfaceProfile *update() override {
        calc_max_nonzero_coef();
        return this;
    }

    double f(const mathlib::Vector3 &p) const override { return p.z - sag(p.x, p.y); }

    mathlib::Vector3 df(const mathlib::Vector3 &p) const override {
        // sphere + conic contribution
        double r2 = p.x * p.x + p.y * p.y;
        double t = 1. - ec() * cv * cv * r2;
        if (t < 0)
            throw exceptions::TraceMissedSurfaceException();
        double e = cv / std::sqrt(t);
        // polynomial asphere contribution
        double r_pow = 1.0;
        double e_asp = 0.0;
        double c_coef = 2.0;
        for (int i = 0; i < max_nonzero_coef; i++) {
            e_asp += c_coef * coefs[static_cast<std::size_t>(i)] * r_pow;
            c_coef += 2.0;
            r_pow *= r2;
        }
        double e_tot = e + e_asp;
        return mathlib::Vector3(-e_tot * p.x, -e_tot * p.y, 1.0);
    }

    double sag(double x, double y) const override {
        double r2 = x * x + y * y;
        // sphere + conic contribution
        double t = 1. - (cc + 1.0) * cv * cv * r2;
        if (t < 0)
            throw exceptions::TraceMissedSurfaceException();
        double z = cv * r2 / (1. + std::sqrt(t));
        // polynomial asphere contribution
        double z_asp = 0.0;
        double r_pow = r2;
        for (int i = 0; i < max_nonzero_coef; i++) {
            z_asp += coefs[static_cast<std::size_t>(i)] * r_pow;
            r_pow *= r2;
        }
        double z_tot = z + z_asp;
        return z_tot;
    }

    double get_by_order(int i) const { return coefs[static_cast<std::size_t>(i / 2 - 1)]; }

    std::vector<mathlib::Vector2> profile(const std::vector<double> &sd, int dir,
                                          int steps) const override {
        (void)sd;
        (void)dir;
        (void)steps;
        return {}; // returns null in the Java
    }

    void apply_scale_factor(double scale_factor) override {
        cv /= scale_factor;
        auto sf = 1.0 / scale_factor;
        for (std::size_t i = 0; i < coefs.size(); i++) {
            coefs[i] = std::pow(sf, 2 * static_cast<double>(i) + 1) * coefs[i];
        }
    }

    void calc_max_nonzero_coef() {
        max_nonzero_coef = -1;
        for (std::size_t i = 0; i < coefs.size(); i++) {
            if (coefs[i] != 0.0)
                max_nonzero_coef = static_cast<int>(i);
        }
        max_nonzero_coef++;
    }

    std::string toString() const override;
};

} // namespace redukti::rayoptics::elem::profiles

#endif // REDUKTI_RAYOPTICS_ELEM_PROFILES_EVENPOLYNOMIAL_H
