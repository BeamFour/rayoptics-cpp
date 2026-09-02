// C++ port of org.redukti.rayoptics.elem.profiles.Spherical
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_ELEM_PROFILES_SPHERICAL_H
#define REDUKTI_RAYOPTICS_ELEM_PROFILES_SPHERICAL_H

#include "redukti/rayoptics/elem/profiles/SurfaceProfile.h"
#include "redukti/rayoptics/exceptions/TraceException.h"

namespace redukti::rayoptics::elem::profiles {

/** Spherical surface profile parameterized by curvature. */
class Spherical : public SurfaceProfile {
public:
    explicit Spherical(double c) { this->cv = c; }
    Spherical() : Spherical(0.0) {}

    SurfaceProfile *update() override { return nullptr; }

    /**
     * surface function for Spherical profile
     *
     * This function implements Spencer's eq 25, with kappa=1 (i.e. spherical).
     *
     * To see this, start with the code:
     *   F = p[2] - 0.5*cv*(np.dot(p, p))
     * Expand np.dot(p, p):
     *   F = p[2] - 0.5*cv*(p[0]*p[0] + p[1]*p[1] + p[2]*p[2])
     * in Spencer's notation rho**2 = p[0]*p[0] + p[1]*p[1] and Z = p[2], so
     *   F = Z - 0.5*cv*(rho**2 + Z**2)
     * which is Spencer's eq 25.
     */
    double f(const mathlib::Vector3 &p) const override {
        return p.z - 0.5 * cv * p.dot(p);
    }

    mathlib::Vector3 df(const mathlib::Vector3 &p) const override {
        return mathlib::Vector3(-cv * p.x, -cv * p.y, 1.0 - cv * p.z);
    }

// MSVC cannot see the isZero guard through the inlined call and warns about a
// potential division by zero. The guard is real: isZero(cv) is false only when
// |cv| > EPSILON.
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4723)
#endif
    double sag(double x, double y) const override {
        if (!mathlib::M::isZero(cv)) {
            double r = 1.0 / cv; // radius = 1/curvature
            double adj = r * r - x * x - y * y;
            if (adj < 0.0)
                throw exceptions::TraceMissedSurfaceException();
            adj = std::sqrt(adj);
            return r * (1.0 - std::abs(adj / r));
        } else {
            return 0.0;
        }
    }
#if defined(_MSC_VER)
#pragma warning(pop)
#endif

    std::vector<mathlib::Vector2> profile(const std::vector<double> &sd, int dir,
                                          int steps) const override {
        (void)sd;
        (void)dir;
        (void)steps;
        return {}; // returns null in the Java
    }

    void apply_scale_factor(double scale_factor) override {
        if (mathlib::M::isZero(scale_factor))
            return;
        cv /= scale_factor;
    }

    /**
     * Intersection with a sphere, starting from an arbitrary point.
     *
     * Substitute expressions equivalent to Welford's 4.8 and 4.9. For the
     * quadratic ax**2 + bx + c = 0: ax2 = 2a, cx2 = 2c.
     */
    surface::IntersectionResult intersect(const mathlib::Vector3 &p,
                                          const mathlib::Vector3 &d, double eps,
                                          util::ZDir z_dir) const override {
        (void)eps;
        double ax2 = cv;
        double cx2 = cv * p.dot(p) - 2.0 * p.z;
        double b = cv * d.dot(p) - d.z;
        double s = 0.0;
        if ((b != 0) || (cx2 != 0) || (ax2 != 0)) {
            double tmp = b * b - ax2 * cx2;
            if (tmp < 0)
                throw exceptions::TraceMissedSurfaceException();
            s = cx2 / (util::value(z_dir) * std::sqrt(tmp) - b);
        }
        // else ax2 = cx2 = b = 0, i.e. ray is tangent to the sphere at p
        mathlib::Vector3 p1 = p.add(d.times(s));
        return surface::IntersectionResult(s, p1);
    }

    std::string toString() const override;
};

} // namespace redukti::rayoptics::elem::profiles

#endif // REDUKTI_RAYOPTICS_ELEM_PROFILES_SPHERICAL_H
