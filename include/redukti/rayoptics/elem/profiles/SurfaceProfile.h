// C++ port of org.redukti.rayoptics.elem.profiles.SurfaceProfile
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_ELEM_PROFILES_SURFACEPROFILE_H
#define REDUKTI_RAYOPTICS_ELEM_PROFILES_SURFACEPROFILE_H

#include "redukti/mathlib/M.h"
#include "redukti/mathlib/Vector2.h"
#include "redukti/mathlib/Vector3.h"
#include "redukti/rayoptics/elem/surface/Aperture.h"
#include "redukti/rayoptics/util/ZDir.h"

#include <string>
#include <vector>

namespace redukti::rayoptics::elem::profiles {

/**
 * Base class for surface profiles.
 *
 * sag() and intersect() can both throw TraceMissedSurfaceException; that is
 * ordinary control flow here, not an error. See the note in TraceException.h.
 */
class SurfaceProfile {
public:
    double cv = 0.0;

    virtual ~SurfaceProfile() = default;

    virtual SurfaceProfile *update() = 0;

    double r() const {
        if (!mathlib::M::isZero(cv))
            return 1.0 / cv;
        else
            return 0.0;
    }

    SurfaceProfile *r(double radius) {
        if (!mathlib::M::isZero(radius))
            cv = 1.0 / radius;
        else
            cv = 0.0;
        return this;
    }

    /** Returns the value of the profile surface function at point p. */
    virtual double f(const mathlib::Vector3 &p) const = 0;

    /** Returns the gradient of the profile surface function at point p. */
    virtual mathlib::Vector3 df(const mathlib::Vector3 &p) const = 0;

    /** Returns the unit normal of the profile at point p. */
    virtual mathlib::Vector3 normal(const mathlib::Vector3 &p) const {
        // TODO check original as that returns p if norm() == 0
        return df(p).normalize();
    }

    /** Returns the sagitta (z coordinate) of the surface at x, y. */
    virtual double sag(double x, double y) const = 0;

    /**
     * Return a 2d polyline approximating the surface profile.
     *
     * @param sd semi-diameter of the profile (array of length 1 or 2)
     * @param dir +1 for profile from neg to positive direction, -1 if otherwise
     * @param steps number of points to generate
     */
    virtual std::vector<mathlib::Vector2> profile(const std::vector<double> &sd, int dir,
                                                  int steps) const = 0;

    /** Apply scale_factor to the profile definition. */
    virtual void apply_scale_factor(double scale_factor) = 0;

    /**
     * Intersect a profile, starting from an arbitrary point.
     *
     * @param p0 start point of the ray in the profile's coordinate system
     * @param d direction cosine of the ray in the profile's coordinate system
     * @param eps numeric tolerance for convergence of any iterative procedure
     * @param z_dir +1 if propagation positive direction, -1 if otherwise
     */
    virtual surface::IntersectionResult intersect(const mathlib::Vector3 &p0,
                                                  const mathlib::Vector3 &d, double eps,
                                                  util::ZDir z_dir) const {
        return intersect_spencer(p0, d, eps, z_dir);
    }

    virtual std::string toString() const { return ""; }

private:
    /**
     * Intersect a profile, starting from an arbitrary point.
     *
     * From Spencer and Murty, General Ray-Tracing Procedure
     * https://doi.org/10.1364/JOSA.52.000672
     */
    surface::IntersectionResult intersect_spencer(const mathlib::Vector3 &p0,
                                                  const mathlib::Vector3 &d, double eps,
                                                  util::ZDir z_dir) const {
        (void)z_dir;
        mathlib::Vector3 p = p0;
        double s1 = -f(p) / d.dot(df(p));
        double delta = std::abs(s1);
        int iter = 0;
        while (delta > eps && iter < 1000) {
            p = p0.add(d.times(s1));
            double s2 = s1 - f(p) / d.dot(df(p));
            delta = std::abs(s2 - s1);
            s1 = s2;
            iter++;
        }
        return surface::IntersectionResult(s1, p);
    }
};

} // namespace redukti::rayoptics::elem::profiles

#endif // REDUKTI_RAYOPTICS_ELEM_PROFILES_SURFACEPROFILE_H
