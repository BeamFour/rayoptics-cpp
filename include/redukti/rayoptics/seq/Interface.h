// C++ port of org.redukti.rayoptics.seq.Interface
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_SEQ_INTERFACE_H
#define REDUKTI_RAYOPTICS_SEQ_INTERFACE_H

#include "redukti/Exceptions.h"
#include "redukti/mathlib/Vector2.h"
#include "redukti/mathlib/Vector3.h"
#include "redukti/rayoptics/elem/profiles/SurfaceProfile.h"
#include "redukti/rayoptics/elem/surface/Aperture.h"
#include "redukti/rayoptics/elem/surface/DecenterData.h"
#include "redukti/rayoptics/seq/Medium.h"
#include "redukti/rayoptics/util/ZDir.h"

#include <cmath>
#include <memory>
#include <optional>
#include <string>

namespace redukti::rayoptics::seq {

/**
 * Basic part of a sequential model.
 *
 * The SequentialModel is a sequence of Interfaces and Gaps. The Interface class
 * is a boundary between two adjacent Gaps and their associated media. It
 * specifies several methods that must be implemented to model the optical
 * behavior of the interface.
 *
 * The Interface class addresses the following use cases:
 *   - support for ray intersection calculation during ray tracing
 *     (interfaces can be tilted and decentered wrt the adjacent gaps)
 *   - support for getting and setting the optical power of the interface
 *   - support for various optical properties, i.e. does it reflect or transmit
 *   - supports a basic idea of size, the max_aperture
 */
class Interface {
public:
    InteractMode interact_mode = InteractMode::TRANSMIT;
    /** refractive index difference across the interface */
    double delta_n = 0.0;
    /** DecenterData for the interface, if specified; null in the Java. */
    std::shared_ptr<elem::surface::DecenterData> decenter;
    /** the maximum aperture radius on the interface */
    double max_aperture = 1.0;
    std::shared_ptr<elem::profiles::SurfaceProfile> profile;

    Interface(std::optional<InteractMode> interact_mode_, double delta_n_, double max_ap,
              std::shared_ptr<elem::surface::DecenterData> decenter_)
        : interact_mode(interact_mode_.has_value() ? *interact_mode_
                                                   : InteractMode::TRANSMIT),
          delta_n(delta_n_), decenter(std::move(decenter_)), max_aperture(max_ap) {
        // TODO phase element
    }

    Interface() = default;

    virtual ~Interface() = default;

    virtual void update() {
        if (decenter)
            decenter->update();
    }

    virtual double profile_cv() const { return 0.0; }

    virtual void set_optical_power(double pwr) { (void)pwr; }

    virtual void set_optical_power(double pwr, double n_before, double n_after) {
        (void)pwr;
        (void)n_before;
        (void)n_after;
    }

    virtual double surface_od() const { throw UnsupportedOperationException(); }

    /**
     * Get a target for ray aiming to aperture boundaries.
     *
     * The main use case is iterating a ray to the internal edge of a surface.
     * Although rel_dir is given as a 2d vector, in practice only the 4 quadrant
     * axes are handled, a 1D directional search along a coordinate axis.
     */
    virtual mathlib::Vector2 edge_pt_target(const mathlib::Vector2 &rel_dir) const {
        return rel_dir.normalize().times(max_aperture);
    }

    /**
     * Returns true if the point (x, y) is inside the clear aperture.
     * `fuzz` is a nullable Double in the Java.
     */
    virtual bool point_inside(double x, double y, std::optional<double> fuzz) const {
        double f = fuzz.has_value() ? *fuzz : 1e-5;
        return std::sqrt(x * x + y * y) <= max_aperture + f;
    }

    virtual void set_max_aperture(double max_ap) { this->max_aperture = max_ap; }

    /** default behavior is returning +/-max_aperture */
    virtual mathlib::Vector2 get_y_aperture_extent() const {
        return mathlib::Vector2(-max_aperture, max_aperture);
    }

    /**
     * Intersect an Interface, starting from an arbitrary point.
     *
     * @param p0 start point of the ray in the interface's coordinate system
     * @param d direction cosine of the ray in the interface's coordinate system
     * @param eps numeric tolerance for convergence of any iterative procedure
     * @param z_dir +1 if propagation positive direction, -1 if otherwise
     *
     * Throws TraceMissedSurfaceException.
     */
    virtual elem::surface::IntersectionResult intersect(const mathlib::Vector3 &p0,
                                                        const mathlib::Vector3 &d,
                                                        double eps,
                                                        util::ZDir z_dir) const {
        (void)p0;
        (void)d;
        (void)eps;
        (void)z_dir;
        throw UnsupportedOperationException();
    }

    elem::surface::IntersectionResult intersect(const mathlib::Vector3 &p0,
                                                const mathlib::Vector3 &d) const {
        return intersect(p0, d, 1.0e-12, util::ZDir::PROPAGATE_RIGHT);
    }

    /** Returns the unit normal of the profile at point p. */
    virtual mathlib::Vector3 normal(const mathlib::Vector3 &p) const {
        (void)p;
        throw UnsupportedOperationException();
    }

    // TODO phase() method

    virtual void apply_scale_factor(double scale_factor) {
        this->max_aperture *= std::abs(scale_factor);
        if (decenter)
            decenter->apply_scale_factor(scale_factor);
    }

    virtual std::string toString() const { return ""; }

    virtual double optical_power() const { return 0.0; }
};

} // namespace redukti::rayoptics::seq

#endif // REDUKTI_RAYOPTICS_SEQ_INTERFACE_H
