// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.tools.DefaultOptimizations
#ifndef REDUKTI_TOOLS_DEFAULTOPTIMIZATIONS_H
#define REDUKTI_TOOLS_DEFAULTOPTIMIZATIONS_H

#include "redukti/spec/Prescription.h"

#include <optional>
#include <vector>

namespace redukti::tools {

/**
 * The two routine optimizations that a prescription imported from a patent
 * almost always needs, wired up so that LensTool2 can run them without the
 * caller having to build a merit function by hand.
 *
 * Both target MTF at the central field, which is what the manual runs this
 * replaces have used and what tends to converge. Effective focal length and
 * f-number are anchored to the prescription by OptimizationBuilder itself, so
 * optimizing an airspace cannot quietly turn the lens into a different one.
 *
 * Only one configuration is optimized per call. There is no support for a
 * variable shared across configurations, so a zoom is optimized one
 * configuration at a time and nothing here tries to hold the back focus common
 * between them - instead the airspaces that have to stay common are left out of
 * the solve, which findVariableThicknesses explains.
 */
class DefaultOptimizations {
public:
    /** Which objective the solve is driven by. */
    enum class Objective { CONTRAST, MTF };

    /** Java's record Result. */
    struct Result {
        std::vector<int> surfaces;
        int status;
        double before;
        double after;

        bool improved() const { return after < before; }
    };

    /**
     * Locates the airspace that acts as the back focus.
     *
     * Normally that is the gap after the last surface. Where the design ends in
     * a cover glass it is the gap in front of the cover glass instead: the short
     * gap between cover glass and image is fixed by the sensor stack and moving
     * it is not what "adjust the back focus" means.
     *
     * @return the surface whose thickness is the back focus, or -1 if the
     *         prescription does not end in an airspace that can be varied
     */
    static int findBackFocusSurface(const spec::Prescription &prescription);

    /**
     * The airspaces a zoom varies between configurations. These are the surfaces
     * whose thickness came from a multi valued row of [variable distances] and
     * whose value actually differs between the configurations being reported.
     *
     * That second half of the test is what keeps the back focus out, and it is
     * deliberately a test on the value rather than on the position. A gap holding
     * the same value at every zoom setting is one that has to stay common - the
     * cover glass to sensor spacing, typically - and since a variable cannot yet
     * be shared across configurations, varying it per configuration would hand
     * each zoom setting its own sensor position.
     *
     * Excluding findBackFocusSurface instead would be wrong on a design
     * that ends in a cover glass: there the airspace in front of the cover glass is
     * the one that moves with the zoom, so that rule points at a genuine zoom
     * airspace, and dropping it would leave the rear group's spacing out of the
     * solve while admitting the fixed sensor gap behind the cover glass.
     */
    static std::vector<int> findVariableThicknesses(const spec::Prescription &prescription);

    /**
     * Runs the solver over the given surfaces' thicknesses for one
     * configuration, targeting central field MTF at the requested frequencies.
     * The prescription is updated in place when the solve improves it.
     */
    static Result optimizeThicknesses(spec::Prescription *prescription,
                                      const std::vector<int> &surfaces,
                                      const std::vector<int> &mtfFrequencies,
                                      int configuration, spec::VigType vigType,
                                      bool dLineOnly, Objective objective);

private:
    DefaultOptimizations() = delete;

    /**
     * True when the thickness row this surface references holds more than one distinct
     * value across the configurations.
     *
     * The array is indexed by configuration, not by the scenario numbering of the
     * input: Prescription fills it from the scenarios the report's
     * scenarios list selects, so the scenarios a patent tabulated but this
     * report does not use are already absent and cannot make an airspace look as though
     * it varies. A value those selected scenarios are missing is rejected while the
     * prescription is built, so nothing here has to allow for one.
     */
    static bool variesByConfiguration(
        const std::optional<std::vector<double>> &thickness_by_scenario);

    static bool isVariableAirspace(const spec::SurfaceType &surface);
};

} // namespace redukti::tools

#endif // REDUKTI_TOOLS_DEFAULTOPTIMIZATIONS_H
