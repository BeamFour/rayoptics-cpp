// C++ port of org.redukti.tools.DefaultOptimizations
#ifndef REDUKTI_TOOLS_DEFAULTOPTIMIZATIONS_H
#define REDUKTI_TOOLS_DEFAULTOPTIMIZATIONS_H

#include "redukti/spec/Prescription.h"

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
 * between them.
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
     * The airspaces a zoom varies between configurations, excluding the back
     * focus. These are the surfaces whose thickness came from a multi valued
     * row of [variable distances].
     */
    static std::vector<int> findVariableThicknesses(const spec::Prescription &prescription,
                                                    int backFocusSurface);

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

    static bool isVariableAirspace(const spec::SurfaceType &surface, int index,
                                   int backFocusSurface);
};

} // namespace redukti::tools

#endif // REDUKTI_TOOLS_DEFAULTOPTIMIZATIONS_H
