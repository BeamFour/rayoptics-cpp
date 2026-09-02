// C++ port of org.redukti.rayoptics.specs.PupilSpec
#ifndef REDUKTI_RAYOPTICS_SPECS_PUPILSPEC_H
#define REDUKTI_RAYOPTICS_SPECS_PUPILSPEC_H

#include "redukti/rayoptics/specs/SpecTypes.h"
#include "redukti/rayoptics/util/Tuples.h"

#include <optional>
#include <string>
#include <vector>

namespace redukti::rayoptics::specs {

class OpticalSpecs;

class PupilSpec {
public:
    /** Back-reference to the owning OpticalSpecs; borrowed, never owned. */
    OpticalSpecs *optical_spec = nullptr;
    SpecKey key;
    double value = 1.0;
    /** Pupil coordinates to trace; each entry has two elements. */
    std::vector<std::vector<double>> pupil_rays;
    std::vector<std::string> ray_labels;

    static const std::vector<std::vector<double>> &default_pupil_rays();
    static const std::vector<std::string> &default_ray_labels();

    PupilSpec(OpticalSpecs *parent, std::optional<util::Pair<ImageKey, ValueKey>> k,
              std::optional<double> value_);

    void update_model();

    void apply_scale_factor(double scale_factor);

    /**
     * The middle element is nullable: Java leaves the derived key null when no
     * branch matches, and FirstOrder relies on that to reject the spec rather
     * than silently taking a wrong branch.
     */
    util::Triple<ImageKey, std::optional<ValueKey>, double> derive_parax_params() const;

    std::string toString() const;

    void list_str(std::string &sb) const;
};

} // namespace redukti::rayoptics::specs

#endif // REDUKTI_RAYOPTICS_SPECS_PUPILSPEC_H
