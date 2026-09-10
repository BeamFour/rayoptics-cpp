// C++ port of org.redukti.rayoptics.specs.PupilSpec
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_SPECS_PUPILSPEC_H
#define REDUKTI_RAYOPTICS_SPECS_PUPILSPEC_H

#include "redukti/rayoptics/specs/SpecTypes.h"
#include "redukti/rayoptics/util/Tuples.h"

#include <optional>
#include <string>
#include <vector>

namespace redukti::rayoptics::specs {

class OpticalSpecs;

/**
 * The PupilSpec class maintains the aperture specification.
 * The PupilSpec can be defined in object or image space.
 * The defining parameters can be pupil, f/# or NA,
 * where pupil is the pupil diameter.
 *
 * Attributes:
 * key: 'aperture', 'object'|'image', 'epd'|'NA'|'f/#'
 * value: size of the pupil
 * pupil_rays: list of relative pupil coordinates for pupil limiting rays
 * ray_labels: list of string labels for pupil_rays
 */
class PupilSpec {
public:
    /** Back-reference to the owning OpticalSpecs; borrowed, never owned. */
    OpticalSpecs *optical_spec = nullptr;
    SpecKey key;
    double value = 1.0;
    /**
     * The PupilSpec class allows rays to be specified as fractions of the pupil dimension.
     * A list of pupil_rays and ray_labels define rays to be used to establish clear aperture
     * dimensions on optical elements and rays to be drawn for the lens layout. A default set of
     * pupil rays is provided that is appropriate for circular pupil systems with plane symmetry.
     */
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
     * return pupil spec as paraxial height or slope value
     */
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
