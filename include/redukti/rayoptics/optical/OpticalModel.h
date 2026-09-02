// C++ port of org.redukti.rayoptics.optical.OpticalModel
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_OPTICAL_OPTICALMODEL_H
#define REDUKTI_RAYOPTICS_OPTICAL_OPTICALMODEL_H

#include "redukti/rayoptics/specs/SpecTypes.h"

#include <memory>
#include <string>

namespace redukti::rayoptics::seq {
class SequentialModel;
}
namespace redukti::rayoptics::specs {
class OpticalSpecs;
}
namespace redukti::rayoptics::parax {
class ParaxModel;
}

namespace redukti::rayoptics::optical {

/**
 * Top level container for the optical model.
 *
 * The OpticalModel owns its four submodels, and each of them holds an
 * `opt_model` back-pointer to this object. Those back-pointers are raw: the
 * owner always outlives them, so there is no lifetime question, and keeping
 * them raw means every method that reaches back through one keeps its Java
 * signature unchanged.
 *
 * A sequential optical model is a sequence of surfaces and gaps. Additionally
 * it includes optical usage information to specify the aperture, field of view,
 * spectrum and focus.
 *
 * NYI in the Java, and so here: specsheet, ele_model.
 */
class OpticalModel {
public:
    std::unique_ptr<seq::SequentialModel> seq_model;
    std::unique_ptr<specs::OpticalSpecs> optical_spec;
    std::unique_ptr<specs::SystemSpec> system_spec;
    std::unique_ptr<parax::ParaxModel> parax_model;
    bool radius_mode = false;
    std::string dimensions = "mm";

    explicit OpticalModel(bool radius_mode_);
    OpticalModel() : OpticalModel(false) {}
    ~OpticalModel();

    /** Model and its constituents are updated. */
    void update_model();

    /** Compute first order and other optical properties. */
    void update_optical_properties();

    double nm_to_sys_units(double nm) const;

    void apply_scale_factor(double scale_factor);
};

} // namespace redukti::rayoptics::optical

#endif // REDUKTI_RAYOPTICS_OPTICAL_OPTICALMODEL_H
