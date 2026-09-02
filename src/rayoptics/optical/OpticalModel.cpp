// C++ port of org.redukti.rayoptics.optical.OpticalModel
#include "redukti/rayoptics/optical/OpticalModel.h"

#include "redukti/rayoptics/parax/ParaxModel.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"

namespace redukti::rayoptics::optical {

// Construction order matters and follows the Java exactly: SequentialModel
// first (its constructor does not reach back through opt_model), then
// OpticalSpecs, then ParaxModel -- which caches opt_model->seq_model -- and
// only then seq_model->update_model(), which reads optical_spec.
OpticalModel::OpticalModel(bool radius_mode_) : radius_mode(radius_mode_) {
    seq_model = std::make_unique<seq::SequentialModel>(this);
    optical_spec = std::make_unique<specs::OpticalSpecs>(this);
    system_spec = std::make_unique<specs::SystemSpec>();
    parax_model = std::make_unique<parax::ParaxModel>(this, 1.0);
    seq_model->update_model();
}

OpticalModel::~OpticalModel() = default;

void OpticalModel::update_model() {
    seq_model->update_model();
    optical_spec->update_model();
    update_optical_properties();
}

void OpticalModel::update_optical_properties() {
    // OpticalSpec maintains first order and ray aiming for fields
    optical_spec->update_optical_properties();
    // Update the ParaxialModel as needed
    parax_model->update_model();
    // Update surface apertures, if requested (do_apertures=True)
    seq_model->update_optical_properties();
}

double OpticalModel::nm_to_sys_units(double nm) const {
    return system_spec->nm_to_sys_units(nm);
}

void OpticalModel::apply_scale_factor(double scale_factor) {
    seq_model->apply_scale_factor(scale_factor);
    optical_spec->apply_scale_factor(scale_factor);
    optical_spec->update_model();
    update_optical_properties();
}

} // namespace redukti::rayoptics::optical
