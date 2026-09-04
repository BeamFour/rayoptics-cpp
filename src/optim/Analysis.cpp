// C++ port of org.redukti.optim.Analysis
#include "redukti/optim/Analysis.h"

#include "redukti/Exceptions.h"
#include "redukti/Text.h"
#include "redukti/optim/ParaxHelper.h"
#include "redukti/rayoptics/parax/ParaxTypes.h"
#include "redukti/rayoptics/specs/FieldSpec.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/rayoptics/specs/PupilSpec.h"
#include "redukti/rayoptics/raytr/RayTypes.h"

#include <cmath>
#include <iostream>

namespace redukti::optim {

using rayoptics::analysis::ContrastAnalysis;
using rayoptics::analysis::ContrastOptions;
using rayoptics::analysis::SpotAnalysis;
using rayoptics::analysis::SpotOptions;
using rayoptics::analysis::TransverseRayAberrationAnalysis;
using rayoptics::optical::OpticalModel;
using spec::VigType;

Analysis::Analysis(spec::Prescription *prescription, std::vector<double> fields,
                   std::vector<int> freqs, int scenario)
    : _prescription(prescription), _fields(std::move(fields)), _freqs(std::move(freqs)),
      _scenario(scenario) {}

Analysis &Analysis::using_gauss_quadrature_pattern(int num_rings, int num_spokes,
                                                   double innerPupilRadius) {
    if (num_rings < 1 || num_spokes < 3)
        throw IllegalArgumentException(
            "Gaussian quadrature requires at least 1 ring and 3 spokes");
    if (!std::isfinite(innerPupilRadius) || innerPupilRadius < 0.0 ||
        innerPupilRadius >= 1.0)
        throw IllegalArgumentException("Inner pupil radius must be finite and in [0, 1)");
    _spot_pattern = SpotOptions::PATTERN_GAUSS_QUADRATURE;
    _num_rings = num_rings;
    _num_spokes = num_spokes;
    _inner_pupil_radius = innerPupilRadius;
    return *this;
}

std::unique_ptr<OpticalModel> Analysis::build_model(VigType vigType) {
    // The Java hands the builder its live prescription; the C++ builder takes a
    // copy, which is equivalent here because it reads it during this call only.
    return spec::RayOpticsModelBuilder(*_prescription)
        .build_optical_model(true, _fields, false, vigType, true, _scenario);
}

std::unique_ptr<OpticalModel> Analysis::build_vignetted_model() {
    if (!_freeze_vignetting)
        return build_model(_vig_type);
    if (!_frozen_vignetting.has_value())
        capture_vignetting();
    // Build without establishing vignetting, then stamp the captured state on.
    auto model = build_model(VigType::None);
    auto &fields = model->optical_spec->fov->fields;
    if (fields.size() != _frozen_vignetting->size())
        throw IllegalStateException(
            "frozen vignetting was captured for " +
            std::to_string(_frozen_vignetting->size()) + " fields but the model has " +
            std::to_string(fields.size()));
    if (_frozen_pupil_value.has_value() &&
        *_frozen_pupil_value != model->optical_spec->pupil->value) {
        model->optical_spec->pupil->value = *_frozen_pupil_value;
        model->update_model();
    }
    for (std::size_t i = 0; i < fields.size(); i++) {
        fields[i]->vux = (*_frozen_vignetting)[i][0];
        fields[i]->vlx = (*_frozen_vignetting)[i][1];
        fields[i]->vuy = (*_frozen_vignetting)[i][2];
        fields[i]->vly = (*_frozen_vignetting)[i][3];
    }
    return model;
}

void Analysis::capture_vignetting() {
    auto reference = build_model(_vig_type);
    auto &fields = reference->optical_spec->fov->fields;
    std::vector<std::vector<double>> captured(fields.size(), std::vector<double>(4, 0.0));
    for (std::size_t i = 0; i < fields.size(); i++) {
        captured[i][0] = fields[i]->vux;
        captured[i][1] = fields[i]->vlx;
        captured[i][2] = fields[i]->vuy;
        captured[i][3] = fields[i]->vly;
    }
    _frozen_vignetting = std::move(captured);
    _frozen_pupil_value = reference->optical_spec->pupil->value;
}

void Analysis::compute() {
    _opt_model = build_vignetted_model();
    _pfo = ParaxHelper::asArray(_opt_model->optical_spec->parax_data->fod);
    if (_compute_spots) {
        SpotOptions options;
        if (_spot_pattern == SpotOptions::PATTERN_GAUSS_QUADRATURE) {
            options.use_gaussian_quadrature()
                .num_rings(_num_rings)
                .num_spokes(_num_spokes)
                .inner_pupil_radius(_inner_pupil_radius)
                .append_failed_rays(_append_failed_spot_rays)
                .check_apertures(_check_spot_apertures);
        } else {
            options.use_hexapolar().num_rays(_num_rays).check_apertures(
                _check_spot_apertures);
        }
        auto spotAnalysis = SpotAnalysis::eval(_opt_model.get(), options);
        _mtfs = _compute_mtf
                    ? std::optional<std::vector<rayoptics::analysis::MTFResultByFreq>>(
                          spotAnalysis.computeMTFs(_freqs))
                    : std::nullopt;
        // Moved, not copied, and moved after computeMTFs which reads it.
        //
        // Each SpotIntercepts borrows a TraceGridByWvl out of its own
        // SpotResultsForField::trace_results. Moving the outer vector steals the
        // element buffer and leaves every element at its old address, so those
        // pointers stay good; copying one would give it a fresh trace_results
        // buffer while the intercepts still pointed at the original.
        _spots = std::move(spotAnalysis.spot_results);
    } else {
        _spots.reset();
        _mtfs.reset();
    }
    // We set append_if_none to retain failed fan rays; their goals apply a penalty.
    if (_compute_ray_aberrations)
        _ray_aberrations = TransverseRayAberrationAnalysis::eval(
            _opt_model.get(), NUM_TRANSVERSE_RAYS, true, rayoptics::raytr::TraceOptions());
    else
        _ray_aberrations.reset();
    _contrasts.clear();
    _contrasts.reserve(_contrast_freqs.size());
    for (std::size_t i = 0; i < _contrast_freqs.size(); i++) {
        ContrastOptions options(_contrast_freqs[i]);
        options.num_rings(_contrast_num_rings)
            .num_spokes(_contrast_num_spokes)
            .calibrate_frequency(_contrast_calibrate_frequency)
            .aim_exit_pupil(_contrast_aim_exit_pupil)
            .center_residuals(_contrast_center_residuals);
        _contrasts.push_back(ContrastAnalysis::eval(_opt_model.get(), options));
    }
}

} // namespace redukti::optim
