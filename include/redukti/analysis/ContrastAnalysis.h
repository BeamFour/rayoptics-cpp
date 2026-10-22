// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.analysis.ContrastOptions,
// ContrastAnalysisResult, ContrastAnalysis, ContrastRayTriplet, ContrastTraceByWvl
// and ContrastTraceCallback.
#ifndef REDUKTI_ANALYSIS_CONTRASTANALYSIS_H
#define REDUKTI_ANALYSIS_CONTRASTANALYSIS_H

#include "redukti/analysis/ExitPupilAiming.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/raytr/RayTrace.h"
#include "redukti/rayoptics/raytr/RayTypes.h"
#include "redukti/rayoptics/raytr/Trace.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"

#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace redukti::analysis {

/** Three rays used by a contrast-optimization pupil sample. */
/** Java's `record ContrastRayTriplet(...)`. */
class ContrastRayTriplet {
public:
    mathlib::Vector2 pupil;
    std::shared_ptr<const rayoptics::raytr::RayPkg> reference;
    std::shared_ptr<const rayoptics::raytr::RayPkg> sagittal;
    std::shared_ptr<const rayoptics::raytr::RayPkg> tangential;
    std::shared_ptr<rayoptics::exceptions::TraceException> referenceError;
    std::shared_ptr<rayoptics::exceptions::TraceException> sagittalError;
    std::shared_ptr<rayoptics::exceptions::TraceException> tangentialError;
    double weight;

    ContrastRayTriplet(const mathlib::Vector2 &pupil_,
                       std::shared_ptr<const rayoptics::raytr::RayPkg> reference_,
                       std::shared_ptr<const rayoptics::raytr::RayPkg> sagittal_,
                       std::shared_ptr<const rayoptics::raytr::RayPkg> tangential_,
                       std::shared_ptr<rayoptics::exceptions::TraceException> referenceError_,
                       std::shared_ptr<rayoptics::exceptions::TraceException> sagittalError_,
                       std::shared_ptr<rayoptics::exceptions::TraceException> tangentialError_,
                       double weight_)
        : pupil(pupil_), reference(std::move(reference_)), sagittal(std::move(sagittal_)),
          tangential(std::move(tangential_)), referenceError(std::move(referenceError_)),
          sagittalError(std::move(sagittalError_)),
          tangentialError(std::move(tangentialError_)), weight(weight_) {}
};

/** Java's `record ContrastTraceByWvl<T>(double wavelength, List<T> samples)`. */
template <typename T> class ContrastTraceByWvl {
public:
    double wavelength;
    std::vector<T> samples;

    ContrastTraceByWvl(double wavelength_, std::vector<T> samples_)
        : wavelength(wavelength_), samples(std::move(samples_)) {}
};

/** Java's `interface ContrastTraceCallback<T>`. */
template <typename T>
using ContrastTraceCallback = std::function<T(
    const ContrastRayTriplet &rays, rayoptics::specs::Field &field, double wavelength, double focus)>;

/** Options for pupil-autocorrelation contrast analysis. */
/** Fluent configuration for a contrast (through-frequency) analysis. */
class ContrastOptions {
public:
    double spatialFrequency;
    int numRings = 3;
    std::optional<int> numSpokes = 6;
    rayoptics::raytr::TraceOptions traceOptions;
    bool calibrateFrequency = false;
    bool aimExitPupil = false;
    bool centerResiduals = false;

    explicit ContrastOptions(double spatialFrequency_);

    ContrastOptions &num_rings(int value);
    ContrastOptions &num_spokes(std::optional<int> value);
    ContrastOptions &trace_options(const rayoptics::raytr::TraceOptions &value);
    ContrastOptions &center_residuals(bool value);
    ContrastOptions &check_apertures(bool value);
    ContrastOptions &calibrate_frequency(bool value);
    ContrastOptions &aim_exit_pupil(bool value);
};

class ContrastAnalysisResult {
public:
    /** Java record Failure(String ray, String exceptionType, int surface). */
    class Failure {
    public:
        std::string ray;
        std::string exceptionType;
        int surface;

        Failure(std::string ray_, std::string exceptionType_, int surface_)
            : ray(std::move(ray_)), exceptionType(std::move(exceptionType_)),
              surface(surface_) {}
    };

    /** Java record Sample(...). */
    class Sample {
    public:
        mathlib::Vector2 pupil;
        double sagittalDifference;
        double tangentialDifference;
        double weight;
        bool valid;
        /** Null when the sample is good. */
        std::optional<Failure> failure;

        Sample(const mathlib::Vector2 &pupil_, double sagittalDifference_,
               double tangentialDifference_, double weight_, bool valid_,
               std::optional<Failure> failure_)
            : pupil(pupil_), sagittalDifference(sagittalDifference_),
              tangentialDifference(tangentialDifference_), weight(weight_),
              valid(valid_), failure(std::move(failure_)) {}

        double sagittalResidual() const;
        double tangentialResidual() const;
    };

    /**
     * One (frequency, field, wavelength) block of samples, with the constant part of the
     * wavefront difference that is to be removed from each orientation.
     *
     * The offsets are zero unless ContrastOptions#center_residuals(boolean) is
     * enabled, in which case they hold the <em>reference wavelength's</em> weighted mean
     * difference for this field. See
     * ContrastAnalysis#center_residuals(ContrastAnalysisResult, int).
     */
    /** Java record WavelengthResult(...). */
    class WavelengthResult {
    public:
        double wavelength;
        double normalizedPupilShift;
        std::vector<Sample> samples;
        double sagittalOffset;
        double tangentialOffset;

        WavelengthResult(double wavelength_, double normalizedPupilShift_,
                         std::vector<Sample> samples_, double sagittalOffset_,
                         double tangentialOffset_)
            : wavelength(wavelength_), normalizedPupilShift(normalizedPupilShift_),
              samples(std::move(samples_)), sagittalOffset(sagittalOffset_),
              tangentialOffset(tangentialOffset_) {}

        WavelengthResult(double wavelength_, double normalizedPupilShift_,
                         std::vector<Sample> samples_)
            : WavelengthResult(wavelength_, normalizedPupilShift_, std::move(samples_),
                               0.0, 0.0) {}

        WavelengthResult withOffsets(double sagittal, double tangential) const;

        double sagittalResidual(int index) const;
        double tangentialResidual(int index) const;
    };

    /** Java record FieldResult(Field field, List<WavelengthResult> wavelengths). */
    class FieldResult {
    public:
        std::shared_ptr<const rayoptics::specs::FieldSnapshot> field;
        std::vector<WavelengthResult> wavelengths;

        FieldResult(rayoptics::specs::Field *field_, std::vector<WavelengthResult> wavelengths_)
            : field(field_ ? std::make_shared<const rayoptics::specs::FieldSnapshot>(*field_) : nullptr),
              wavelengths(std::move(wavelengths_)) {}
    };

    /**
     * @param spatialFrequency image-space spatial frequency in cycles per
     *                         optical-system length unit (normally cycles/mm)
     */
    double spatialFrequency;
    std::vector<FieldResult> fields;

    explicit ContrastAnalysisResult(double spatialFrequency_)
        : spatialFrequency(spatialFrequency_) {}
};

/** Wavefront-difference analysis used by contrast optimization. */
class ContrastAnalysis {
public:
    /**
     * Trace contrast-optimization ray triplets and apply the analysis callback
     * while the wavelength-specific chief ray and reference sphere are active.
     */
    template <typename T>
    static std::vector<ContrastTraceByWvl<T>> trace_by_wavelength(
        rayoptics::optical::OpticalModel *opt_model, const ContrastTraceCallback<T> &callback, int fi, std::optional<int> wl,
        int num_rings, std::optional<int> num_spokes,
        const mathlib::Vector2 &sagittal_shift, const mathlib::Vector2 &tangential_shift,
        const rayoptics::raytr::TraceOptions &trace_options) {
        return trace_by_wavelength<T>(opt_model, callback, fi, wl, num_rings, num_spokes,
                                      sagittal_shift, tangential_shift, 0.0, trace_options,
                                      false);
    }

    template <typename T>
    static std::vector<ContrastTraceByWvl<T>> trace_by_wavelength(
        rayoptics::optical::OpticalModel *opt_model, const ContrastTraceCallback<T> &callback, int fi, std::optional<int> wl,
        int num_rings, std::optional<int> num_spokes,
        const mathlib::Vector2 &sagittal_shift, const mathlib::Vector2 &tangential_shift,
        double spatial_frequency, const rayoptics::raytr::TraceOptions &trace_options,
        bool aim_exit_pupil) {
        auto osp = opt_model->optical_spec.get();
        const auto &wavelengths = osp->wvls->wavelengths;
        std::vector<double> wavelengthList;
        if (!wl.has_value())
            wavelengthList = wavelengths;
        else
            wavelengthList = {wavelengths[static_cast<std::size_t>(*wl)]};
        rayoptics::specs::Field &field = *osp->fov->fields[static_cast<std::size_t>(fi)];
        auto focus = osp->defocus()->get_focus();
        auto referenceWavelength = opt_model->seq_model->central_wavelength();
        auto referenceCoordinates = rayoptics::raytr::Trace::setup_pupil_coords(
            opt_model, field, referenceWavelength, focus, std::nullopt, std::nullopt);
        auto referenceImagePoint = referenceCoordinates.ref_sphere->image_pt.project_xy();

        std::vector<ContrastTraceByWvl<T>> result;
        rayoptics::raytr::TraceRingsDef definition;
        definition.num_rings = num_rings;
        for (double wavelength : wavelengthList) {
            auto coordinates = rayoptics::raytr::Trace::setup_pupil_coords(
                opt_model, field, wavelength, focus, referenceImagePoint, std::nullopt);
            field.chief_ray =
                std::const_pointer_cast<rayoptics::raytr::ChiefRayPkg>(coordinates.chief_ray_pkg);
            field.ref_sphere =
                std::const_pointer_cast<rayoptics::raytr::ReferenceSphere>(coordinates.ref_sphere);
            std::optional<mathlib::Vector2> sagittalExitShift;
            std::optional<mathlib::Vector2> tangentialExitShift;
            if (aim_exit_pupil) {
                // Convert image-space frequency to a physical separation on this field's
                // wavelength-specific reference sphere: Delta = |R| lambda |nu| / |n_img|,
                // with lambda in system length units and nu in reciprocal length units.
                // The vectors below are target offsets from each reference ray's sphere
                // coordinate. trace_triplets() inverse-aims the partners to those targets;
                // the supplied normalized entrance-pupil shifts provide the initial guesses.
                double physicalShift = ExitPupilAiming::referenceSphereShift(
                    opt_model, field, wavelength, spatial_frequency);
                sagittalExitShift = mathlib::Vector2(physicalShift, 0.0);
                tangentialExitShift = mathlib::Vector2(0.0, physicalShift);
            }
            auto rays = trace_triplets(
                opt_model, definition, num_spokes, sagittal_shift, tangential_shift,
                sagittalExitShift, tangentialExitShift, field, wavelength, trace_options,
                aim_exit_pupil);
            std::vector<T> samples;
            samples.reserve(rays.size());
            for (const auto &ray : rays)
                samples.push_back(callback(ray, field, wavelength, focus));
            result.push_back(ContrastTraceByWvl<T>(wavelength, std::move(samples)));
        }
        return result;
    }

    /**
     * Trace the reference and the two displaced rays required by contrast
     * optimization. Pupil displacements are expressed in relative pupil
     * coordinates of the physically vignetted pupil. Every quadrature sample
     * is returned so clients can keep a stable residual layout; a ray is null
     * when the trace fails.
     */
    static std::vector<ContrastRayTriplet> trace_triplets(
        rayoptics::optical::OpticalModel *opt_model, const rayoptics::raytr::TraceRingsDef &grid_rng,
        std::optional<int> num_spokes, const mathlib::Vector2 &sagittal_shift,
        const mathlib::Vector2 &tangential_shift, rayoptics::specs::Field &fld, double wvl,
        const rayoptics::raytr::TraceOptions &trace_options);

    /**
     * Trace contrast triplets, optionally inverse-aiming each displaced partner at its
     * requested coordinate separation on the exit-pupil reference sphere.
     */
    static std::vector<ContrastRayTriplet> trace_triplets(
        rayoptics::optical::OpticalModel *opt_model, const rayoptics::raytr::TraceRingsDef &grid_rng,
        std::optional<int> num_spokes, const mathlib::Vector2 &sagittal_shift,
        const mathlib::Vector2 &tangential_shift,
        std::optional<mathlib::Vector2> sagittal_exit_shift,
        std::optional<mathlib::Vector2> tangential_exit_shift, rayoptics::specs::Field &fld,
        double wvl, const rayoptics::raytr::TraceOptions &trace_options, bool aim_exit_pupil);

    /**
     * Generate fixed-count quadrature samples in the physical vignetted pupil.
     * The complete pattern is contracted about the overlap centre until every
     * sample and both of its displaced partners are valid. A common contraction
     * preserves relative quadrature weights and avoids frequency-dependent ray
     * loss at the pupil boundary.
     */
    static std::vector<rayoptics::raytr::GaussianQuadraturePoint> generate_quadrature(
        const rayoptics::raytr::TraceRingsDef &grid_rng, std::optional<int> num_spokes,
        const mathlib::Vector2 &sagittal_shift,
        const mathlib::Vector2 &tangential_shift, rayoptics::specs::Field &fld);

    static bool inside_vignetted_pupil(const mathlib::Vector2 &pupil,
                                       const rayoptics::specs::Field &fld);

    static ContrastAnalysisResult eval(rayoptics::optical::OpticalModel *opticalModel,
                                       const ContrastOptions &options);

    static double normalized_entry_pupil_shift(rayoptics::optical::OpticalModel *opticalModel,
                                               double wavelength,
                                               double spatialFrequency);

    static double exit_pupil_frequency_calibration(rayoptics::optical::OpticalModel *opticalModel,
                                                   rayoptics::specs::Field &field, double wavelength,
                                                   double shift, int axis,
                                                   const ContrastOptions &options);

    /**
     * Remove the constant part of the wavefront difference from every residual.
     *
     * |OTF| = |<exp(i.Phi)>| ~ 1 - Var(Phi)/2 depends on the <em>variance</em>
     * of the phase difference. A constant dW across the pupil is wavefront tilt,
     * which is an image displacement: it moves the phase transfer function and leaves the
     * modulus alone. Left in, it contributes mean^2 to the sum of squares - up to
     * 57% of an outer-field tangential block - and, being reducible by adding tilt, offers
     * the solver merit reduction that corresponds to no optical improvement at all.
     *
     * The mean is subtracted per (field, orientation), <em>not</em> per wavelength. A
     * tilt common to every wavelength is a harmless image shift, but one that differs
     * between wavelengths is lateral colour, and that genuinely does reduce polychromatic
     * MTF because the per-wavelength complex OTFs acquire different phases and partly
     * cancel. Subtracting the reference wavelength's mean from all of them discards the
     * common part and preserves the difference.
     *
     * Defocus is untouched: it makes dW linear in the shear direction rather
     * than constant, so its mean over a symmetric pupil is already zero. The sagittal mean
     * is identically zero by symmetry on a rotationally symmetric system, so in practice
     * only tangential residuals move.
     */
    static void center_residuals(ContrastAnalysisResult &result,
                                 int referenceWavelengthIndex);

private:
    /** Quadrature-weighted mean difference over the valid samples of one block. */
    static double weighted_mean(const ContrastAnalysisResult::WavelengthResult &wavelength,
                                int orientation);

    /** Direction cosine of the traced ray in image space, or null if it did not get there. */
    static std::optional<double> image_direction(rayoptics::optical::OpticalModel *opticalModel,
                                                 rayoptics::specs::Field &field, double wavelength,
                                                 const rayoptics::raytr::TraceOptions &traceOptions,
                                                 int axis, const mathlib::Vector2 &pupil);

    static ContrastAnalysisResult::Sample sample(rayoptics::optical::OpticalModel *opticalModel,
                                                 const ContrastRayTriplet &rays,
                                                 rayoptics::specs::Field &field, double wavelength,
                                                 double focus,
                                                 const ContrastOptions &options);

    static ContrastAnalysisResult::Failure failure(const ContrastRayTriplet &rays);

    static ContrastAnalysisResult::Failure failure(
        const std::string &ray, const rayoptics::exceptions::TraceException &error);

    static double opd(rayoptics::optical::OpticalModel *opticalModel,
                      const std::shared_ptr<const rayoptics::raytr::RayPkg> &ray,
                      rayoptics::specs::Field &field, double wavelength, double focus,
                      const std::optional<mathlib::Vector2> &pupil);
};

} // namespace redukti::analysis

#endif // REDUKTI_ANALYSIS_CONTRASTANALYSIS_H
