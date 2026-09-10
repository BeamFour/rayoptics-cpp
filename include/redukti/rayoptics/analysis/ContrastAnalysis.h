// C++ port of org.redukti.rayoptics.analysis.ContrastOptions,
// ContrastAnalysisResult and ContrastAnalysis.
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_ANALYSIS_CONTRASTANALYSIS_H
#define REDUKTI_RAYOPTICS_ANALYSIS_CONTRASTANALYSIS_H

#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/raytr/RayTrace.h"
#include "redukti/rayoptics/raytr/RayTypes.h"

#include <optional>
#include <string>
#include <vector>

namespace redukti::rayoptics::analysis {

/** Options for pupil-autocorrelation contrast analysis. */
/** Fluent configuration for a contrast (through-frequency) analysis. */
class ContrastOptions {
public:
    double spatialFrequency;
    int numRings = 3;
    std::optional<int> numSpokes = 6;
    raytr::TraceOptions traceOptions;
    bool calibrateFrequency = false;
    bool aimExitPupil = false;
    bool centerResiduals = false;

    explicit ContrastOptions(double spatialFrequency_);

    ContrastOptions &num_rings(int value);
    ContrastOptions &num_spokes(std::optional<int> value);
    ContrastOptions &trace_options(const raytr::TraceOptions &value);
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
        std::shared_ptr<const specs::FieldSnapshot> field;
        std::vector<WavelengthResult> wavelengths;

        FieldResult(specs::Field *field_, std::vector<WavelengthResult> wavelengths_)
            : field(field_ ? std::make_shared<const specs::FieldSnapshot>(*field_) : nullptr),
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
    static ContrastAnalysisResult eval(optical::OpticalModel *opticalModel,
                                       const ContrastOptions &options);

    static double normalized_entry_pupil_shift(optical::OpticalModel *opticalModel,
                                               double wavelength,
                                               double spatialFrequency);

    static double exit_pupil_frequency_calibration(optical::OpticalModel *opticalModel,
                                                   specs::Field &field, double wavelength,
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
    static std::optional<double> image_direction(optical::OpticalModel *opticalModel,
                                                 specs::Field &field, double wavelength,
                                                 const raytr::TraceOptions &traceOptions,
                                                 int axis, const mathlib::Vector2 &pupil);

    static ContrastAnalysisResult::Sample sample(optical::OpticalModel *opticalModel,
                                                 const raytr::ContrastRayTriplet &rays,
                                                 specs::Field &field, double wavelength,
                                                 double focus,
                                                 const ContrastOptions &options);

    static ContrastAnalysisResult::Failure failure(const raytr::ContrastRayTriplet &rays);

    static ContrastAnalysisResult::Failure failure(
        const std::string &ray, const exceptions::TraceException &error);

    static double opd(optical::OpticalModel *opticalModel,
                      const std::shared_ptr<const raytr::RayPkg> &ray,
                      specs::Field &field, double wavelength, double focus,
                      const std::optional<mathlib::Vector2> &pupil);
};

} // namespace redukti::rayoptics::analysis

#endif // REDUKTI_RAYOPTICS_ANALYSIS_CONTRASTANALYSIS_H
