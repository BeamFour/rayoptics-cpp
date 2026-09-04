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

    double spatialFrequency;
    std::vector<FieldResult> fields;

    explicit ContrastAnalysisResult(double spatialFrequency_)
        : spatialFrequency(spatialFrequency_) {}
};

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

    static void center_residuals(ContrastAnalysisResult &result,
                                 int referenceWavelengthIndex);

private:
    static double weighted_mean(const ContrastAnalysisResult::WavelengthResult &wavelength,
                                int orientation);

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
