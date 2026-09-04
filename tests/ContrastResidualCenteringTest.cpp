// C++ port of org.redukti.rayoptics.analysis.ContrastResidualCenteringTest.
//
// Centring must remove exactly the constant part of a block and nothing else:
// sum(w.dW^2) = Var + mean^2 becomes Var.
//
// Built from synthetic samples rather than a traced lens so the mean is known
// exactly and the identity can be asserted rather than approximated.
#include "TestHarness.h"

#include "redukti/mathlib/Vector2.h"
#include "redukti/rayoptics/analysis/ContrastAnalysis.h"

#include <cmath>
#include <vector>

namespace {

using redukti::mathlib::Vector2;
using redukti::rayoptics::analysis::ContrastAnalysis;
using redukti::rayoptics::analysis::ContrastAnalysisResult;

// Weights sum to 1, as the quadrature normalizes them, and the shape is chosen
// with weighted mean exactly zero so the block's mean is the injected tilt and
// nothing else. WEIGHTS and SHAPE are checked against that premise below.
const std::vector<double> WEIGHTS{0.1, 0.2, 0.4, 0.2, 0.1};
const std::vector<double> SHAPE{-0.315, -0.115, 0.035, 0.085, 0.235};

/** A block whose tangential differences carry a known tilt on top of a known shape. */
ContrastAnalysisResult blockWithTangentialTilt(double tilt) {
    std::vector<ContrastAnalysisResult::Sample> samples;
    for (std::size_t i = 0; i < WEIGHTS.size(); i++) {
        samples.push_back(ContrastAnalysisResult::Sample(
            Vector2(0.1 * static_cast<double>(i), 0.0),
            // Sagittal deliberately carries no tilt: on a rotationally
            // symmetric system its mean is identically zero.
            SHAPE[i], SHAPE[i] + tilt, WEIGHTS[i], true, std::nullopt));
    }
    ContrastAnalysisResult result(40.0);
    std::vector<ContrastAnalysisResult::WavelengthResult> wavelengths;
    wavelengths.push_back(
        ContrastAnalysisResult::WavelengthResult(587.5618, 0.0674, samples));
    result.fields.push_back(ContrastAnalysisResult::FieldResult(nullptr, wavelengths));
    return result;
}

double sumOfSquares(const ContrastAnalysisResult &result, bool tangential) {
    const auto &block = result.fields[0].wavelengths[0];
    double sos = 0.0;
    for (std::size_t i = 0; i < block.samples.size(); i++) {
        double r = tangential ? block.tangentialResidual(static_cast<int>(i))
                              : block.sagittalResidual(static_cast<int>(i));
        sos += r * r;
    }
    return sos;
}

/** One unit-weight sample whose tangential difference is exactly `value`. */
ContrastAnalysisResult::WavelengthResult block(double wavelength, double value) {
    std::vector<ContrastAnalysisResult::Sample> samples{ContrastAnalysisResult::Sample(
        Vector2(0.0, 0.0), 0.0, value, 1.0, true, std::nullopt)};
    return ContrastAnalysisResult::WavelengthResult(wavelength, 0.0674, samples);
}

} // namespace

TEST(contrast_centering_fixture_shape_has_zero_weighted_mean) {
    double mean = 0.0;
    for (std::size_t i = 0; i < WEIGHTS.size(); i++)
        mean += WEIGHTS[i] * SHAPE[i];
    // the other tests read the block mean as the injected tilt alone
    CHECK_CLOSE(mean, 0.0, 1.0e-15);
}

TEST(contrast_centering_removes_exactly_the_mean_square_term) {
    double tilt = 0.22;
    auto uncentred = blockWithTangentialTilt(tilt);
    double before = sumOfSquares(uncentred, true);

    auto centred = blockWithTangentialTilt(tilt);
    ContrastAnalysis::center_residuals(centred, 0);
    double after = sumOfSquares(centred, true);

    // The shape has weighted mean zero, so the block's mean is exactly the tilt
    // and the removed term is exactly mean^2 (weights sum to 1).
    const auto &blk = centred.fields[0].wavelengths[0];
    CHECK_CLOSE(blk.tangentialOffset, tilt, 1.0e-12);
    CHECK_CLOSE(after, before - tilt * tilt, 1.0e-12);
    CHECK(before != after);
}

TEST(contrast_centred_residuals_have_zero_weighted_mean) {
    auto centred = blockWithTangentialTilt(0.22);
    ContrastAnalysis::center_residuals(centred, 0);
    const auto &blk = centred.fields[0].wavelengths[0];

    double weightedMean = 0.0;
    for (std::size_t i = 0; i < blk.samples.size(); i++)
        weightedMean += std::sqrt(blk.samples[i].weight) *
                        blk.tangentialResidual(static_cast<int>(i));
    CHECK_CLOSE(weightedMean, 0.0, 1.0e-12);
}

/** A block with no tilt must come through untouched, as sagittal blocks always do. */
TEST(contrast_centering_leaves_a_block_with_no_constant_part_unchanged) {
    auto result = blockWithTangentialTilt(0.0);
    double before = sumOfSquares(result, false);
    ContrastAnalysis::center_residuals(result, 0);
    CHECK_CLOSE(result.fields[0].wavelengths[0].sagittalOffset, 0.0, 1.0e-12);
    CHECK_CLOSE(sumOfSquares(result, false), before, 1.0e-12);
}

/**
 * Every wavelength gets the reference wavelength's offset, not its own: a tilt
 * common to all wavelengths is a harmless image shift, while one that differs
 * between them is lateral colour and does cost polychromatic MTF. Subtracting
 * per-wavelength means would discard that difference.
 */
TEST(contrast_all_wavelengths_share_the_reference_wavelength_offset) {
    ContrastAnalysisResult result(40.0);
    std::vector<ContrastAnalysisResult::WavelengthResult> wavelengths;
    wavelengths.push_back(block(486.13, 0.10));
    wavelengths.push_back(block(587.56, 0.22)); // reference
    wavelengths.push_back(block(656.27, 0.31));
    result.fields.push_back(ContrastAnalysisResult::FieldResult(nullptr, wavelengths));

    ContrastAnalysis::center_residuals(result, 1);

    const auto &centred = result.fields[0].wavelengths;
    for (const auto &blk : centred)
        CHECK_CLOSE(blk.tangentialOffset, 0.22, 1.0e-12);

    // The colour differences survive: F and C keep their offset from the reference.
    CHECK_CLOSE(centred[0].tangentialResidual(0), 0.10 - 0.22, 1.0e-12);
    CHECK_CLOSE(centred[2].tangentialResidual(0), 0.31 - 0.22, 1.0e-12);
}
