// C++ port of org.redukti.rayoptics.analysis.WeightedSpotAnalysisTest.
//
// The per-ray weights the Gaussian-quadrature pattern attaches to its grid
// items have to reach the centroid and the RMS radius. Two hand-built rays with
// known weights make that arithmetic checkable without tracing anything.
#include "TestHarness.h"

#include "redukti/mathlib/Matrix3.h"
#include "redukti/mathlib/Vector2.h"
#include "redukti/mathlib/Vector3.h"
#include "redukti/rayoptics/analysis/SpotAnalysis.h"
#include "redukti/rayoptics/analysis/SpotIntercepts.h"
#include "redukti/rayoptics/math/Tfm3d.h"
#include "redukti/rayoptics/raytr/RayTypes.h"
#include "redukti/rayoptics/specs/Field.h"

#include <cmath>
#include <memory>
#include <vector>

namespace {

using redukti::mathlib::Matrix3;
using redukti::mathlib::Vector2;
using redukti::mathlib::Vector3;
using redukti::rayoptics::analysis::SpotAnalysisResult;
using redukti::rayoptics::analysis::SpotIntercepts;
using redukti::rayoptics::analysis::SpotOptions;
using redukti::rayoptics::math::Tfm3d;
using redukti::rayoptics::raytr::GridItem;
using redukti::rayoptics::raytr::ReferenceSphere;
using redukti::rayoptics::raytr::TraceGridByWvl;
using redukti::rayoptics::specs::Field;

} // namespace

TEST(spot_computes_weighted_centroid_and_rms_radius) {
    std::vector<GridItem> items{GridItem(Vector2(0.0, 0.0), nullptr).withWeight(0.75),
                                GridItem(Vector2(2.0, 0.0), nullptr).withWeight(0.25)};
    TraceGridByWvl trace(550.0, items);
    SpotIntercepts intercepts(trace);

    CHECK_CLOSE(intercepts.compute_centroid().x, 0.5, 1.0e-15);

    Field field(nullptr);
    // The Java passes null for the last-surface transform; it is never read on
    // this path, and a Tfm3d is a value here, so an identity stands in.
    field.ref_sphere = std::make_shared<ReferenceSphere>(
        Vector3::ZERO, Vector3::ZERO, 1.0, Tfm3d(Matrix3::IDENTITY, Vector3::ZERO));
    SpotAnalysisResult::SpotResultsForField result(
        &field, std::vector<TraceGridByWvl>{trace}, 550.0, true);

    CHECK_CLOSE(result.get_mean_radius(), std::sqrt(0.75) * 1000.0, 1.0e-12);
}

TEST(spot_selecting_grid_and_quadrature_patterns_is_unambiguous) {
    SpotOptions options(true);
    CHECK(options.is_gauss_quadrature());
    CHECK(!options.is_grid());
    CHECK(!options.is_hexapolar());

    options.use_grid();
    CHECK(options.is_grid());
    CHECK(!options.is_gauss_quadrature());
    CHECK(!options.is_hexapolar());
}
