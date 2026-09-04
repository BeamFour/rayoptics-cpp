// C++ port of org.redukti.rayoptics.raytr.GaussianQuadraturePatternTest.
#include "TestHarness.h"

#include "redukti/mathlib/Vector2.h"
#include "redukti/rayoptics/raytr/RayTypes.h"
#include "redukti/rayoptics/raytr/Trace.h"
#include "redukti/rayoptics/specs/Field.h"

#include <cmath>
#include <limits>
#include <vector>

namespace {

using redukti::mathlib::Vector2;
using redukti::rayoptics::raytr::GaussianQuadraturePoint;
using redukti::rayoptics::raytr::Trace;
using redukti::rayoptics::raytr::TraceRingsDef;
using redukti::rayoptics::specs::Field;

double totalWeight(const std::vector<GaussianQuadraturePoint> &points) {
    double sum = 0.0;
    for (const auto &point : points)
        sum += point.weight;
    return sum;
}

double weightedMoment(const std::vector<GaussianQuadraturePoint> &points, int xPower,
                      int yPower) {
    double sum = 0.0;
    for (const auto &point : points)
        sum += point.weight * std::pow(point.pupil.x, xPower) *
               std::pow(point.pupil.y, yPower);
    return sum;
}

double vignettedY(const Field &field, const Vector2 &pupil) {
    return field.apply_vignetting(std::vector<double>{pupil.x, pupil.y})[1];
}

Field asymmetricVignettedField() {
    Field field(nullptr);
    field.vux = 0.216;
    field.vlx = 0.216;
    field.vuy = 0.636;
    field.vly = 0.687;
    return field;
}

} // namespace

TEST(gq_generates_requested_number_of_rings_and_spokes) {
    TraceRingsDef definition;

    auto points = Trace::generate_gaussian_quadrature(definition, 3, 8);

    CHECK_EQ(static_cast<int>(points.size()), 24);
    CHECK_CLOSE(totalWeight(points), 1.0, 1.0e-14);
}

TEST(gq_uses_optiland_default_spoke_count) {
    TraceRingsDef definition;

    auto points = Trace::generate_gaussian_quadrature(definition, 3, std::nullopt);

    CHECK_EQ(static_cast<int>(points.size()), 3 * 4 * (3 + 1));
}

TEST(gq_integrates_low_order_disk_moments) {
    TraceRingsDef definition;
    auto points = Trace::generate_gaussian_quadrature(definition, 3, std::nullopt);

    double meanX = weightedMoment(points, 1, 0);
    double meanY = weightedMoment(points, 0, 1);
    double meanX2 = weightedMoment(points, 2, 0);
    double meanY2 = weightedMoment(points, 0, 2);
    double meanR4 = 0.0;
    for (const auto &point : points) {
        double r2 = point.pupil.x * point.pupil.x + point.pupil.y * point.pupil.y;
        meanR4 += point.weight * r2 * r2;
    }

    CHECK_CLOSE(meanX, 0.0, 1.0e-14);
    CHECK_CLOSE(meanY, 0.0, 1.0e-14);
    CHECK_CLOSE(meanX2, 0.25, 1.0e-14);
    CHECK_CLOSE(meanY2, 0.25, 1.0e-14);
    CHECK_CLOSE(meanR4, 1.0 / 3.0, 1.0e-14);
}

TEST(gq_integrates_low_order_annular_moments) {
    TraceRingsDef definition;
    definition.min_radius = 0.5;
    auto points = Trace::generate_gaussian_quadrature(definition, 3, 12);

    double meanR2 = 0.0;
    for (const auto &point : points) {
        double r2 = point.pupil.x * point.pupil.x + point.pupil.y * point.pupil.y;
        meanR2 += point.weight * r2;
    }
    CHECK_CLOSE(meanR2, (1.0 + 0.25) / 2.0, 1.0e-14);
    for (const auto &point : points)
        CHECK(point.pupil.x * point.pupil.x + point.pupil.y * point.pupil.y > 0.25);
}

TEST(gq_applies_pupil_scale_and_offset) {
    TraceRingsDef definition;
    definition.cx = 2.0;
    definition.cy = -3.0;
    definition.max_radius = 0.5;

    auto points = Trace::generate_gaussian_quadrature(definition, 2, 4);

    for (const auto &point : points) {
        double x = point.pupil.x - definition.cx;
        double y = point.pupil.y - definition.cy;
        double radius = std::sqrt(x * x + y * y);
        CHECK(radius < definition.max_radius);
    }
}

TEST(gq_rejects_invalid_dimensions) {
    using redukti::IllegalArgumentException;
    TraceRingsDef definition;

    CHECK_THROWS(Trace::generate_gaussian_quadrature(definition, 0, std::nullopt),
                 IllegalArgumentException);
    CHECK_THROWS(Trace::generate_gaussian_quadrature(definition, 1, 0),
                 IllegalArgumentException);
    CHECK_THROWS(Trace::generate_gaussian_quadrature(definition, 1, 2),
                 IllegalArgumentException);
    definition.min_radius = 1.0;
    CHECK_THROWS(Trace::generate_gaussian_quadrature(definition, 1, 3),
                 IllegalArgumentException);
}

TEST(gq_contrast_pattern_and_both_partners_stay_inside_the_vignetted_pupil) {
    TraceRingsDef definition;
    definition.num_rings = 3;
    Field field = asymmetricVignettedField();
    Vector2 sagittalShift(0.06743, 0.0);
    Vector2 tangentialShift(0.0, 0.06743);

    auto points = Trace::generate_contrast_quadrature(definition, 6, sagittalShift,
                                                      tangentialShift, field);

    CHECK_EQ(static_cast<int>(points.size()), 18);
    CHECK_CLOSE(totalWeight(points), 1.0, 1.0e-14);
    for (const auto &point : points) {
        const Vector2 &pupil = point.pupil;
        CHECK(Trace::inside_vignetted_pupil(pupil, field));
        CHECK(Trace::inside_vignetted_pupil(pupil.plus(sagittalShift), field));
        CHECK(Trace::inside_vignetted_pupil(pupil.plus(tangentialShift), field));
    }
}

/**
 * The samples are absolute pupil coordinates, so the tracer must be told not to
 * apply Field vignetting a second time. Re-applying it would rescale the
 * displaced rays and turn the requested shear into a smaller, field- and
 * direction-dependent one - the defect this pattern exists to avoid.
 */
TEST(gq_reapplying_field_vignetting_would_corrupt_the_requested_shear) {
    TraceRingsDef definition;
    definition.num_rings = 3;
    Field field = asymmetricVignettedField();
    Vector2 tangentialShift(0.0, 0.06743);

    auto points = Trace::generate_contrast_quadrature(definition, 6, Vector2(0.0, 0.0),
                                                      tangentialShift, field);

    for (const auto &point : points) {
        const Vector2 &pupil = point.pupil;
        double reapplied =
            vignettedY(field, pupil.plus(tangentialShift)) - vignettedY(field, pupil);
        // 1-vuy=0.364 and 1-vly=0.313, so a second application shrinks the
        // 40 cycle/mm shear to somewhere near a third of its intended size.
        CHECK(reapplied < 0.5 * tangentialShift.y);
    }
}

/**
 * A base sample and its displaced partner can sit on opposite sides of the y
 * axis, where the vignetting factors differ. Re-applying vignetting there would
 * not even be a rigid translation.
 */
TEST(gq_shear_straddling_the_vignetting_axis_is_not_a_rigid_translation) {
    TraceRingsDef definition;
    definition.num_rings = 3;
    Field field = asymmetricVignettedField();
    Vector2 tangentialShift(0.0, 0.06743);

    auto points = Trace::generate_contrast_quadrature(definition, 6, Vector2(0.0, 0.0),
                                                      tangentialShift, field);

    const GaussianQuadraturePoint *crossing = nullptr;
    const GaussianQuadraturePoint *below = nullptr;
    for (const auto &point : points) {
        if (crossing == nullptr && point.pupil.y < 0.0 &&
            point.pupil.y + tangentialShift.y > 0.0)
            crossing = &point;
        if (below == nullptr && point.pupil.y + tangentialShift.y < 0.0)
            below = &point;
    }
    CHECK(crossing != nullptr);
    CHECK(below != nullptr);
    if (crossing == nullptr || below == nullptr)
        return;

    double crossingShear = vignettedY(field, crossing->pupil.plus(tangentialShift)) -
                           vignettedY(field, crossing->pupil);
    double belowShear = vignettedY(field, below->pupil.plus(tangentialShift)) -
                        vignettedY(field, below->pupil);
    // a straddling pair and a wholly-below pair should disagree once vignetting
    // is re-applied
    CHECK(std::abs(crossingShear - belowShear) > 1.0e-6);
}

TEST(gq_rejects_shear_that_leaves_too_little_pupil_to_sample) {
    using redukti::IllegalArgumentException;
    TraceRingsDef definition;
    definition.num_rings = 3;
    Field field(nullptr);

    // Well inside the pupil: sampled normally.
    Trace::generate_contrast_quadrature(definition, 6, Vector2(0.5, 0.0),
                                        Vector2(0.0, 0.5), field);

    // The overlap centre is still inside the pupil here, but the pattern would
    // have to collapse onto it, reporting zero wavefront difference.
    CHECK_THROWS(Trace::generate_contrast_quadrature(definition, 6, Vector2(1.41, 0.0),
                                                     Vector2(0.0, 1.41), field),
                 IllegalArgumentException);

    // Beyond the point where even the overlap centre is valid.
    CHECK_THROWS(Trace::generate_contrast_quadrature(definition, 6, Vector2(1.9, 0.0),
                                                     Vector2(0.0, 1.9), field),
                 IllegalArgumentException);
}

TEST(gq_contrast_pattern_samples_are_distinct) {
    TraceRingsDef definition;
    definition.num_rings = 3;
    Field field = asymmetricVignettedField();
    Vector2 shift(0.06743, 0.0);

    auto points = Trace::generate_contrast_quadrature(definition, 6, shift,
                                                      Vector2(0.0, 0.06743), field);

    double closest = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < points.size(); i++)
        for (std::size_t j = i + 1; j < points.size(); j++)
            closest = std::min(closest, points[i].pupil.minus(points[j].pupil).len());
    // samples collapsed towards the overlap centre if this fails
    CHECK(closest > 0.01);
}
