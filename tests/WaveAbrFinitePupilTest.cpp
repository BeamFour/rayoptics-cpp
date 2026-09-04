// C++ port of org.redukti.rayoptics.raytr.WaveAbrFinitePupilTest.
//
// Regression tests for the finite-reference-sphere correction in WaveAbr, and
// in particular the sign of the Hopkins discriminant: upstream ray-optics
// carried F^2 + J/R until commit 2de1e18 ("Fix sign error in OPD closing
// equation. Fixes #221"), and the test below shows the two are not equivalent
// conventions -- the '+' form does not put the chord point on the sphere.
#include "TestHarness.h"

#include "integration/IntegrationModels.h"

#include "redukti/mathlib/Vector3.h"
#include "redukti/rayoptics/elem/transform/Transform.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/parax/ParaxTypes.h"
#include "redukti/rayoptics/raytr/RayTypes.h"
#include "redukti/rayoptics/raytr/Trace.h"
#include "redukti/rayoptics/raytr/WaveAbr.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/specs/FieldSpec.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/rayoptics/util/Lists.h"

#include <cmath>
#include <memory>
#include <vector>

namespace {

using redukti::mathlib::Vector3;
using redukti::rayoptics::elem::transform::Transform;
using redukti::rayoptics::optical::OpticalModel;
using redukti::rayoptics::raytr::ChiefRayPkg;
using redukti::rayoptics::raytr::RayData;
using redukti::rayoptics::raytr::RayPkg;
using redukti::rayoptics::raytr::ReferenceSphere;
using redukti::rayoptics::raytr::RefSphereCR;
using redukti::rayoptics::raytr::Trace;
using redukti::rayoptics::raytr::TraceOptions;
using redukti::rayoptics::raytr::WaveAbr;
using redukti::rayoptics::specs::Field;
using redukti::rayoptics::util::Lists::get;

/** Java's private record SphereGeometry. */
struct SphereGeometry {
    double radius;
    double f;
    double j;
    double solutionSign;
    double ekp;
};

/** Java's private record Fixture; the model owns the field, so it is kept alive. */
struct Fixture {
    std::unique_ptr<OpticalModel> model;
    Field *field;
    double wavelength;
    double focus;
    RefSphereCR setup;
    std::shared_ptr<const RayPkg> ray;
};

Fixture tracedOffAxisRay() {
    auto model = integration::build_US003549241Example05();
    Field *field = model->optical_spec->fov->fields[2].get();
    double wavelength = model->seq_model->central_wavelength();
    double focus = model->optical_spec->defocus()->get_focus();
    auto setup = Trace::setup_pupil_coords(model.get(), *field, wavelength, focus,
                                           std::nullopt, std::nullopt);
    field->chief_ray =
        std::const_pointer_cast<ChiefRayPkg>(setup.chief_ray_pkg);
    field->ref_sphere = std::const_pointer_cast<ReferenceSphere>(setup.ref_sphere);
    auto ray = Trace::trace_base(model.get(), std::vector<double>{0.65, 0.45}, *field,
                                 wavelength, TraceOptions());
    return Fixture{std::move(model), field, wavelength, focus, setup, ray};
}

SphereGeometry sphereGeometry(const std::shared_ptr<const RayPkg> &rayPkg,
                              const std::shared_ptr<const ChiefRayPkg> &chiefRayPkg,
                              const std::shared_ptr<const ReferenceSphere> &sphere) {
    int lastRealSurface = -2;
    const auto &ray = rayPkg->ray;
    const auto &chiefRay = chiefRayPkg->chief_ray->ray;
    double ekp = WaveAbr::eic_distance(
        RayData(get(ray, lastRealSurface).p, get(ray, lastRealSurface).d),
        RayData(get(chiefRay, lastRealSurface).p, get(chiefRay, lastRealSurface).d));
    auto after = Transform::transform_after_surface(
        *chiefRayPkg->cr_exp_seg->ifc,
        RayData(get(ray, lastRealSurface).p, get(ray, lastRealSurface).d));
    double distance = ekp - chiefRayPkg->cr_exp_seg->exp_dst;
    auto chordCoordinate = after.pt.minus(after.dir.times(distance))
                               .minus(chiefRayPkg->cr_exp_seg->exp_pt);
    double radius = sphere->ref_sphere_radius;
    double f = sphere->ref_dir.dot(after.dir) - after.dir.dot(chordCoordinate) / radius;
    double j = chordCoordinate.dot(chordCoordinate) / radius -
               2.0 * sphere->ref_dir.dot(chordCoordinate);
    double solutionSign = sphere->ref_dir.z * get(chiefRay, -1).d.z < 0.0 ? -1.0 : 1.0;
    return SphereGeometry{radius, f, j, solutionSign, ekp};
}

double sphereDistance(const SphereGeometry &geometry, bool upstreamPlusSign) {
    double signedTerm = geometry.j / geometry.radius;
    double discriminant =
        geometry.f * geometry.f + (upstreamPlusSign ? signedTerm : -signedTerm);
    // ray must intersect the reference sphere
    CHECK(discriminant >= 0.0);
    double denominator = geometry.f + geometry.solutionSign * std::sqrt(discriminant);
    return denominator == 0.0 ? 0.0 : geometry.j / denominator;
}

double sphereResidual(const SphereGeometry &geometry, double ep) {
    return ep * ep - 2.0 * geometry.radius * geometry.f * ep + geometry.radius * geometry.j;
}

} // namespace

/**
 * Small debugger entry point for following one off-axis sample through
 * wave_abr_full_calc_finite_pup().
 */
TEST(waveabr_calculates_finite_pupil_opd_for_a_single_off_axis_sample_ray) {
    auto fixture = tracedOffAxisRay();
    const auto &fod = fixture.model->optical_spec->parax_data->fod;

    // wave_abr_full_calc() dispatches to wave_abr_full_calc_finite_pup()
    // because this fixture has a finite reference-sphere radius.
    double opd = WaveAbr::wave_abr_full_calc(fod, *fixture.field, fixture.wavelength,
                                             fixture.focus, fixture.ray,
                                             fixture.setup.chief_ray_pkg,
                                             fixture.setup.ref_sphere);

    CHECK(std::isfinite(fixture.setup.ref_sphere->ref_sphere_radius));
    CHECK(std::isfinite(opd));
}

TEST(waveabr_hopkins_discriminant_places_chord_point_on_reference_sphere) {
    auto fixture = tracedOffAxisRay();
    auto geometry = sphereGeometry(fixture.ray, fixture.setup.chief_ray_pkg,
                                   fixture.setup.ref_sphere);

    double ep = sphereDistance(geometry, false);
    double residual = sphereResidual(geometry, ep);

    // F^2 - J/R must put the displaced chord point on the reference sphere
    CHECK_CLOSE(residual, 0.0, 1.0e-11);

    // The '+' form is what this calculation used before the sign fix, and what
    // upstream ray-optics carried until commit 2de1e18. It is kept here to
    // demonstrate that the two are not equivalent sign conventions: the
    // resulting point is not on the sphere.
    double upstreamEp = sphereDistance(geometry, true);
    double upstreamResidual = sphereResidual(geometry, upstreamEp);
    // F^2 + J/R should expose a measurable sphere-equation residual
    CHECK(std::abs(upstreamResidual) > 1.0e-6);
}

TEST(waveabr_finite_pupil_opd_uses_the_geometric_reference_sphere_intersection) {
    auto fixture = tracedOffAxisRay();
    const auto &chief = fixture.setup.chief_ray_pkg;
    const auto &sphere = fixture.setup.ref_sphere;
    const auto &fod = fixture.model->optical_spec->parax_data->fod;
    auto geometry = sphereGeometry(fixture.ray, chief, sphere);

    double e1 = WaveAbr::eic_distance(
        RayData(fixture.ray->ray[1].p, fixture.ray->ray[0].d),
        RayData(chief->chief_ray->ray[1].p, chief->chief_ray->ray[0].d));
    double ep = sphereDistance(geometry, false);
    double expected = -std::abs(fod.n_obj) * e1 - fixture.ray->op_delta +
                      std::abs(fod.n_img) * geometry.ekp + chief->chief_ray->op_delta -
                      std::abs(fod.n_img) * ep;

    double actual = WaveAbr::wave_abr_full_calc(fod, *fixture.field, fixture.wavelength,
                                                fixture.focus, fixture.ray, chief, sphere);

    // This connects the geometric invariant above to the public OPD calculation:
    // it passes only while WaveAbr terminates the ray on the reference sphere,
    // and would fail if the discriminant reverted to F^2 + J/R.
    CHECK_CLOSE(actual, expected, 1.0e-12);
}
