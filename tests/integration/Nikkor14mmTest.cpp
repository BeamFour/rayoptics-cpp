// C++ port of org.redukti.rayoptics.integration.Nikkor14mmTest.
//
// A wide-angle zoom (JP2019-008031 Example 1) with four even-aspheric
// surfaces. Checks first-order data against the published values, then the
// wide-angle chief ray, reference sphere and exit-pupil aiming.
//
// The prescription was transcribed by scratchpad/xnikkor.py; the assertions
// are written out by hand because the Java writes them by hand too.
#include "../TestHarness.h"

#include "redukti/rayoptics/analysis/ContrastAnalysis.h"
#include "redukti/rayoptics/elem/profiles/EvenPolynomial.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/parax/ParaxTypes.h"
#include "redukti/rayoptics/raytr/ExitPupilAiming.h"
#include "redukti/rayoptics/raytr/RayTrace.h"
#include "redukti/rayoptics/raytr/Trace.h"
#include "redukti/rayoptics/raytr/VigCalc.h"
#include "redukti/rayoptics/raytr/Wideangle.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/rayoptics/util/Lists.h"

#include <array>
#include <cmath>
#include <memory>
#include <vector>

namespace {

using namespace redukti::rayoptics;
using redukti::mathlib::Vector2;
using redukti::mathlib::Vector3;
using util::Pair;
namespace profiles = elem::profiles;

using Separations = std::array<double, 4>;

Separations exitPupilSeparations(const raytr::ContrastRayTriplet &rays,
                                 specs::Field &field, double wavelength, double focus) {
    (void)wavelength;
    (void)focus;
    CHECK(rays.referenceError == nullptr);
    CHECK(rays.sagittalError == nullptr);
    CHECK(rays.tangentialError == nullptr);
    auto reference = raytr::ExitPupilAiming::sphere_coord(rays.reference, field.chief_ray,
                                                          field.ref_sphere);
    auto sagittal = raytr::ExitPupilAiming::sphere_coord(rays.sagittal, field.chief_ray,
                                                         field.ref_sphere);
    auto tangential = raytr::ExitPupilAiming::sphere_coord(
        rays.tangential, field.chief_ray, field.ref_sphere);
    CHECK(reference.has_value());
    CHECK(sagittal.has_value());
    CHECK(tangential.has_value());
    return Separations{sagittal->x - reference->x, sagittal->y - reference->y,
                       tangential->x - reference->x, tangential->y - reference->y};
}

double maxSeparationError(const std::vector<Separations> &separations, double requested) {
    double maximum = 0.0;
    for (const auto &s : separations) {
        maximum = std::fmax(maximum, std::hypot(s[0] - requested, s[1]));
        maximum = std::fmax(maximum, std::hypot(s[2], s[3] - requested));
    }
    return maximum;
}

/** The Java compare(RaySeg, RaySeg). */
bool compare(const raytr::RaySeg &s1, const raytr::RaySeg &s2) {
    return s1.p.effectivelyEqual(s2.p) && s1.d.effectivelyEqual(s2.d) &&
           std::abs(s1.dst - s2.dst) < 1e-13;
}

std::unique_ptr<optical::OpticalModel> build() {
    auto opm = std::make_unique<optical::OpticalModel>();
    auto sm = opm->seq_model.get();
    auto osp = opm->optical_spec.get();
    osp->pupil = std::make_unique<specs::PupilSpec>(
        osp, Pair<specs::ImageKey, specs::ValueKey>(specs::ImageKey::Image,
                                                    specs::ValueKey::Fnum),
        4.0);
    osp->fov = std::make_unique<specs::FieldSpec>(
        osp, Pair<specs::ImageKey, specs::ValueKey>(specs::ImageKey::Object,
                                                    specs::ValueKey::Angle),
        std::vector<double>{0., 57.68}, true);
    osp->wvls = std::make_unique<specs::WvlSpec>(
        std::vector<specs::WvlWt>{specs::WvlWt(486.1327, 0.5), specs::WvlWt(587.5618, 1.0), specs::WvlWt(656.2725, 0.5)}, 1);
    opm->system_spec->title = "JP2019-008031 Example 1 (Nikon Nikkor Z 14-30mm f/4 S)";
    opm->system_spec->dimensions = "MM";
    opm->radius_mode = true;
    sm->gaps[0]->thi = 1e10;
    {
        seq::SurfaceData sd(190.7535, 3.0);
        sd.rindex(1.6937, 53.32)->max_aperture(29.285);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(18.8098, 9.5);
        sd.max_aperture(22.485);
        sm->add_surface(sd);
    }
    {
        auto prof = std::make_shared<profiles::EvenPolynomial>();
        prof->r(18.8098)->setCc(-1.0)->setCoefs(std::vector<double>{0.0, -1.33157E-5, -3.07345E-8, 6.9126E-11, -3.76684E-14, 0.0, 0.0});
        sm->ifcs[static_cast<std::size_t>(*sm->cur_surface)]->profile = prof;
    }
    {
        seq::SurfaceData sd(51.563, 2.9);
        sd.rindex(1.6937, 53.32)->max_aperture(19.205);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(22.702, 9.7);
        sd.max_aperture(14.475);
        sm->add_surface(sd);
    }
    {
        auto prof = std::make_shared<profiles::EvenPolynomial>();
        prof->r(22.702)->setCc(-1.0)->setCoefs(std::vector<double>{0.0, 3.67009E-5, 1.37031E-7, -5.20756E-10, 3.14884E-12, -5.6153E-15, 0.0});
        sm->ifcs[static_cast<std::size_t>(*sm->cur_surface)]->profile = prof;
    }
    {
        seq::SurfaceData sd(-71.0651, 1.9);
        sd.rindex(1.49782, 82.57)->max_aperture(15.05);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(44.4835, 0.1);
        sd.max_aperture(15.05);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(32.608, 4.5);
        sd.rindex(1.90265, 35.73)->max_aperture(15.05);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(296.5863, 28.616);
        sd.max_aperture(15.05);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(63.0604, 2.0);
        sd.rindex(1.59349, 67.0)->max_aperture(9.04);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(499.8755, 0.1);
        sd.max_aperture(9.04);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(24.0057, 1.2);
        sd.rindex(1.883, 40.66)->max_aperture(9.605);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(13.347, 4.5);
        sd.rindex(1.56883, 56.0)->max_aperture(8.54);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(333.9818, 2.5);
        sd.max_aperture(8.54);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(0.0, 7.483);
        sd.max_aperture(5.6335);
        sm->add_surface(sd);
    }
    sm->set_stop();
    {
        seq::SurfaceData sd(36.3784, 1.1);
        sd.rindex(1.816, 46.59)->max_aperture(8.59);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(14.0097, 4.71);
        sd.rindex(1.51612, 64.08)->max_aperture(8.42);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(61.0448, 0.2);
        sd.max_aperture(8.42);
        sm->add_surface(sd);
    }
    {
        auto prof = std::make_shared<profiles::EvenPolynomial>();
        prof->r(61.0448)->setCc(0)->setCoefs(std::vector<double>{0.0, 1.75905E-5, -6.64635E-8, 2.26551E-10, -4.40763E-12, 0.0, 0.0});
        sm->ifcs[static_cast<std::size_t>(*sm->cur_surface)]->profile = prof;
    }
    {
        seq::SurfaceData sd(27.9719, 3.15);
        sd.rindex(1.49782, 82.57)->max_aperture(8.55);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(-75.3921, 0.25);
        sd.max_aperture(8.55);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(91.9654, 3.05);
        sd.rindex(1.49782, 82.57)->max_aperture(8.915);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(-29.3923, 1.579);
        sd.max_aperture(8.915);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(72.093, 1.0);
        sd.rindex(1.795, 45.31)->max_aperture(9.065);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(20.9929, 5.766);
        sd.max_aperture(9.065);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(-538.2301, 4.8);
        sd.rindex(1.49782, 82.57)->max_aperture(10.935);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(-20.1257, 0.1);
        sd.max_aperture(10.935);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(-38.9341, 1.4);
        sd.rindex(1.76546, 46.75)->max_aperture(11.06);
        sm->add_surface(sd);
    }
    {
        auto prof = std::make_shared<profiles::EvenPolynomial>();
        prof->r(-38.9341)->setCc(-1.0)->setCoefs(std::vector<double>{0.0, -2.67902E-5, -3.34364E-8, -1.13765E-10, -1.88017E-13, 0.0, 0.0});
        sm->ifcs[static_cast<std::size_t>(*sm->cur_surface)]->profile = prof;
    }
    {
        seq::SurfaceData sd(154.832, 21.36);
        sd.max_aperture(11.815);
        sm->add_surface(sd);
    }
    sm->do_apertures = false;
    opm->update_model();
    opm->update_model();
    sm->do_apertures = false;
    opm->update_model();
    return opm;
}

} // namespace

TEST(integration_nikkor14mm_first_order_and_wide_angle_aiming) {
    auto opm = build();
    auto sm = opm->seq_model.get();
    auto osp = opm->optical_spec.get();
    const auto &fod = osp->parax_data->fod;

    CHECK_CLOSE(fod.efl, 14.42, 0.001);
    CHECK_CLOSE(fod.ffl, 15.5, 0.01);
    CHECK_CLOSE(fod.pp1, 29.92, 0.01);
    CHECK_CLOSE(fod.ppk, 6.939, 0.01);
    CHECK_CLOSE(fod.bfl, 21.36, 0.001);
    CHECK_CLOSE(fod.fno, 4.0, 0.001);
    CHECK_CLOSE(fod.m, -1.442e-09, 1e-9);
    CHECK_CLOSE(fod.obj_ang, 57.68, 0.01);
    CHECK_CLOSE(fod.enp_dist, 19.84, 0.01);
    CHECK_CLOSE(fod.enp_radius, 1.803, 0.01);
    CHECK_CLOSE(fod.obj_na, 1.803e-10, 1e-9);
    CHECK_CLOSE(fod.n_obj, 1.000, 1e-9);
    CHECK_CLOSE(fod.img_dist, 21.36, 0.01);
    CHECK_CLOSE(fod.img_ht, 22.79, 0.01);
    CHECK_CLOSE(fod.exp_dist, -26.59, 0.01);
    CHECK_CLOSE(fod.exp_radius, 5.993, 0.01);
    CHECK_CLOSE(fod.n_img, 1.000, 0.001);
    CHECK_CLOSE(fod.opt_inv, 2.849, 0.001);

    raytr::VigCalc::set_vig(opm.get());
    opm->update_model();

    const double cr_expected_op_delta[] = {129.68211720000002, 144.5443500746337};
    const double expected_z_enp[] = {19.84067682425449, 19.127560250656217};
    const double expected_ref_sphere_radius[] = {47.948425307468156, 45.491290263961105};
    const raytr::RaySeg cr_expected_final_rayseg[] = {
        raytr::RaySeg(Vector3(0., 0., 0.), Vector3(0., 0., 1.), 0.0,
                      Vector3(-0., -0., 1.)),
        raytr::RaySeg(Vector3(0., 20.362505270545146, 0.),
                      Vector3(0., 0.4476133316595901, 0.8942272112391799), 0.0,
                      Vector3(-0., -0., 1.))};

    for (std::size_t fi = 1; fi < osp->fov->fields.size(); fi++) {
        specs::Field &fld = *osp->fov->fields[fi];
        auto wvl = sm->central_wavelength();
        auto foc = osp->defocus()->get_focus();
        auto t = raytr::Trace::setup_pupil_coords(opm.get(), fld, wvl, foc, std::nullopt,
                                                  std::nullopt);
        CHECK(fld.z_enp.has_value());
        CHECK_CLOSE(*fld.z_enp, expected_z_enp[fi], 1e-5);
        CHECK_CLOSE(t.ref_sphere->ref_sphere_radius, expected_ref_sphere_radius[fi], 1e-5);
        CHECK_CLOSE(t.chief_ray_pkg->chief_ray->op_delta, cr_expected_op_delta[fi], 1e-5);
        CHECK(compare(util::Lists::get(t.chief_ray_pkg->chief_ray->ray, -1),
                      cr_expected_final_rayseg[fi]));
        fld.chief_ray = std::const_pointer_cast<raytr::ChiefRayPkg>(t.chief_ray_pkg);
        fld.ref_sphere = std::const_pointer_cast<raytr::ReferenceSphere>(t.ref_sphere);
    }

    int fieldIndex = 1;
    int wavelengthIndex = osp->wvls->reference_wvl;
    double wavelength = osp->wvls->wavelengths[static_cast<std::size_t>(wavelengthIndex)];
    double normalizedShift = analysis::ContrastAnalysis::normalized_entry_pupil_shift(
        opm.get(), wavelength, 40.0);
    Vector2 sagittalShift(normalizedShift, 0.0);
    Vector2 tangentialShift(0.0, normalizedShift);

    raytr::ContrastTraceCallback<Separations> sep = exitPupilSeparations;
    raytr::TraceOptions options;
    auto unaimed = sm->trace_contrast<Separations>(sep, fieldIndex, wavelengthIndex, 1, 6,
                                                   sagittalShift, tangentialShift, 40.0,
                                                   options, false);
    auto aimed = sm->trace_contrast<Separations>(sep, fieldIndex, wavelengthIndex, 1, 6,
                                                 sagittalShift, tangentialShift, 40.0,
                                                 options, true);
    double requestedShift = raytr::ExitPupilAiming::referenceSphereShift(
        opm.get(), *osp->fov->fields[static_cast<std::size_t>(fieldIndex)], wavelength,
        40.0);
    double unaimedError = maxSeparationError(unaimed[0].samples, requestedShift);
    double aimedError = maxSeparationError(aimed[0].samples, requestedShift);
    // wide-angle unaimed shear should expose pupil mapping error
    CHECK(unaimedError > 0.5);
    // aimed shear should reach the reference-sphere target
    CHECK(aimedError < 2.0e-6);
    // aiming should improve wide-angle shear by at least five orders of magnitude
    CHECK(aimedError * 100000.0 < unaimedError);

    auto result =
        raytr::Wideangle::eval_real_image_ht(opm.get(), *osp->fov->fields[1], 587.5618);
    Vector3 expect_pt(0.0, -2866312975.4227800369262695, 419590299.1519107818603516);
    Vector3 expect_dir(-0., 0.2866312938130761, 0.9580409706307147);
    double expect_z_enp = 130.10449270101637;
    CHECK_CLOSE(result.z_enp, expect_z_enp, 1e-5);
    CHECK(expect_pt.isEqual(result.ray_data.pt, 1e-7));
    CHECK(expect_dir.isEqual(result.ray_data.dir, 1e-7));
}
