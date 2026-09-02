// Expectations dumped from the Java classes in rayoptics/target/classes on
// JDK 25. Asserted exactly, as elsewhere in this port.
//
// Note on toString: Glass, Medium, EvenPolynomial and RadialPolynomial never
// override no-arg toString() in the Java -- only toString(StringBuilder) --
// so printing one directly yields an identity hash there. The port maps the
// StringBuilder formatter onto toString(), and the expectations below come
// from calling that formatter.
#include "TestHarness.h"

#include "redukti/Exceptions.h"
#include "redukti/rayoptics/elem/profiles/EvenPolynomial.h"
#include "redukti/rayoptics/elem/profiles/RadialPolynomial.h"
#include "redukti/rayoptics/elem/profiles/Spherical.h"
#include "redukti/rayoptics/elem/surface/Aperture.h"
#include "redukti/rayoptics/elem/surface/DecenterData.h"
#include "redukti/rayoptics/elem/surface/Surface.h"
#include "redukti/rayoptics/exceptions/TraceException.h"
#include "redukti/rayoptics/seq/Gap.h"
#include "redukti/rayoptics/seq/Glass.h"
#include "redukti/rayoptics/util/Lists.h"

#include <memory>
#include <vector>

using namespace redukti::rayoptics;
using redukti::IllegalArgumentException;
using redukti::mathlib::Vector2;
using redukti::mathlib::Vector3;

// ---------------------------------------------------------------------------
// Glass and the catalog
// ---------------------------------------------------------------------------

TEST(glass_catalog_size) {
    seq::Glass::ensureCatalogLoaded();
    // The Java static initialiser loads nine catalogs; a lost entry here means
    // the transliterated catalog data dropped a line.
    CHECK_EQ(static_cast<int>(seq::Glass::glasses().size()), 2481);
}

TEST(glass_lookup_by_name) {
    for (const char *name : {"N-BK7", "S-LAH64", "SF11", "FK50", "F16", "BAF12-Old"}) {
        auto g = seq::Glass::glass_by_name(name);
        CHECK(g != nullptr);
    }
    auto bk7 = seq::Glass::glass_by_name("N-BK7");
    CHECK(bk7 != nullptr);
    CHECK_CLOSE(bk7->nd, 1.5168, 0.0);
    CHECK_CLOSE(bk7->vd, 64.17, 0.0);
    CHECK(bk7->catalog_name.has_value());
    CHECK_STR_EQ(*bk7->catalog_name, "Schott");

    // glass_by_catalog_name normalises the catalog name case-insensitively.
    auto viaCatalog = seq::Glass::glass_by_catalog_name(std::string("schott"), "N-BK7");
    CHECK(viaCatalog != nullptr);
    CHECK_STR_EQ(viaCatalog->toString(), bk7->toString());

    CHECK(seq::Glass::get_catalog_name(std::string("hOyA")).has_value());
    CHECK_STR_EQ(*seq::Glass::get_catalog_name(std::string("hOyA")), "Hoya");
    CHECK(!seq::Glass::get_catalog_name(std::string("nope")).has_value());
    CHECK(!seq::Glass::get_catalog_name(std::nullopt).has_value());
    CHECK(seq::Glass::glass_by_name("does-not-exist") == nullptr);
}

TEST(glass_rindex_from_catalog) {
    auto bk7 = seq::Glass::glass_by_name("N-BK7");
    CHECK_CLOSE(bk7->rindex(seq::Glass::d), 1.5168, 0.0);
    CHECK_CLOSE(bk7->rindex(seq::Glass::C), 1.51432, 0.0);
    CHECK_CLOSE(bk7->rindex(seq::Glass::F), 1.52238, 0.0);
    CHECK_CLOSE(bk7->rindex(seq::Glass::e), 1.51872, 0.0);
    CHECK_CLOSE(bk7->rindex(seq::Glass::g), 1.52668, 0.0);
    // A catalog glass has no data for an arbitrary wavelength.
    CHECK_THROWS(bk7->rindex(500.0), IllegalArgumentException);
}

TEST(glass_rindex_computed_when_unnamed) {
    // The (nd, vd, dpgf) constructor leaves label null, which selects the
    // computed dispersion path rather than the catalog one.
    seq::Glass synth(1.5168, 64.17, 0.0);
    CHECK(!synth.label.has_value());
    CHECK_STR_EQ(synth.toString(), "Glass(nd=1.5168, vd=64.17)");
    CHECK_CLOSE(synth.rindex(seq::Glass::d), 1.5168, 0.0);
    CHECK_CLOSE(synth.rindex(486.1327), 1.5223971858679333, 0.0);
    CHECK_CLOSE(synth.rindex(656.2725), 1.5143349088394855, 0.0);
    CHECK_CLOSE(synth.compute_index_from_nd_vd(500.0), 1.5214324682827298, 0.0);
}

TEST(glass_find_glasses) {
    auto matches = seq::Glass::find_glasses(1.5168, 64.17);
    CHECK_EQ(static_cast<int>(matches.size()), 3);
    CHECK_STR_EQ(*matches[0].glass->catalog_name, "Hoya");
    CHECK_STR_EQ(*matches[0].glass->label, "BSC7");
    CHECK_CLOSE(matches[0].nd_difference, 0.0, 0.0);
    CHECK_CLOSE(matches[0].vd_difference, 0.030000000000001137, 0.0);
    CHECK_CLOSE(matches[0].score, 9.000000000000682E-4, 0.0);
    CHECK(matches[0].exact);

    CHECK_STR_EQ(*matches[1].glass->catalog_name, "Schott");
    CHECK_STR_EQ(*matches[1].glass->label, "BK7");
    CHECK_CLOSE(matches[1].score, 0.0, 0.0);

    CHECK_STR_EQ(*matches[2].glass->catalog_name, "Schott");
    CHECK_STR_EQ(*matches[2].glass->label, "N-BK7");

    // Rejected inputs return empty.
    CHECK_EQ(static_cast<int>(
                 seq::Glass::find_glasses(std::nan(""), 1.0).size()), 0);
    CHECK_EQ(static_cast<int>(
                 seq::Glass::find_glasses(1.5, 64.0, 0.0, 1.0, 3).size()), 0);
    CHECK_EQ(static_cast<int>(
                 seq::Glass::find_glasses(1.5, 64.0, 0.005, 1.0, 0).size()), 0);
}

// ---------------------------------------------------------------------------
// Medium, Air, Gap
// ---------------------------------------------------------------------------

TEST(medium_air_gap) {
    CHECK_STR_EQ(seq::Air::INSTANCE()->toString(), "Air()");
    CHECK_STR_EQ(seq::Medium(1.5).toString(), "Medium(n=1.5)");
    CHECK_STR_EQ(seq::Medium("X", 1.7, "Cat").toString(), "Cat(X)");
    CHECK_CLOSE(seq::Medium(1.5).rindex(400.0), 1.5, 0.0);

    seq::Gap g(2.5, seq::Air::INSTANCE());
    CHECK_STR_EQ(g.toString(), "Gap(t=2.5, medium=Air())");
    g.apply_scale_factor(2.0);
    CHECK_STR_EQ(g.toString(), "Gap(t=5.0, medium=Air())");
    // The default Gap shares the Air singleton.
    CHECK(seq::Gap().medium == seq::Air::INSTANCE());
}

// ---------------------------------------------------------------------------
// Profiles
// ---------------------------------------------------------------------------

TEST(spherical_profile) {
    elem::profiles::Spherical sph(0.02);
    CHECK_STR_EQ(sph.toString(), "Spherical(c=0.02)");
    CHECK_CLOSE(sph.r(), 50.0, 0.0);
    CHECK_CLOSE(sph.f(Vector3(1.0, 2.0, 0.05)), -2.5000000000004186E-5, 0.0);
    CHECK_STR_EQ(sph.df(Vector3(1.0, 2.0, 0.05)).toString(), "[-0.02,-0.04,0.999]");
    CHECK_STR_EQ(sph.normal(Vector3(1.0, 2.0, 0.05)).toString(),
                 "[-0.0199999900000075,-0.039999980000015,0.9989995005003746]");
    CHECK_CLOSE(sph.sag(3.0, 4.0), 0.25062814466900174, 0.0);

    auto ir = sph.intersect(Vector3(0.5, 0.5, -5.0), Vector3(0, 0, 1), 1e-12,
                            util::ZDir::PROPAGATE_RIGHT);
    CHECK_CLOSE(ir.distance, 5.005000250025002, 0.0);
    CHECK_STR_EQ(ir.intersection_point.toString(), "[0.5,0.5,0.005000250025002373]");

    elem::profiles::Spherical flat(0.0);
    CHECK_CLOSE(flat.sag(3.0, 4.0), 0.0, 0.0);
    auto ir2 = flat.intersect(Vector3(0.5, 0.5, -5.0), Vector3(0, 0, 1), 1e-12,
                              util::ZDir::PROPAGATE_RIGHT);
    CHECK_CLOSE(ir2.distance, 5.0, 0.0);

    elem::profiles::Spherical sc(0.02);
    sc.apply_scale_factor(2.0);
    CHECK_CLOSE(sc.cv, 0.01, 0.0);
}

TEST(spherical_sag_throws_when_off_surface) {
    elem::profiles::Spherical sph(0.02);
    // x beyond the radius: the sag square root goes negative.
    CHECK_THROWS(sph.sag(60.0, 0.0), exceptions::TraceMissedSurfaceException);
    // and it must be catchable as the base type the ray tracer catches.
    bool caught = false;
    try {
        sph.sag(60.0, 0.0);
    } catch (const exceptions::TraceException &) {
        caught = true;
    }
    CHECK(caught);
    // and as RuntimeException, which is what ConstraintEdgeThickness catches.
    bool caughtRuntime = false;
    try {
        sph.sag(60.0, 0.0);
    } catch (const redukti::RuntimeException &) {
        caughtRuntime = true;
    }
    CHECK(caughtRuntime);
}

TEST(even_polynomial_profile) {
    std::vector<double> coefs = {1e-6, -2e-9, 3e-12, 0, 0, 0, 0, 0, 0, 0};
    elem::profiles::EvenPolynomial ep(0.02, -0.5, std::nullopt, std::nullopt, coefs);
    CHECK_EQ(ep.max_nonzero_coef, 3);
    CHECK_CLOSE(ep.ec(), 0.5, 0.0);
    CHECK_CLOSE(ep.sag(3.0, 4.0), 0.25033708057498333, 0.0);
    CHECK_CLOSE(ep.f(Vector3(3.0, 4.0, 0.3)), 0.04966291942501666, 0.0);
    CHECK_STR_EQ(ep.df(Vector3(3.0, 4.0, 0.3)).toString(),
                 "[-0.06015599860405027,-0.0802079981387337,1.0]");
    CHECK_STR_EQ(ep.toString(),
                 "EvenPolynomial(c=0.02, cc=-0.5, coefs=[1.0E-6, -2.0E-9, 3.0E-12, "
                 "0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])");
    CHECK_THROWS(ep.sag(200.0, 0.0), exceptions::TraceMissedSurfaceException);

    elem::profiles::EvenPolynomial ep2(0.02, -0.5, std::nullopt, std::nullopt, coefs);
    ep2.apply_scale_factor(2.0);
    CHECK_CLOSE(ep2.cv, 0.01, 0.0);
    CHECK_CLOSE(ep2.coefs[0], 5.0E-7, 0.0);
    CHECK_CLOSE(ep2.coefs[1], -2.5E-10, 0.0);
}

TEST(radial_polynomial_profile) {
    std::vector<double> rcoefs = {0, 1e-5, 0, -2e-8, 0, 0, 0, 0, 0, 0};
    elem::profiles::RadialPolynomial rp(0.02, -0.5, std::nullopt, std::nullopt, rcoefs);
    CHECK_EQ(rp.max_nonzero_coef, 4);
    CHECK_CLOSE(rp.ec, 0.5, 0.0);
    CHECK_CLOSE(rp.cc(), -0.5, 0.0);
    CHECK_CLOSE(rp.sag(3.0, 4.0), 0.25055078369998335, 0.0);
    CHECK_CLOSE(rp.f(Vector3(3.0, 4.0, 0.3)), 0.04944921630001664, 0.0);
    CHECK_STR_EQ(rp.df(Vector3(3.0, 4.0, 0.3)).toString(),
                 "[-0.06020456485405027,-0.0802727531387337,1.0]");
    // On axis r == 0 takes the r_pow = 1.0 branch.
    CHECK_STR_EQ(rp.df(Vector3(0.0, 0.0, 0.0)).toString(), "[-0.0,-0.0,1.0]");
    CHECK_STR_EQ(rp.toString(),
                 "RadialPolynomial(c=0.02, ec=0.5, coefs=[0.0, 1.0E-5, 0.0, -2.0E-8, "
                 "0.0, 0.0, 0.0, 0.0, 0.0, 0.0])");
}

// ---------------------------------------------------------------------------
// Apertures
// ---------------------------------------------------------------------------

TEST(circular_aperture) {
    elem::surface::Circular circ(0.5, -0.25, 0.0, 3.0, false);
    CHECK_STR_EQ(circ.dimension().toString(), "[3.0,3.0]");
    CHECK_CLOSE(circ.max_dimension(), 3.0, 0.0);
    CHECK(circ.point_inside(1.0, 1.0, std::nullopt));
    CHECK(!circ.point_inside(10.0, 0.0, std::nullopt));
    CHECK_STR_EQ(circ.edge_pt_target(Vector2(1.0, 1.0)).toString(),
                 "[2.1213203435596424,2.1213203435596424]");
    auto bb = circ.bounding_box();
    CHECK_STR_EQ(bb.first.toString(), "[-2.5,-3.25]");
    CHECK_STR_EQ(bb.second.toString(), "[3.5,2.75]");

    // An obscuration inverts the inside test.
    elem::surface::Circular obs(0.0, 0.0, 0.0, 3.0, true);
    CHECK(!obs.point_inside(1.0, 1.0, std::nullopt));

    elem::surface::Circular cs(0.5, -0.25, 0.0, 3.0, false);
    cs.apply_scale_factor(2.0);
    CHECK_CLOSE(cs.radius, 6.0, 0.0);
    CHECK_CLOSE(cs.x_offset, 1.0, 0.0);
    CHECK_CLOSE(cs.y_offset, -0.5, 0.0);
}

// ---------------------------------------------------------------------------
// Surface
// ---------------------------------------------------------------------------

TEST(surface_basics) {
    elem::surface::Surface s(seq::InteractMode::TRANSMIT, 0.5, 4.0, nullptr, "L1",
                             std::make_shared<elem::profiles::Spherical>(0.02));
    CHECK_STR_EQ(s.toString(),
                 "Surface(lbl=L1, profile=Spherical(c=0.02), interact_mode='TRANSMIT')");
    CHECK_CLOSE(s.optical_power(), 0.01, 0.0);
    CHECK_CLOSE(s.profile_cv(), 0.02, 0.0);
    CHECK_CLOSE(s.surface_od(), 4.0, 0.0);
    CHECK(s.point_inside(1.0, 1.0, std::nullopt));

    s.set_optical_power(0.02);
    CHECK_CLOSE(s.profile->cv, 0.04, 0.0);
    s.set_optical_power(0.03, 1.0, 1.5);
    CHECK_CLOSE(s.delta_n, 0.5, 0.0);
    CHECK_CLOSE(s.profile->cv, 0.06, 0.0);

    CHECK_STR_EQ(elem::surface::Surface().toString(),
                 "Surface(profile=Spherical(c=0.0), interact_mode='TRANSMIT')");
}

TEST(surface_with_clear_aperture) {
    elem::surface::Surface s(seq::InteractMode::TRANSMIT, 0.5, 4.0, nullptr, "L2",
                             std::make_shared<elem::profiles::Spherical>(0.02));
    s.clear_apertures.push_back(
        std::make_shared<elem::surface::Circular>(0.0, 0.0, 0.0, 2.0, false));

    // The clear aperture, not max_aperture, now governs.
    CHECK_CLOSE(s.surface_od(), 2.0, 0.0);
    CHECK(s.point_inside(1.5, 0.0, std::nullopt));
    CHECK(!s.point_inside(3.0, 0.0, std::nullopt));
    CHECK_STR_EQ(s.edge_pt_target(Vector2(1.0, 0.0)).toString(), "[2.0,0.0]");

    s.set_max_aperture(5.0);
    CHECK_CLOSE(s.max_aperture, 5.0, 0.0);
    CHECK_CLOSE(
        std::static_pointer_cast<elem::surface::Circular>(s.clear_apertures[0])->radius,
        5.0, 0.0);
}

TEST(surface_apply_scale_factor) {
    elem::surface::Surface s(seq::InteractMode::TRANSMIT, 0.5, 4.0, nullptr, "L3",
                             std::make_shared<elem::profiles::Spherical>(0.02));
    s.apply_scale_factor(2.0);
    CHECK_CLOSE(s.max_aperture, 8.0, 0.0);
    CHECK_CLOSE(s.profile->cv, 0.01, 0.0);
}

// ---------------------------------------------------------------------------
// DecenterData
// ---------------------------------------------------------------------------

TEST(decenter_data) {
    elem::surface::DecenterData dd("decenter", 0.5, -0.25, 1.0, 2.0, 3.0);
    dd.update();
    CHECK(dd.rot_mat.has_value());
    auto before = dd.tform_before_surf();
    CHECK(before.rt.has_value());
    // Element-wise rather than by rendered string: this matrix comes from
    // euler2mat, so it inherits the JVM/libm sin-cos disagreement measured in
    // TrigProbeTest. One entry (m01) differs by a single ulp. The entries are
    // all <= 1 in magnitude, so one ulp at magnitude 1 is the right bound.
    const double ONE_ULP_AT_1 = 2.220446049250313e-16;
    const auto &rt = *before.rt;
    CHECK_CLOSE(rt.m00, 0.9980211966240684, ONE_ULP_AT_1);
    CHECK_CLOSE(rt.m01, -0.05230407459247085, ONE_ULP_AT_1);
    CHECK_CLOSE(rt.m02, -0.03489949670250097, ONE_ULP_AT_1);
    CHECK_CLOSE(rt.m10, 0.0529362307009748, ONE_ULP_AT_1);
    CHECK_CLOSE(rt.m11, 0.9984455618447168, ONE_ULP_AT_1);
    CHECK_CLOSE(rt.m12, 0.017441774902830158, ONE_ULP_AT_1);
    CHECK_CLOSE(rt.m20, 0.0339329716976837, ONE_ULP_AT_1);
    CHECK_CLOSE(rt.m21, -0.0192547088685617, ONE_ULP_AT_1);
    CHECK_CLOSE(rt.m22, 0.9992386149554826, ONE_ULP_AT_1);
    CHECK_STR_EQ(before.t.toString(), "[0.5,-0.25,0.0]");

    // dtype "decenter" falls through to the null-rotation branch afterwards.
    auto after = dd.tform_after_surf();
    CHECK(!after.rt.has_value());
    CHECK_STR_EQ(after.t.toString(), "[0.0,0.0,0.0]");

    elem::surface::DecenterData rev("reverse", 0.5, -0.25, 1.0, 2.0, 3.0);
    rev.update();
    auto revBefore = rev.tform_before_surf();
    CHECK(!revBefore.rt.has_value());
    CHECK_STR_EQ(revBefore.t.toString(), "[0.0,0.0,0.0]");

    // No euler rotation leaves rot_mat null.
    elem::surface::DecenterData bend("bend", 0.0, 0.0, 0.0, 0.0, 0.0);
    bend.update();
    CHECK(!bend.rot_mat.has_value());
    CHECK(!bend.tform_after_surf().rt.has_value());
}

// ---------------------------------------------------------------------------
// Lists (Python-style slicing)
// ---------------------------------------------------------------------------

TEST(lists_slice) {
    std::vector<int> v = {0, 1, 2, 3, 4, 5};
    auto all = util::Lists::slice(v, std::nullopt, std::nullopt, std::nullopt);
    CHECK_EQ(static_cast<int>(all.size()), 6);

    auto from2 = util::Lists::from(v, 2);
    CHECK_EQ(static_cast<int>(from2.size()), 4);
    CHECK_EQ(from2[0], 2);

    auto upto3 = util::Lists::upto(v, 3);
    CHECK_EQ(static_cast<int>(upto3.size()), 3);
    CHECK_EQ(upto3[2], 2);

    auto every2 = util::Lists::step(v, 2);
    CHECK_EQ(static_cast<int>(every2.size()), 3);
    CHECK_EQ(every2[1], 2);

    auto rev = util::Lists::step(v, -1);
    CHECK_EQ(static_cast<int>(rev.size()), 6);
    CHECK_EQ(rev[0], 5);

    // Negative indices count from the end.
    auto tail = util::Lists::from(v, -2);
    CHECK_EQ(static_cast<int>(tail.size()), 2);
    CHECK_EQ(tail[0], 4);

    CHECK_EQ(util::Lists::get(v, -1), 5);
    util::Lists::set(v, -1, 99);
    CHECK_EQ(util::Lists::get(v, 5), 99);

    CHECK_THROWS(util::Lists::slice(v, std::nullopt, std::nullopt, 0),
                 IllegalArgumentException);
}
