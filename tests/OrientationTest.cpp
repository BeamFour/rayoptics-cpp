// C++ port of org.redukti.rayoptics.util.OrientationTest and
// org.redukti.rayoptics.analysis.AnalysisTraceOptionsTest.
//
// Both are small guards on constants and defaults that other code indexes by,
// so they live together rather than in a file each.
#include "TestHarness.h"

#include "redukti/rayoptics/analysis/ContrastAnalysis.h"
#include "redukti/rayoptics/analysis/SpotAnalysis.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/rayoptics/specs/WvlSpec.h"
#include "redukti/rayoptics/util/Orientation.h"

#include <string>

namespace {

namespace Orientation = redukti::rayoptics::util::Orientation;
using redukti::rayoptics::analysis::ContrastOptions;
using redukti::rayoptics::analysis::SpotOptions;

} // namespace

// ===========================================================================
// OrientationTest
// ===========================================================================

TEST(orientation_checked_returns_the_two_meridians_and_rejects_everything_else) {
    using redukti::IllegalArgumentException;
    CHECK_EQ(Orientation::checked(Orientation::SAGITTAL), Orientation::SAGITTAL);
    CHECK_EQ(Orientation::checked(Orientation::TANGENTIAL), Orientation::TANGENTIAL);

    bool threw = false;
    std::string message;
    try {
        Orientation::checked(2);
    } catch (const IllegalArgumentException &e) {
        threw = true;
        message = e.what();
    }
    CHECK(threw);
    CHECK(message.find("SAGITTAL") != std::string::npos);
    CHECK_THROWS(Orientation::checked(-1), IllegalArgumentException);
}

/**
 * The ray-coordinate aliases and the MTF orientations index the same analysis
 * arrays from different callers - a spot deviation goal takes the x intercepts
 * for Orientation::X while a contrast goal takes the sagittal residual for
 * Orientation::SAGITTAL. They have to stay equal or the two would disagree
 * about which meridian index 0 is.
 */
TEST(orientation_ray_coordinate_aliases_match_the_meridians_they_name) {
    CHECK_EQ(Orientation::X, Orientation::SAGITTAL);
    CHECK_EQ(Orientation::Y, Orientation::TANGENTIAL);
}

/** for (xy = 0; xy < COUNT; xy++) must cover both meridians and nothing else. */
TEST(orientation_count_bounds_a_loop_over_both_meridians) {
    CHECK_EQ(Orientation::COUNT, 2);
    CHECK(Orientation::SAGITTAL < Orientation::COUNT);
    CHECK(Orientation::TANGENTIAL < Orientation::COUNT);
}

TEST(orientation_names_the_meridian_for_labels_and_descriptions) {
    CHECK_STR_EQ(Orientation::name(Orientation::SAGITTAL), "sag");
    CHECK_STR_EQ(Orientation::name(Orientation::TANGENTIAL), "tan");
}

// ===========================================================================
// AnalysisTraceOptionsTest
// ===========================================================================

TEST(analysis_spot_checks_apertures_by_default_and_can_disable_them) {
    SpotOptions options;
    CHECK(options._trace_options.check_apertures);
    CHECK(!options.check_apertures(false)._trace_options.check_apertures);
}

TEST(analysis_contrast_skips_aperture_checks_by_default_and_can_enable_them) {
    ContrastOptions options(40.0);
    CHECK(!options.traceOptions.check_apertures);
    CHECK(options.check_apertures(true).traceOptions.check_apertures);
}

// ===========================================================================
// Spectral line lookup
// ===========================================================================

/**
 * No Java counterpart: this pins a deliberate divergence from it.
 *
 * Sodium D and helium d are different lines. The Java uppercases every key into
 * a single map before looking anything up, so the two collide and it answers
 * 587.5618 for both -- verified against the JVM. That silently returns the
 * wrong line for "D", so this port matches case exactly and keeps them apart.
 *
 * OpticalSpecs defaults its WvlSpec to WvlWt("d", 1.0), so the lowercase value
 * is what every default-wavelength model traces at; it was 589.2938 until this
 * was found, putting those models 1.7 nm off.
 */
TEST(wvlspec_distinguishes_sodium_D_from_helium_d) {
    using redukti::rayoptics::specs::WvlSpec;
    CHECK_EQ(WvlSpec::get_wavelength("d"), 587.5618);
    CHECK_EQ(WvlSpec::get_wavelength("D"), 589.2938);
    // Other spellings still resolve case-insensitively, as the Java does.
    CHECK_EQ(WvlSpec::get_wavelength("He-Ne"), 632.8);
    CHECK_EQ(WvlSpec::get_wavelength("HE-NE"), 632.8);
    CHECK_EQ(WvlSpec::get_wavelength("F"), 486.1327);
    CHECK_THROWS(WvlSpec::get_wavelength("nosuchline"),
                 redukti::IllegalArgumentException);
}

/** The default model traces the helium d line, not sodium D. */
TEST(wvlspec_default_optical_spec_uses_the_helium_d_line) {
    redukti::rayoptics::optical::OpticalModel opm;
    const auto &wvls = opm.optical_spec->wvls->wavelengths;
    CHECK_EQ(static_cast<int>(wvls.size()), 1);
    CHECK_EQ(wvls[0], 587.5618);
}
