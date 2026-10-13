// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.spec.WavelengthSetTest
#include "TestHarness.h"

#include "redukti/Exceptions.h"
#include "redukti/rayoptics/seq/Glass.h"
#include "redukti/rayoptics/specs/WvlSpec.h"
#include "redukti/spec/WavelengthSet.h"

#include <string>
#include <vector>

namespace {

using redukti::IllegalArgumentException;
using redukti::rayoptics::seq::Glass;
using redukti::rayoptics::specs::WvlSpec;
using redukti::spec::WavelengthSet;

} // namespace

TEST(wavelength_set_the_standard_sets_are_their_documented_lines_and_weights) {
    CHECK_STR_EQ(WavelengthSet::dLine().letters(), "d");
    CHECK(WavelengthSet::dLine().weights() == (std::vector<double>{1.0}));

    CHECK_STR_EQ(WavelengthSet::standardTriColor().letters(), "d F C");
    CHECK(WavelengthSet::standardTriColor().hasEqualWeights());

    auto five = WavelengthSet::standardFiveColor();
    CHECK_STR_EQ(five.letters(), "d C e F g");
    CHECK(five.weights() == (std::vector<double>{1.0, 0.475, 0.98, 0.49, 0.15}));
    CHECK(five.hasStandardWeights());

    // The measurement weighting is the same five lines, weighted differently: it is
    // greener than the design set, so e carries more than d.
    auto k9 = WavelengthSet::k9Tester();
    CHECK_STR_EQ(k9.letters(), "d C e F g");
    CHECK(k9.weights() == (std::vector<double>{1.0, 0.330, 1.046, 0.654, 0.235}));
    CHECK(!k9.hasStandardWeights());
    CHECK(!k9.hasEqualWeights());
}

/** Wavelengths come from the ported spectral line table, not from a copy of it. */
TEST(wavelength_set_lines_resolve_through_wvl_spec) {
    CHECK(WavelengthSet::standardFiveColor().wavelengths() ==
          (std::vector<double>{WvlSpec::get_wavelength("d"), WvlSpec::get_wavelength("C"),
                               WvlSpec::get_wavelength("e"), WvlSpec::get_wavelength("F"),
                               WvlSpec::get_wavelength("g")}));
    CHECK_EQ(WavelengthSet::eLine().wavelengths()[0], WvlSpec::get_wavelength("e"));
}

TEST(wavelength_set_letters_round_trip_and_unknown_ones_are_rejected) {
    auto e = WavelengthSet::letterOf(WvlSpec::get_wavelength("e"));
    CHECK(e.has_value() && *e == "e");
    CHECK(!WavelengthSet::letterOf(500.0).has_value());

    // Case matters: C and c are not the same line, and only one of them is a line.
    CHECK_EQ(WavelengthSet::wavelengthOf("C"), WvlSpec::get_wavelength("C"));
    bool threw = false;
    std::string message;
    try {
        WavelengthSet::wavelengthOf("c");
    } catch (const IllegalArgumentException &unknown) {
        threw = true;
        message = unknown.what();
    }
    CHECK(threw);
    CHECK(message.find("unknown wavelength 'c'") != std::string::npos);

    // Resolvable in the ported table, but not a Fraunhofer line a trial may name.
    CHECK_THROWS(WavelengthSet::wavelengthOf("He-Ne"), IllegalArgumentException);
}

TEST(wavelength_set_rejects_mismatched_or_negative_weights) {
    CHECK_THROWS(WavelengthSet::ofLetters({"d", "F"}, std::vector<double>{1.0}),
                 IllegalArgumentException);
    CHECK_THROWS(WavelengthSet::ofLetters({"d"}, std::vector<double>{-1.0}),
                 IllegalArgumentException);
    CHECK_THROWS(WavelengthSet::of({0.0}, {1.0}), IllegalArgumentException);
    CHECK_THROWS(WavelengthSet::ofLetters({}), IllegalArgumentException);
}

/**
 * A catalog glass has measured indices at seven lines and throws for anything else, so
 * every set offered here has to stay within them or a trial using it cannot be traced.
 */
TEST(wavelength_set_every_built_in_set_is_traceable_with_a_catalog_glass) {
    auto bk7 = Glass::glass_by_catalog_name("Schott", "N-BK7");
    CHECK(bk7 != nullptr);
    if (bk7 == nullptr)
        return;
    const std::vector<WavelengthSet> sets = {
        WavelengthSet::dLine(), WavelengthSet::eLine(), WavelengthSet::standardTriColor(),
        WavelengthSet::standardFiveColor(), WavelengthSet::k9Tester()};
    for (const auto &set : sets)
        for (double wavelength : set.wavelengths())
            CHECK(bk7->rindex(wavelength) > 1.0);

    // The lines a trial may name are wider than that on purpose: they serve a glass
    // given as index and Abbe number, whose index comes from a dispersion model.
    CHECK_THROWS(bk7->rindex(WavelengthSet::wavelengthOf("i")), IllegalArgumentException);
    Glass modelled(1.5168, 64.17, -0.0009);
    CHECK(modelled.rindex(WavelengthSet::wavelengthOf("i")) > 1.0);
}

// The Java also compares hashCode(), which has no counterpart here.
TEST(wavelength_set_a_set_is_equal_to_another_built_the_same_way) {
    CHECK(WavelengthSet::standardFiveColor() == WavelengthSet::standardFiveColor());
    CHECK(WavelengthSet::standardFiveColor() != WavelengthSet::k9Tester());
}
