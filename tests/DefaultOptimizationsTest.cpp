// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.tools.DefaultOptimizationsTest.
//
// Which airspaces the routine optimization picks out of a prescription.
#include "TestHarness.h"

#include "redukti/Exceptions.h"
#include "redukti/importers/OpticalBenchDataImporter.h"
#include "redukti/spec/Prescription.h"
#include "redukti/spec/SurfaceType.h"
#include "redukti/tools/DefaultOptimizations.h"
#include "redukti/tools/LensTool2.h"

#include <string>
#include <vector>

namespace {

using redukti::IllegalArgumentException;
using redukti::importers::OpticalBenchDataImporter;
using redukti::spec::Prescription;
using redukti::tools::DefaultOptimizations;
using redukti::tools::LensTool2;

/**
 * A two configuration zoom shaped like the patents that model the sensor stack:
 * two airspaces that move with the zoom, then a cover glass, then the fixed gap
 * between cover glass and image. 'Bf' is a multi valued row holding one value, as
 * these files write it.
 */
const std::string ZOOM_WITH_COVER_GLASS = "[descriptive data]\n"
                                          "title\tTest zoom\n"
                                          "[variable distances]\n"
                                          "Focal Length\t50\t60\n"
                                          "F-Number\t4\t5\n"
                                          "Angle of View\t40\t35\n"
                                          "d2\t10\t6\n"
                                          "d4\t20\t26\n"
                                          "Bf\t1\t1\n"
                                          "[lens data]\n"
                                          "1\t50\t4\t1.5\t20\t60\n"
                                          "2\t-50\td2\t\t20\n"
                                          "3\t80\t3\t1.6\t20\t50\n"
                                          "4\t-80\td4\t\t20\n"
                                          "5\tCG\t1.6\t1.5168\t44\t64\tBK7\n"
                                          "6\tCG\tBf\t\t44\n"
                                          "[report data]\n"
                                          "lens name\tTest zoom\n"
                                          "scenarios\t0\t1\n"
                                          "names\tWide\tLong\n";

/**
 * A three scenario zoom whose optical rows are all filled in, so that only the
 * thickness rows decide anything. 'Bf' holds the same value at the outer two
 * scenarios and a placeholder at the middle one, which is how these files write a
 * scenario the patent never tabulated.
 */
const std::string ZOOM_WITH_THIRD_SCENARIO = "[descriptive data]\n"
                                             "title\tTest zoom\n"
                                             "[variable distances]\n"
                                             "Focal Length\t50\t55\t60\n"
                                             "F-Number\t4\t4.5\t5\n"
                                             "Angle of View\t40\t37\t35\n"
                                             "d2\t10\t8\t6\n"
                                             "d4\t20\t23\t26\n"
                                             "Bf\t1\tMIDDLE\t1\n"
                                             "[lens data]\n"
                                             "1\t50\t4\t1.5\t20\t60\n"
                                             "2\t-50\td2\t\t20\n"
                                             "3\t80\t3\t1.6\t20\t50\n"
                                             "4\t-80\td4\t\t20\n"
                                             "5\tCG\t1.6\t1.5168\t44\t64\tBK7\n"
                                             "6\tCG\tBf\t\t44\n"
                                             "[report data]\n"
                                             "lens name\tTest zoom\n"
                                             "scenarios\tSELECTED\n"
                                             "names\tWide\tLong\n";

/** Java's String.replace: every occurrence. */
std::string replace(std::string text, const std::string &from, const std::string &to) {
    std::size_t pos = 0;
    while ((pos = text.find(from, pos)) != std::string::npos) {
        text.replace(pos, from.size(), to);
        pos += to.size();
    }
    return text;
}

Prescription prescription(const std::string &text) {
    OpticalBenchDataImporter::LensSpecifications specs;
    specs.parse_buffer(text);
    return LensTool2::createPrescription(specs, true, false);
}

std::string zoomSelecting(const std::string &middleValue, const std::string &selected) {
    return replace(replace(ZOOM_WITH_THIRD_SCENARIO, "MIDDLE", middleValue), "SELECTED",
                   selected);
}

bool contains(const std::string &text, const std::string &part) {
    return text.find(part) != std::string::npos;
}

} // namespace

TEST(zoom_varies_the_airspaces_that_move_with_the_zoom_and_not_the_fixed_sensor_gap) {
    auto p = prescription(ZOOM_WITH_COVER_GLASS);
    CHECK_EQ(p.get_num_configurations(), 2);
    // The airspace in front of the cover glass is the back focus by position, but on
    // this zoom it is also one of the airspaces that moves, so it has to be varied.
    CHECK_EQ(DefaultOptimizations::findBackFocusSurface(p), 3);
    CHECK(DefaultOptimizations::findVariableThicknesses(p) == (std::vector<int>{1, 3}));
}

TEST(an_airspace_holding_one_value_in_every_configuration_is_not_a_zoom_airspace) {
    // Bf is referenced from a multi valued row, so it carries a thickness per
    // configuration; holding the same value in both is what keeps it out.
    auto p = prescription(ZOOM_WITH_COVER_GLASS);
    const auto &surfaces = p.get_surfaces();
    CHECK(surfaces[5]._thickness_by_scenario == (std::vector<double>{1.0, 1.0}));
    // the fixed cover glass to image gap must not be varied
    for (int surface : DefaultOptimizations::findVariableThicknesses(p))
        CHECK(surface != 5);
}

TEST(a_zoom_without_a_cover_glass_varies_its_rear_airspace) {
    // With no cover glass the rear airspace is both the back focus and the airspace
    // that moves with the zoom; leaving it out would drop the strongest variable.
    auto p = prescription(replace(
        replace(ZOOM_WITH_COVER_GLASS,
                "5\tCG\t1.6\t1.5168\t44\t64\tBK7\n"
                "6\tCG\tBf\t\t44\n",
                ""),
        "Bf\t1\t1\n", ""));
    CHECK_EQ(p.get_surfaces().size(), std::size_t{4});
    CHECK_EQ(DefaultOptimizations::findBackFocusSurface(p), 3);
    CHECK(DefaultOptimizations::findVariableThicknesses(p) == (std::vector<int>{1, 3}));
}

TEST(only_the_scenarios_the_report_selects_decide_whether_an_airspace_varies) {
    // The middle scenario is not one of the two being reported, so whatever stands in
    // it has no say. Bf reads as one value across the selected pair either way.
    for (const std::string middle : {"undefined", "", "17"}) {
        auto p = prescription(zoomSelecting(middle, "0\t2"));
        CHECK_EQ(p.get_num_configurations(), 2);
        const auto &surfaces = p.get_surfaces();
        CHECK(surfaces[1]._thickness_by_scenario == (std::vector<double>{10.0, 6.0}));
        CHECK(surfaces[3]._thickness_by_scenario == (std::vector<double>{20.0, 26.0}));
        // the middle scenario must not leak into the selected pair
        CHECK(surfaces[5]._thickness_by_scenario == (std::vector<double>{1.0, 1.0}));
        CHECK(DefaultOptimizations::findVariableThicknesses(p) ==
              (std::vector<int>{1, 3}));
    }
}

TEST(a_missing_value_in_a_selected_scenario_is_rejected_before_any_airspace_is_picked) {
    // Selecting the scenario the patent left blank is an error in the input, caught
    // while the prescription is built - so a missing value can never reach the
    // selection above as the zero that a lenient parse would produce, where it would
    // read as a thickness that varies.
    for (const std::string middle : {"undefined", "", "nonsense", "NaN"}) {
        std::string message;
        bool thrown = false;
        try {
            prescription(zoomSelecting(middle, "0\t1"));
        } catch (const IllegalArgumentException &e) {
            thrown = true;
            message = e.what();
        }
        // the middle scenario must not be accepted
        CHECK(thrown);
        CHECK(contains(message, "'Bf'"));
        CHECK(contains(message, "scenario 1"));
    }
}

TEST(a_prescription_with_one_configuration_has_no_zoom_airspaces) {
    auto p = prescription(replace(ZOOM_WITH_COVER_GLASS,
                                  "scenarios\t0\t1\nnames\tWide\tLong",
                                  "scenarios\t0\nnames\tWide"));
    CHECK_EQ(p.get_num_configurations(), 1);
    CHECK(DefaultOptimizations::findVariableThicknesses(p).empty());
    // A prime is optimized on its back focus instead, which is still located.
    CHECK_EQ(DefaultOptimizations::findBackFocusSurface(p), 3);
}
