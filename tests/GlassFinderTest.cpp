// Port of org.redukti.tools.GlassFinderTest.
#include "TestHarness.h"

#include "redukti/rayoptics/seq/Glass.h"
#include "redukti/tools/GlassFinder.h"

#include <string>

namespace {

using redukti::rayoptics::seq::Glass;
using redukti::tools::GlassFinder;

bool contains(const std::string &s, const std::string &part) {
    return s.find(part) != std::string::npos;
}

/** Occurrences of `part`: Java's s.split(part, -1).length - 1. */
int countOf(const std::string &s, const std::string &part) {
    int n = 0;
    for (std::size_t at = s.find(part); at != std::string::npos;
         at = s.find(part, at + part.size()))
        n++;
    return n;
}

} // namespace

TEST(glassfinder_selects_an_exact_match_and_preserves_other_sections) {
    std::string input = "[descriptive data]\n"
                        "title\tExample\n"
                        "[lens data]\n"
                        "1\t50\t4\t1.58267\t20\t46.48\n"
                        "[notes]\n"
                        "unchanged\n";

    GlassFinder::EnrichmentResult result = GlassFinder::enrich(input, false);
    auto expected = Glass::find_glasses(1.58267, 46.48).at(0).glass;

    CHECK_EQ(result.selected, 1);
    CHECK_EQ(result.ambiguous, 0);
    CHECK(contains(result.text, "1\t50\t4\t1.58267\t20\t46.48\t" + *expected->label + "\t" +
                                    *expected->catalog_name));
    CHECK(contains(result.text, "[notes]\nunchanged\n"));
}

TEST(glassfinder_appends_candidates_without_assigning_an_ambiguous_approximate_match) {
    std::string input = "[lens data]\n"
                        "1\t50\t4\t1.58270\t20\t46.50\n";

    GlassFinder::EnrichmentResult result = GlassFinder::enrich(input, false);

    CHECK_EQ(result.selected, 0);
    CHECK_EQ(result.ambiguous, 1);
    CHECK(contains(result.text, "\t\tcandidate="));
    CHECK(countOf(result.text, "candidate=") >= 2);
}

TEST(glassfinder_leaves_existing_glass_assignments_untouched) {
    std::string input = "[lens data]\n"
                        "1\t50\t4\t1.58267\t20\t46.19\tBAF3\tSchott\n";

    GlassFinder::EnrichmentResult result = GlassFinder::enrich(input, false);

    CHECK_STR_EQ(result.text, input);
    CHECK_EQ(result.selected, 0);
    CHECK_EQ(result.ambiguous, 0);
}

TEST(glassfinder_enrichment_is_idempotent_for_candidate_rows) {
    std::string input = "[lens data]\n"
                        "1\t50\t4\t1.58270\t20\t46.50\n";

    std::string first = GlassFinder::enrich(input, false).text;
    std::string second = GlassFinder::enrich(first, false).text;

    CHECK_STR_EQ(second, first);
}

TEST(glassfinder_promotes_existing_candidates_when_nd_matches_and_vd_is_within_one_decimal) {
    std::string input =
        "[lens data]\n"
        "5\t90.458\t10.27\t1.497\t69.92\t81.6\t\t"
        "\tcandidate=Hoya/FCD1,nd=1.49700,vd=81.61,dnd=0.00000,dvd=0.01"
        "\tcandidate=Ohara/FPL51,nd=1.49700,vd=81.61,dnd=0.00000,dvd=0.01"
        "\tcandidate=Schott/N-PK52A,nd=1.49700,vd=81.61,dnd=0.00000,dvd=0.01\n";

    GlassFinder::EnrichmentResult result = GlassFinder::enrich(input, false);

    CHECK_EQ(result.selected, 1);
    CHECK_EQ(result.ambiguous, 0);
    CHECK(contains(result.text, "5\t90.458\t10.27\t1.497\t69.92\t81.6\tFCD1\tHoya\t"));
    CHECK_EQ(countOf(result.text, "candidate="), 3);
}

TEST(glassfinder_reports_and_preserves_unmatched_glass_data) {
    std::string input = "[lens data]\n"
                        "1\t50\t4\t1.10\t20\t10.0\n";

    GlassFinder::EnrichmentResult result = GlassFinder::enrich(input, false);

    CHECK_EQ(result.unmatched, 1);
    CHECK_STR_EQ(result.text, input);
}
