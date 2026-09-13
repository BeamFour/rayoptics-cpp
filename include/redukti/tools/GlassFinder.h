// C++ port of org.redukti.tools.GlassFinder
#ifndef REDUKTI_TOOLS_GLASSFINDER_H
#define REDUKTI_TOOLS_GLASSFINDER_H

#include "redukti/rayoptics/seq/Glass.h"
#include "redukti/util/Args.h"

#include <string>

namespace redukti::tools {

/** Adds catalog glass suggestions to OpticalBench lens-data rows. */
class GlassFinder {
public:
    /** Java's record EnrichmentResult. */
    struct EnrichmentResult {
        std::string text;
        int selected;
        int ambiguous;
        int unmatched;
    };

    static EnrichmentResult enrich(const std::string &input, bool force);

    /**
     * @param indexLine which line the prescription's refractive index column is
     *            quoted at. With Glass::IndexLine::E the index is matched against
     *            ne, and a matched surface additionally has its index and Abbe
     *            columns rewritten to the catalog's d line values, so the file
     *            is left consistently on the d line rather than half converted.
     */
    static EnrichmentResult enrich(const std::string &input, bool force,
                                   rayoptics::seq::Glass::IndexLine indexLine);

    /**
     * @param indexLine which line the refractive index column is quoted at
     * @param abbeLine  which line the Abbe number column is quoted at. The two are
     *                  independent: Leica quotes ne with ve, whereas ne paired with
     *                  vd turns up as a transcription slip. Whenever either is the e
     *                  line, a matched surface has its index and Abbe columns
     *                  rewritten to the catalog's d line values, so the file is left
     *                  consistently on the d line rather than half converted.
     */
    static EnrichmentResult enrich(const std::string &input, bool force,
                                   rayoptics::seq::Glass::IndexLine indexLine,
                                   rayoptics::seq::Glass::IndexLine abbeLine);

    /** The body of the Java `main`, minus argv parsing and the usage banner. */
    static void run(const util::Args &arguments);

    static void usage();

private:
    GlassFinder() = delete;
};

} // namespace redukti::tools

#endif // REDUKTI_TOOLS_GLASSFINDER_H
