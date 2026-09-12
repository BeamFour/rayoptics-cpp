// Command-line entry point for LensTool2, matching the Java `main`.
//
// The Java's --patent / --example (fetch the prescription from the
// PhotonsToPhotos Optical Bench) are not ported: they need an HTTP client, and
// this port has no third-party dependencies.
#include "redukti/tools/LensTool2.h"
#include "redukti/util/Args.h"

#include <exception>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char **argv) {
    using redukti::util::Args;
    std::vector<std::string> args;
    for (int i = 1; i < argc; i++)
        args.push_back(argv[i]);
    Args arguments;
    try {
        arguments = Args::parseArguments(args);
    } catch (const std::exception &e) {
        // The Java lets the exception escape main, which exits non-zero.
        std::cerr << e.what() << std::endl;
        return 1;
    }
    if (!arguments.specfile.has_value()) {
        // The trailing backslash-newline pairs mirror the Java usage text,
        // which wraps the synopsis the way a shell continuation would.
        std::cerr << "Usage: --specfile inputfile [--outdir dir] \\\n";
        std::cerr << "       [--only-d-line] [--dont-use-glass-types] \\\n";
        std::cerr << "       [--output-ray-aberration-plots] [--output-wavelength-mtfs] "
                     "[--auto-size-spot-diagrams] \\\n";
        std::cerr << "       [--use-spot-pattern " << Args::spot_pattern_names()
                  << "] [--spot-grid-size count] [--vig-type " << Args::vig_type_names()
                  << "] \\\n";
        std::cerr << "       [--real-ray-aiming|--paraxial-ray-aiming] [--mtf "
                     "freq,freq,...] \\\n";
        std::cerr << "       [--assign-glass-types [--index-line d|e] [--force] "
                     "[--update-specfile]] [--optimize [--optimize-goal contrast|mtf] | "
                     "--optimize trial]\n";
        std::cerr << "       --assign-glass-types matches each surface's nd/vd to a "
                     "catalog glass for this run;\n";
        std::cerr << "         --force re-matches surfaces that already name a glass, "
                     "--update-specfile writes the result back to the specfile\n";
        std::cerr << "         --index-line e when the prescription quotes the refractive "
                     "index at the e line rather than the d line\n";
        std::cerr << "       --optimize varies the back focus on a prime, or the other "
                     "variable airspaces on a zoom, at the central field\n";
        std::cerr << "       --optimize-goal defaults to contrast; mtf uses the geometric "
                     "MTF directly, which stalls more easily\n";
        std::cerr << "       --optimize n runs the specfile's [trial n] section, writes "
                     "the result as <specfile>-trial<n>.txt and reports on it\n";
        std::cerr << "         a [pipeline n] section runs its trials in order, each "
                     "starting from the last result, and writes "
                     "<specfile>-pipeline<n>.txt\n";
        std::cerr << "       --mtf takes spatial frequencies in cycles/mm and defaults to "
                     "10,30,50, which is what the reports under Examples/ use\n";
        std::cerr << "       --real-ray-aiming aims the chief ray by tracing a real ray at "
                     "the entrance pupil, --paraxial-ray-aiming uses paraxial aiming; real "
                     "is the default\n";
        std::cerr << "       Output files are created alongside the specfile unless "
                     "--outdir is given\n";
        return 1;
    }
    if (arguments.update_specfile && !arguments.assign_glass_types) {
        std::cerr << "--update-specfile only applies with --assign-glass-types" << std::endl;
        return 1;
    }
    try {
        redukti::tools::LensTool2::run(arguments, redukti::tools::LensTool2::today());
    } catch (const std::exception &e) {
        // The Java prints the message and a stack trace, then falls through to
        // a normal exit; the message is what matters to a caller.
        std::cerr << "Failed due to: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
