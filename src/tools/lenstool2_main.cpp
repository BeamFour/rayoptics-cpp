// Command-line entry point for LensTool2, matching the Java `main`.
#include "redukti/tools/LensTool2.h"
#include "redukti/util/Args.h"

#include <exception>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char **argv) {
    std::vector<std::string> args;
    for (int i = 1; i < argc; i++)
        args.push_back(argv[i]);
    auto arguments = redukti::util::Args::parseArguments(args);
    if (!arguments.specfile.has_value()) {
        std::cerr
            << "Usage: --specfile inputfile [--scenario num] [--dump-system] "
               "[--only-d-line] [-o outfilename] [--dont-use-glass-types] \\n";
        std::cerr << "       [--output-ray-aberration-plots] [--output-wavelength-mtfs] "
                     "[--auto-size-spot-diagrams] [--do-wideangle-layout] \\n";
        std::cerr << "       [--use-spot-pattern "
                  << redukti::util::Args::spot_pattern_names() << "] [--vig-type "
                  << redukti::util::Args::vig_type_names()
                  << "] [--wide-angle|--no-wide-angle] \\n";
        std::cerr << "       [--mtf freq,freq,...]\n";
        std::cerr << "       --scenario defaults to 0\n";
        std::cerr << "       --mtf takes spatial frequencies in cycles/mm and defaults "
                     "to 10,30,50, which is what the reports under Examples/ use\n";
        std::cerr << "       Output file will be created in the same location as the "
                     "specfile\n";
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
