// Command-line entry point for GlassFinder, matching the Java `main`.
#include "redukti/tools/GlassFinder.h"
#include "redukti/util/Args.h"

#include <exception>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char **argv) {
    std::vector<std::string> args;
    for (int i = 1; i < argc; i++)
        args.push_back(argv[i]);
    try {
        auto arguments = redukti::util::Args::parseArguments(args);
        if (!arguments.specfile.has_value()) {
            redukti::tools::GlassFinder::usage();
            return 2;
        }
        redukti::tools::GlassFinder::run(arguments);
    } catch (const std::exception &e) {
        // The Java lets the exception escape main, which exits non-zero.
        std::cerr << e.what() << std::endl;
        return 1;
    }
    return 0;
}
