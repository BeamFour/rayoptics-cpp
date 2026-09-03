// C++ port of org.redukti.util.Args and org.redukti.util.Helper
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_UTIL_ARGS_H
#define REDUKTI_UTIL_ARGS_H

#include "redukti/spec/Prescription.h"

#include <optional>
#include <string>
#include <vector>

namespace redukti::util {

/** Command-line options for the tools. */
class Args {
public:
    int scenario = 0;
    /** Null until --specfile is given; the tools treat that as a usage error. */
    std::optional<std::string> specfile;
    std::string outputType = "layout";
    std::optional<std::string> outputFile;
    std::optional<std::string> outdir;
    bool dumpSystem = false;
    bool use_glass_types = true;
    bool include_lost_rays = false;
    bool only_d_line = false;
    bool do_ray_aberrations = false;
    bool do_mono_chrome_mtfs = false;
    std::vector<int> mtf_freqs = default_mtf_freqs();
    int spot_pattern = 1; // SpotOptions::PATTERN_HEXAPOLAR
    bool auto_size_spots = false;
    bool do_wideangle_layout = false;
    bool force = false;
    spec::VigType vig_type = spec::VigType::SetPupil;
    /** Java uses a Boolean here so "unset" differs from false. */
    std::optional<bool> wide_angle;
    bool generate_java = false;
    bool legacy_notebook = false;
    std::optional<std::string> reference_file;

    static Args parseArguments(const std::vector<std::string> &args);

    static std::vector<int> default_mtf_freqs() { return {10, 30, 50}; }

    static std::vector<int> parse_mtf_freqs(const std::optional<std::string> &value);
    static spec::VigType parse_vig_type(const std::optional<std::string> &value);
    static std::string vig_type_names();
    static int parse_spot_pattern(const std::optional<std::string> &value);
    static std::string spot_pattern_names();
};

/** Path and file helpers, matching org.redukti.util.Helper. */
class Helper {
public:
    static std::string getOutputPath(const Args &arguments);
    static std::string getOutputPath(const Args &arguments, const std::string &extension);
    static std::string getOutputFileWithPath(const std::string &specfile,
                                             const std::string &outputFile,
                                             const std::optional<std::string> &outdir);
    static std::string getFilename(const std::string &specfile);
    static std::string replaceExtension(const std::string &fileName,
                                        const std::string &extension);
    static std::string getOutputPathChangeExt(const std::string &specfile,
                                              const std::string &extension);
    static void createOutputFile(const std::string &outpath, const std::string &content);
};

} // namespace redukti::util

#endif // REDUKTI_UTIL_ARGS_H
