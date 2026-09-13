// C++ port of org.redukti.util.Args and org.redukti.util.Helper
#ifndef REDUKTI_UTIL_ARGS_H
#define REDUKTI_UTIL_ARGS_H

#include "redukti/rayoptics/analysis/PupilMapAnalysis.h"
#include "redukti/rayoptics/seq/Glass.h"
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
    // The Java's patent / example (fetching a prescription from the
    // PhotonsToPhotos Optical Bench) are deliberately not ported: they need an
    // HTTP client, and this port has no third-party dependencies.
    std::optional<std::string> outputFile;
    std::optional<std::string> outdir;
    bool use_glass_types = true;
    bool only_d_line = false;
    bool do_ray_aberrations = false;
    bool do_mono_chrome_mtfs = false;
    /**
     * Spatial frequencies in cycles/mm at which the MTF is reported. The
     * default is what every report under Examples/ uses, so that lenses stay
     * comparable across reports; override it only when checking a design
     * against a manufacturer's own choice of frequencies.
     */
    std::vector<int> mtf_freqs = default_mtf_freqs();
    /**
     * Ray pattern used for spot diagrams, one of the
     * SpotOptions.PATTERN_* constants. Defaults to hexapolar, which is
     * also SpotOptions' own default.
     */
    int spot_pattern = 1; // SpotOptions::PATTERN_HEXAPOLAR
    /** Number of samples along each dimension of a rectangular spot grid. */
    int spot_grid_size = 64;
    bool auto_size_spots = false;
    /**
     * Run the glass type matcher over the prescription before analysing it, so
     * that surfaces quoting only nd and vd pick up a real catalog glass and its
     * full dispersion curve.
     */
    bool assign_glass_types = false;
    /**
     * With assign_glass_types, also write the enriched prescription back over
     * the input file. Off by default: assigning glass types is an analysis
     * choice, and overwriting the author's input should be asked for.
     */
    bool update_specfile = false;
    /**
     * Line the prescription's refractive index column is quoted at, "d" or "e".
     * Some patents tabulate the index at the e line while still quoting the Abbe
     * number as vd; matching on the wrong line finds nothing.
     */
    std::string index_line = "d";
    /**
     * Line the prescription's Abbe number column is quoted at, "d" or "e".
     * Independent of index_line: Leica patents quote ne with ve, while a
     * prescription pairing ne with vd is usually a transcription slip.
     */
    std::string abbe_line = "d";
    /**
     * Run the routine airspace optimization before reporting: the back focus on
     * a prime, the variable airspaces other than the back focus on a zoom.
     */
    bool optimize = false;
    /**
     * Objective for `optimize`: "contrast" or "mtf". Contrast is the default
     * because the geometric MTF merit surface is rough at the scale the solver
     * steps, so a solve driven by it stalls in a local minimum.
     */
    std::string optimize_goal = "contrast";
    /**
     * Number of the [trial n] or [pipeline n] section to run before reporting, from
     * --optimize n. Empty when none was asked for.
     */
    std::optional<int> optimize_trial;
    bool force = false;
    /**
     * Vignetting calculation applied once the model is built. Defaults to the
     * value every tool in this module already hard-codes, so wiring an existing
     * tool up to this field is a no-op.
     *
     * In LensTool2 this settles the models the ANALYSIS outputs are computed from - the
     * spot diagrams, the MTF, the ray aberration fans and the vignetting and paraxial
     * dumps. The layout diagrams and the routine optimization set their own; see
     * LensTool2::LAYOUT_VIG_TYPE and LensTool2::OPTIMIZATION_VIG_TYPE.
     */
    spec::VigType vig_type = spec::VigType::SetPupil;
    /**
     * Write the measured pupil maps beside the other reports: which part of each field's
     * pupil the lens passes, and which surface blocks the rest. See
     * rayoptics::analysis::PupilMapAnalysis.
     */
    bool output_pupil_maps = false;
    /** Samples per axis in a pupil map; the cost is the square of this. */
    int pupil_map_samples = rayoptics::analysis::PupilMapAnalysis::DEFAULT_NUM_SAMPLES;
    /**
     * Selects the chief ray aiming algorithm. TRUE aims with a real ray trace
     * at the entrance pupil (what the model calls a wide angle system), FALSE
     * uses paraxial aiming. Null leaves it derived from the half angle of view.
     *
     * Real ray aiming is slower but is what makes very wide angle lenses trace
     * correctly, so it is the default in LensTool2.
     */
    /** The Java uses a Boolean here so "unset" differs from false; std::optional plays that part. */
    std::optional<bool> real_ray_aiming;
    /** Emit Java model building code rather than Python. */
    bool generate_java = false;
    /** Emit the original plotting notebook script rather than a comparison model. */
    bool legacy_notebook = false;
    /**
     * Path to the upstream reference values produced by dump_reference.py. When
     * set, the exporter emits a JUnit regression test instead of a model builder.
     */
    std::optional<std::string> reference_file;

    static Args parseArguments(const std::vector<std::string> &args);

    /** The org.redukti.rayoptics.seq.Glass.IndexLine this maps to. */
    rayoptics::seq::Glass::IndexLine index_line_value() const;

    /** The org.redukti.rayoptics.seq.Glass.IndexLine the Abbe column maps to. */
    rayoptics::seq::Glass::IndexLine abbe_line_value() const;

    /** Accepts d or e, rejecting anything else rather than defaulting. */
    static std::string parse_index_line(const std::optional<std::string> &value);

    /** Accepts contrast or mtf, rejecting anything else rather than defaulting. */
    static std::string parse_optimize_goal(const std::optional<std::string> &value);

    /** Accepts the number of a [trial n] section. */
    static int parse_trial_number(const std::string &value);

    /** The MTF frequencies used by every report under Examples/. */
    static std::vector<int> default_mtf_freqs() { return {10, 30, 50}; }

    /**
     * Parses a comma separated list of spatial frequencies in cycles/mm, e.g.
     * "10,20,40". As with the other typed options a bad value is rejected
     * rather than defaulted, since silently reporting the standard 10/30/50
     * would look like a valid answer to a question that was never asked.
     */
    static std::vector<int> parse_mtf_freqs(const std::optional<std::string> &value);
    /**
     * Accepts either the enum constant (SetPupil) or its kebab-case spelling
     * (set-pupil), ignoring case and any - or _ separators.
     *
     * An unrecognized value is rejected rather than quietly falling back to the
     * default: the vignetting type changes the model that gets built, so a typo
     * would otherwise shift every number downstream with nothing to show for it.
     */
    static spec::VigType parse_vig_type(const std::optional<std::string> &value);
    /** The accepted --vig-type spellings, for usage and error messages. */
    static std::string vig_type_names();
    /** The Java enum constant name, as VigType.name() gives it: SetPupil, SetVig and so on. */
    static std::string vig_type_name(spec::VigType value);
    /**
     * Accepts hex, grid or gaussian, returning the matching
     * SpotOptions.PATTERN_* constant. As with --vig-type an
     * unrecognized value is rejected rather than defaulted, since the sampling
     * pattern changes every spot number it produces.
     */
    static int parse_spot_pattern(const std::optional<std::string> &value);
    /** The accepted --use-spot-pattern spellings, for usage and error messages. */
    static std::string spot_pattern_names();

private:
    static int parse_positive_int(const std::string &option,
                                  const std::optional<std::string> &value);
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
