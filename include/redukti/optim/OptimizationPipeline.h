// C++ port of org.redukti.optim.OptimizationPipeline
#ifndef REDUKTI_OPTIM_OPTIMIZATIONPIPELINE_H
#define REDUKTI_OPTIM_OPTIMIZATIONPIPELINE_H

#include <optional>
#include <string>
#include <vector>

namespace redukti::optim {

/**
 * A `[pipeline n]` section: trials run one after another, each starting from the design
 * the one before it produced.
 *
 * A zoom is the case it was made for. A configuration is optimized at a time, so the
 * wide end is optimized first and the tele end then starts from that result; since the
 * two share every radius and glass thickness, the second stage partly undoes the first,
 * and the pipeline is written to alternate - `trials 1 2 1 2` - until the design
 * settles. Nothing is specific to zooms, though: a coarse stage followed by a fine one,
 * or spot goals followed by contrast, chain the same way.
 *
 * Trial and pipeline numbers share one numbering, so `--optimize 3` runs whichever of
 * them the file defines.
 */
class OptimizationPipeline {
public:
    OptimizationPipeline(int number, std::optional<std::string> description,
                         std::optional<std::string> outdir, std::vector<int> trials);

    int number() const { return _number; }

    /** The pipeline's description, or empty when it has none. */
    const std::optional<std::string> &description() const { return _description; }

    /**
     * Where LensTool2 puts the optimized prescription and the report, relative to the
     * prescription's file unless absolute, or empty when the pipeline does not say. The
     * stages' own outdirs are not used: a pipeline writes one result.
     */
    const std::optional<std::string> &outdir() const { return _outdir; }

    /** The trials to run, in order; a trial may appear more than once. */
    const std::vector<int> &trials() const { return _trials; }

    /** The stages as the trials line lists them, "1 2 1 2", for reporting. */
    std::string trialsText() const;

    /** The trials the pipeline names, each once, in order of number. */
    std::vector<int> distinctTrials() const;

    /** This pipeline as a `[pipeline n]` section, which OptimizationTrial::readPipeline reads back. */
    std::string toPipeline() const;

private:
    int _number;
    std::optional<std::string> _description;
    std::optional<std::string> _outdir;
    std::vector<int> _trials;
};

} // namespace redukti::optim

#endif // REDUKTI_OPTIM_OPTIMIZATIONPIPELINE_H
