// C++ port of org.redukti.optim.OptimizationTrial
#ifndef REDUKTI_OPTIM_OPTIMIZATIONTRIAL_H
#define REDUKTI_OPTIM_OPTIMIZATIONTRIAL_H

#include "redukti/Exceptions.h"
#include "redukti/optim/OptimizationBuilder.h"
#include "redukti/optim/OptimizationPipeline.h"
#include "redukti/optim/Var.h"
#include "redukti/spec/Prescription.h"

#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace redukti::optim {

/**
 * Reads a `[trial n]` section of a prescription file into an OptimizationBuilder;
 * OptimizationBuilder::toTrial(int) writes one back. The format is documented in
 * Documentation/OPTIMIZER.md.
 *
 * Values mean what they mean to the builder: surfaces are its zero-based positions in
 * [lens data], aspheric coefficients its indices into the coefficient array, and so on.
 * Everything the text alone can settle - syntax, per-field value counts, surface numbers -
 * is checked before the builder is made, and a problem is reported with its line number.
 */
class OptimizationTrial {
public:
    /** A problem with a trial, reported against the line it was found on. */
    class TrialException : public IllegalArgumentException {
    public:
        explicit TrialException(std::string message)
            : IllegalArgumentException(std::move(message)) {}
    };

    /**
     * A trial as read: the prescription built from the trial's own text, and the builder
     * that borrows it. The Java returns the builder alone, its prescription kept alive by
     * the collector; here the caller has to hold the two together.
     */
    struct Trial {
        std::unique_ptr<spec::Prescription> prescription;
        OptimizationBuilder builder;
    };

    /**
     * Reads trial `number` from the text of a prescription file into a builder for that
     * prescription, built from the same text with the trial's wavelength settings.
     *
     * @throws TrialException if the file has no such trial, or the trial has a problem
     */
    static Trial read(const std::string &text, int number, bool useGlassTypes);

    /**
     * Reads `[pipeline number]`, or returns nothing when the number names a trial rather
     * than a pipeline: the two share one numbering.
     *
     * @throws TrialException if the file defines neither, if both claim the number, or if
     *                        the pipeline has a problem
     */
    static std::optional<OptimizationPipeline> readPipeline(const std::string &text,
                                                            int number);

    /** How a variable is named when a trial runs: by surface position, as in the trial. */
    static std::string describe(const Var &variable);

    // ------------------------------------------------------------------
    // Shared with OptimizationBuilder::toTrial; package-private in the Java
    // ------------------------------------------------------------------

    static std::string paraxialName(int paraxId);

    /** A number as a trial writes it: whole numbers without a fraction, others as Java does. */
    static std::string format(double value);
    static std::string format(const std::vector<double> &values);
    static std::string format(const std::vector<int> &values);

    /** One line of a trial or pipeline: the keyword padded to a column, then its values. */
    static void line(std::string &sb, const std::string &key, const std::string &values);

    /** SetPupil as set-pupil: the spelling --vig-type and a trial use. */
    static std::string kebab(const std::string &name);

    /** Splits a line on tabs exactly as the prescription reader does. */
    static std::vector<std::string> splitTabs(std::string line);

private:
    OptimizationTrial() = delete;
};

} // namespace redukti::optim

#endif // REDUKTI_OPTIM_OPTIMIZATIONTRIAL_H
