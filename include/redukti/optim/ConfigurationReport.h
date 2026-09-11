// C++ port of org.redukti.optim.ConfigurationReport
#ifndef REDUKTI_OPTIM_CONFIGURATIONREPORT_H
#define REDUKTI_OPTIM_CONFIGURATIONREPORT_H

#include "redukti/optim/Analysis.h"
#include "redukti/spec/Prescription.h"

#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace redukti::optim {

/**
 * Evaluates every configuration of a prescription, so an optimization aimed at
 * one of them can be checked against the rest.
 *
 * Only the varying spaces of a zoom carry a value per configuration. Curvatures,
 * aspheric terms and the fixed spaces are shared, so optimizing any of them
 * against a single configuration is a bet that the others will hold up. Nothing
 * in the merit checks that bet - the solve never evaluates the configurations it
 * was not pointed at - which makes the check an after-the-fact measurement
 * rather than something the builder could enforce.
 *
 * Usage is capture, solve, capture, compare:
 *
 *   auto before = ConfigurationReport::capture(prescription, fields, frequencies);
 *   solver->solve();
 *   auto after = ConfigurationReport::capture(prescription, fields, frequencies);
 *   std::cout << ConfigurationReport::compare(before, after);
 *
 * A single-configuration prescription reports one row, so the same call is
 * harmless where there is nothing to compare across.
 */
class ConfigurationReport {
public:
    /** What one configuration measured at one moment. Java record Configuration. */
    class Configuration {
    public:
        int scenario;
        std::string name;
        double focalLength;
        double fNumber;
        std::vector<int> frequencies;
        /** Null in the Java when the configuration failed to evaluate. */
        std::optional<std::vector<std::vector<double>>> sagittal;
        std::optional<std::vector<std::vector<double>>> tangential;
        std::optional<std::vector<double>> spotRms;
        /** Null in the Java when the configuration evaluated cleanly. */
        std::optional<std::string> failure;
    };

    /** Every configuration of a prescription, measured together. Java record Snapshot. */
    class Snapshot {
    public:
        std::vector<double> fields;
        std::vector<Configuration> configurations;
    };

    /** Applied to each configuration's analysis before it computes. */
    using Configure = std::function<void(Analysis &)>;

    static Snapshot capture(spec::Prescription *prescription,
                            const std::vector<double> &fields,
                            const std::vector<int> &frequencies);

    /**
     * @param configure applied to each configuration's analysis before it
     *                  computes; use it to match whatever the optimization ran
     *                  under, since a report taken on different sampling is not
     *                  comparable with the merit
     */
    static Snapshot capture(spec::Prescription *prescription,
                            const std::vector<double> &fields,
                            const std::vector<int> &frequencies,
                            const Configure &configure);

    /**
     * Side-by-side before and after for every configuration, with a verdict per
     * configuration so the question "did anything get worse" has a one-line
     * answer.
     */
    static std::string compare(const Snapshot &before, const Snapshot &after);

private:
    ConfigurationReport() = delete;

    /** Anything worse than this counts as a regression rather than numerical noise. */
    static constexpr double MTF_REGRESSION = 0.002;
    static constexpr double SPOT_REGRESSION_FRACTION = 0.01;

    static Configuration measure(spec::Prescription *prescription,
                                 const std::vector<double> &fields,
                                 const std::vector<int> &frequencies, int scenario,
                                 const Configure &configure);

    /** One MTF row as before -> after, returning how many values regressed. */
    static int appendRow(std::string &sb, const std::string &label,
                         const std::vector<double> &was, const std::vector<double> &now);

    static void appendAbsolute(std::string &sb, const Configuration &configuration);
};

} // namespace redukti::optim

#endif // REDUKTI_OPTIM_CONFIGURATIONREPORT_H
