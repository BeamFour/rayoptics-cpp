// C++ port of org.redukti.optim.ConfigurationReport
#include "redukti/optim/ConfigurationReport.h"

#include "redukti/Exceptions.h"
#include "redukti/Text.h"
#include "redukti/optim/ParaxHelper.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <typeinfo>

namespace redukti::optim {

namespace {

/**
 * Java's String.format("%<width>.<precision>f"): formatF supplies the digits
 * (Java rounds the shortest round-tripping decimal, not the exact binary value,
 * which is why this cannot be printf) and this pads to the field width.
 */
std::string formatFixed(double value, int width, int precision) {
    std::string s = formatF(value, precision);
    if (static_cast<int>(s.size()) < width)
        s.insert(s.begin(), static_cast<std::size_t>(width - static_cast<int>(s.size())),
                 ' ');
    return s;
}

/** Java's String.format("%-<width>s"). */
std::string padRight(const std::string &s, int width) {
    std::string out = s;
    if (static_cast<int>(out.size()) < width)
        out.append(static_cast<std::size_t>(width - static_cast<int>(out.size())), ' ');
    return out;
}

/**
 * The Java writes line breaks as both "\n" and "%n". On Windows the latter is
 * "\r\n", so the two disagree there; the rest of this port emits bare newlines
 * everywhere and so does this.
 */
constexpr const char *NL = "\n";

const double NaN = std::numeric_limits<double>::quiet_NaN();

/**
 * Java's e.getClass().getSimpleName(). typeid gives an implementation-defined
 * spelling -- "class redukti::IllegalStateException" on MSVC, a mangled name on
 * gcc/clang -- so take the text after the last "::" and drop any leading
 * "class "/"struct ". Diagnostic text only; nothing parses it.
 */
std::string simpleName(const std::exception &e) {
    std::string name = typeid(e).name();
    auto pos = name.rfind("::");
    if (pos != std::string::npos)
        name = name.substr(pos + 2);
    for (const char *prefix : {"class ", "struct "}) {
        std::string p(prefix);
        if (name.compare(0, p.size(), p) == 0)
            name = name.substr(p.size());
    }
    return name;
}

} // namespace

ConfigurationReport::Snapshot ConfigurationReport::capture(
    spec::Prescription *prescription, const std::vector<double> &fields,
    const std::vector<int> &frequencies) {
    return capture(prescription, fields, frequencies, [](Analysis &analysis) {
        analysis.vignetting(spec::VigType::SetPupil)
            .using_gauss_quadrature_pattern(6, 12);
    });
}

ConfigurationReport::Snapshot ConfigurationReport::capture(
    spec::Prescription *prescription, const std::vector<double> &fields,
    const std::vector<int> &frequencies, const Configure &configure) {
    int count = std::max(1, prescription->get_num_configurations());
    Snapshot snapshot;
    snapshot.fields = fields;
    snapshot.configurations.reserve(static_cast<std::size_t>(count));
    for (int scenario = 0; scenario < count; scenario++)
        snapshot.configurations.push_back(
            measure(prescription, fields, frequencies, scenario, configure));
    return snapshot;
}

ConfigurationReport::Configuration ConfigurationReport::measure(
    spec::Prescription *prescription, const std::vector<double> &fields,
    const std::vector<int> &frequencies, int scenario, const Configure &configure) {
    std::string name =
        prescription->_configuration_names.has_value() &&
                scenario < static_cast<int>(prescription->_configuration_names->size())
            ? (*prescription->_configuration_names)[static_cast<std::size_t>(scenario)]
            : "config " + intToString(scenario);
    Configuration configuration;
    configuration.scenario = scenario;
    configuration.name = name;
    configuration.frequencies = frequencies;
    try {
        Analysis analysis(prescription, fields, frequencies, scenario);
        configure(analysis);
        analysis.required_analyses(true, false, true);
        analysis.compute();

        std::vector<std::vector<double>> sagittal(frequencies.size());
        std::vector<std::vector<double>> tangential(frequencies.size());
        for (std::size_t i = 0; i < frequencies.size(); i++) {
            sagittal[i] = (*analysis._mtfs)[i].sag_mtf_by_field;
            tangential[i] = (*analysis._mtfs)[i].tan_mtf_by_field;
        }
        std::vector<double> spot(analysis._spots->size(), 0.0);
        for (std::size_t i = 0; i < spot.size(); i++)
            spot[i] = (*analysis._spots)[i].get_mean_radius();

        configuration.focalLength =
            analysis._pfo[ParaxHelper::Effective_focal_length];
        configuration.fNumber = analysis._pfo[ParaxHelper::Fno];
        configuration.sagittal = std::move(sagittal);
        configuration.tangential = std::move(tangential);
        configuration.spotRms = std::move(spot);
        return configuration;
    } catch (const Exception &e) {
        configuration.focalLength = NaN;
        configuration.fNumber = NaN;
        configuration.sagittal.reset();
        configuration.tangential.reset();
        configuration.spotRms.reset();
        configuration.failure = simpleName(e) + ": " + e.what();
        return configuration;
    }
}

std::string ConfigurationReport::compare(const Snapshot &before, const Snapshot &after) {
    std::string sb;
    sb += "Configuration report";
    sb += NL;
    sb += "fields: ";
    for (double field : after.fields)
        sb += formatF(field, 2) + " ";
    sb += NL;

    for (std::size_t c = 0; c < after.configurations.size(); c++) {
        const auto &now = after.configurations[c];
        const Configuration *was =
            c < before.configurations.size() ? &before.configurations[c] : nullptr;
        sb += NL;
        sb += std::string(78, '-');
        sb += NL;
        sb += "[" + intToString(now.scenario) + "] " + now.name + "   efl " +
              formatF(was == nullptr ? NaN : was->focalLength, 3) + " -> " +
              formatF(now.focalLength, 3) + "   f/# " +
              formatF(was == nullptr ? NaN : was->fNumber, 3) + " -> " +
              formatF(now.fNumber, 3);
        sb += NL;

        if (now.failure.has_value()) {
            sb += "    FAILED: " + *now.failure;
            sb += NL;
            continue;
        }
        if (was == nullptr || was->failure.has_value()) {
            sb += "    no comparable baseline";
            sb += NL;
            appendAbsolute(sb, now);
            continue;
        }

        int regressions = 0;
        double worstMtf = 0.0;
        std::string worstWhere;
        for (std::size_t f = 0; f < now.frequencies.size(); f++) {
            regressions += appendRow(sb, intToString(now.frequencies[f]) + " sag",
                                     (*was->sagittal)[f], (*now.sagittal)[f]);
            regressions += appendRow(sb, intToString(now.frequencies[f]) + " tan",
                                     (*was->tangential)[f], (*now.tangential)[f]);
            for (std::size_t i = 0; i < (*now.sagittal)[f].size(); i++) {
                double dSag = (*was->sagittal)[f][i] - (*now.sagittal)[f][i];
                double dTan = (*was->tangential)[f][i] - (*now.tangential)[f][i];
                if (dSag > worstMtf) {
                    worstMtf = dSag;
                    worstWhere = intToString(now.frequencies[f]) + " sag field " +
                                 intToString(static_cast<int>(i));
                }
                if (dTan > worstMtf) {
                    worstMtf = dTan;
                    worstWhere = intToString(now.frequencies[f]) + " tan field " +
                                 intToString(static_cast<int>(i));
                }
            }
        }

        sb += "    spot RMS  ";
        int spotRegressions = 0;
        for (std::size_t i = 0; i < now.spotRms->size(); i++) {
            double delta = (*now.spotRms)[i] - (*was->spotRms)[i];
            bool worse = delta > (*was->spotRms)[i] * SPOT_REGRESSION_FRACTION;
            if (worse)
                spotRegressions++;
            sb += formatFixed((*now.spotRms)[i], 8, 3) + (worse ? "*" : " ");
        }
        sb += NL;

        if (regressions == 0 && spotRegressions == 0) {
            sb += "    OK - nothing worse than the baseline";
            sb += NL;
        } else {
            sb += "    REGRESSED - " + intToString(regressions) +
                  " MTF value(s) and " + intToString(spotRegressions) +
                  " spot value(s) worse; worst MTF drop " + formatF(worstMtf, 4) +
                  " at " + worstWhere;
            sb += NL;
        }
    }
    return sb;
}

int ConfigurationReport::appendRow(std::string &sb, const std::string &label,
                                   const std::vector<double> &was,
                                   const std::vector<double> &now) {
    int regressions = 0;
    sb += "    " + padRight(label, 8);
    for (std::size_t i = 0; i < now.size(); i++) {
        bool worse = was[i] - now[i] > MTF_REGRESSION;
        if (worse)
            regressions++;
        sb += formatFixed(now[i], 7, 3) + (worse ? "*" : " ");
    }
    sb += NL;
    return regressions;
}

void ConfigurationReport::appendAbsolute(std::string &sb,
                                         const Configuration &configuration) {
    for (std::size_t f = 0; f < configuration.frequencies.size(); f++) {
        sb += "    " + padRight(intToString(configuration.frequencies[f]) + " sag", 8);
        for (double v : (*configuration.sagittal)[f])
            sb += formatFixed(v, 7, 3) + " ";
        sb += NL;
        sb += "    " + padRight(intToString(configuration.frequencies[f]) + " tan", 8);
        for (double v : (*configuration.tangential)[f])
            sb += formatFixed(v, 7, 3) + " ";
        sb += NL;
    }
}

} // namespace redukti::optim
