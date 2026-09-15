// C++ port of org.redukti.optim.OptimizationTrial
#include "redukti/optim/OptimizationTrial.h"

#include "redukti/Exceptions.h"
#include "redukti/Text.h"
#include "redukti/importers/OpticalBenchDataImporter.h"
#include "redukti/optim/OptimizationConfiguration.h"
#include "redukti/optim/OptimizationValidation.h"
#include "redukti/optim/ParaxHelper.h"
#include "redukti/util/Args.h"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <set>
#include <utility>

namespace redukti::optim {

using TrialException = OptimizationTrial::TrialException;

namespace {

std::string lower(const std::string &s) {
    std::string out;
    for (char c : s)
        out.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    return out;
}

bool isSpace(char c) {
    return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f' || c == '\v';
}

/** Java's String.trim(). */
std::string trim(const std::string &s) {
    std::size_t b = 0;
    while (b < s.size() && isSpace(s[b]))
        b++;
    std::size_t e = s.size();
    while (e > b && isSpace(s[e - 1]))
        e--;
    return s.substr(b, e - b);
}

/** Java's `text.split("\\s+")` on text that has already been trimmed. */
std::vector<std::string> words(const std::string &text) {
    std::vector<std::string> result;
    std::size_t i = 0;
    while (i < text.size()) {
        while (i < text.size() && isSpace(text[i]))
            i++;
        std::size_t start = i;
        while (i < text.size() && !isSpace(text[i]))
            i++;
        if (i > start)
            result.push_back(text.substr(start, i - start));
    }
    return result;
}

/** Java's `text.split("\\r?\\n", -1)`: trailing empty lines are kept. */
std::vector<std::string> splitLines(const std::string &text) {
    std::vector<std::string> lines;
    std::string current;
    for (std::size_t i = 0; i < text.size(); i++) {
        if (text[i] == '\n') {
            if (!current.empty() && current.back() == '\r')
                current.pop_back();
            lines.push_back(current);
            current.clear();
        } else
            current.push_back(text[i]);
    }
    lines.push_back(current);
    return lines;
}

/** Java's String.join(" ", ...). */
std::string join(const std::vector<std::string> &parts, std::size_t from, std::size_t to) {
    std::string result;
    for (std::size_t i = from; i < to && i < parts.size(); i++) {
        if (!result.empty())
            result += " ";
        result += parts[i];
    }
    return result;
}

/**
 * The Pattern `\[\s*trial\s+(\d+)\s*\]` and its pipeline twin, matched against a trimmed
 * line. Returns false when the line is some other section; `digitsOut` is set only on a match.
 */
bool matchHeader(const std::string &trimmed, const char *keyword, std::string &digitsOut) {
    std::size_t i = 0;
    if (i >= trimmed.size() || trimmed[i] != '[')
        return false;
    i++;
    while (i < trimmed.size() && isSpace(trimmed[i]))
        i++;
    std::size_t word = std::strlen(keyword);
    if (trimmed.size() < i + word || lower(trimmed.substr(i, word)) != keyword)
        return false;
    i += word;
    if (i >= trimmed.size() || !isSpace(trimmed[i]))
        return false;
    while (i < trimmed.size() && isSpace(trimmed[i]))
        i++;
    std::size_t digits = i;
    while (i < trimmed.size() && std::isdigit(static_cast<unsigned char>(trimmed[i])))
        i++;
    if (i == digits)
        return false;
    std::string text = trimmed.substr(digits, i - digits);
    while (i < trimmed.size() && isSpace(trimmed[i]))
        i++;
    if (i + 1 != trimmed.size() || trimmed[i] != ']')
        return false;
    digitsOut = text;
    return true;
}

/** Java's Integer.parseInt: an optional sign then digits, nothing else. */
bool parseInt(const std::string &text, int &value) {
    if (text.empty())
        return false;
    std::size_t from = (text[0] == '+' || text[0] == '-') ? 1 : 0;
    if (from == text.size())
        return false;
    for (std::size_t i = from; i < text.size(); i++)
        if (!std::isdigit(static_cast<unsigned char>(text[i])))
            return false;
    errno = 0;
    long long parsed = std::strtoll(text.c_str(), nullptr, 10);
    if (errno == ERANGE || parsed < INT_MIN || parsed > INT_MAX)
        return false;
    value = static_cast<int>(parsed);
    return true;
}

/** Java's Double.parseDouble, which accepts the whole string or nothing. */
bool parseDouble(const std::string &text, double &value) {
    if (text.empty())
        return false;
    // strtod also takes hex and "inf"/"nan" spellings Java would too; what it takes and
    // Java does not is leading whitespace, which a word from a split cannot carry anyway.
    const char *start = text.c_str();
    char *end = nullptr;
    errno = 0;
    double parsed = std::strtod(start, &end);
    if (end != start + text.size())
        return false;
    value = parsed;
    return true;
}

/**
 * Enough of java.math.BigDecimal for `fields a to b step s`: a decimal is an integer
 * unscaled value and a scale, so 0.1 is exactly 0.1 and the steps of a field list land on
 * the values that were written rather than drifting the way repeated double addition does.
 */
struct Decimal {
    long long unscaled = 0;
    int scale = 0;

    static bool parse(const std::string &text, Decimal &result) {
        std::size_t i = 0;
        bool negative = false;
        if (i < text.size() && (text[i] == '+' || text[i] == '-')) {
            negative = text[i] == '-';
            i++;
        }
        std::string digits;
        int fraction = 0;
        bool any = false;
        while (i < text.size() && std::isdigit(static_cast<unsigned char>(text[i]))) {
            digits.push_back(text[i++]);
            any = true;
        }
        if (i < text.size() && text[i] == '.') {
            i++;
            while (i < text.size() && std::isdigit(static_cast<unsigned char>(text[i]))) {
                digits.push_back(text[i++]);
                fraction++;
                any = true;
            }
        }
        if (!any)
            return false;
        int exponent = 0;
        if (i < text.size() && (text[i] == 'e' || text[i] == 'E')) {
            i++;
            std::string exponentText;
            if (i < text.size() && (text[i] == '+' || text[i] == '-'))
                exponentText.push_back(text[i++]);
            bool anyExponent = false;
            while (i < text.size() && std::isdigit(static_cast<unsigned char>(text[i]))) {
                exponentText.push_back(text[i++]);
                anyExponent = true;
            }
            if (!anyExponent || !parseInt(exponentText, exponent))
                return false;
        }
        if (i != text.size() || digits.size() > 17)
            return false;
        long long value = 0;
        for (char c : digits)
            value = value * 10 + (c - '0');
        result.unscaled = negative ? -value : value;
        result.scale = fraction - exponent;
        return true;
    }

    static long long rescale(long long unscaled, int by) {
        for (int i = 0; i < by; i++)
            unscaled *= 10;
        return unscaled;
    }

    Decimal add(const Decimal &other) const {
        int target = std::max(scale, other.scale);
        return Decimal{rescale(unscaled, target - scale) +
                           rescale(other.unscaled, target - other.scale),
                       target};
    }

    int compareTo(const Decimal &other) const {
        int target = std::max(scale, other.scale);
        long long a = rescale(unscaled, target - scale);
        long long b = rescale(other.unscaled, target - other.scale);
        return a < b ? -1 : (a > b ? 1 : 0);
    }

    int signum() const { return unscaled < 0 ? -1 : (unscaled > 0 ? 1 : 0); }

    /**
     * BigDecimal.doubleValue(), which rounds the decimal value to the nearest double;
     * strtod on the same digits rounds identically.
     */
    double toDouble() const {
        std::string text = std::to_string(unscaled);
        if (scale > 0) {
            bool negative = !text.empty() && text[0] == '-';
            std::string digits = negative ? text.substr(1) : text;
            while (static_cast<int>(digits.size()) <= scale)
                digits.insert(digits.begin(), '0');
            digits.insert(digits.end() - scale, '.');
            text = (negative ? "-" : "") + digits;
        } else
            for (int i = 0; i < -scale; i++)
                text += "0";
        return std::strtod(text.c_str(), nullptr);
    }
};

/** First-order quantities a paraxial goal can name, in ParaxHelper ids. */
const std::pair<const char *, int> PARAXIAL_QUANTITIES[] = {
    {"efl", ParaxHelper::Effective_focal_length},
    {"bfl", ParaxHelper::Back_focal_length},
    {"ffl", ParaxHelper::Ffl},
    {"fno", ParaxHelper::Fno},
    {"img-dist", ParaxHelper::Image_distance},
    {"obj-dist", ParaxHelper::Object_distance},
    {"pp1", ParaxHelper::Pp1},
    {"ppk", ParaxHelper::Ppk},
    {"enp-dist", ParaxHelper::Enp_dist},
    {"enp-radius", ParaxHelper::Enp_radius},
    {"exp-dist", ParaxHelper::Exp_dist},
    {"exp-radius", ParaxHelper::Exp_radius},
    {"img-ht", ParaxHelper::Img_ht},
    {"obj-ang", ParaxHelper::Obj_ang},
    {"obj-na", ParaxHelper::Obj_na},
    {"img-na", ParaxHelper::Img_na},
    {"red", ParaxHelper::Red},
    {"power", ParaxHelper::Power},
    {"opt-inv", ParaxHelper::Optical_invariant},
};

std::string paraxialQuantityNames() {
    std::string result;
    for (const auto &entry : PARAXIAL_QUANTITIES) {
        if (!result.empty())
            result += ", ";
        result += entry.first;
    }
    return result;
}

TrialException pipelineError(int number, int line, const std::string &message) {
    return TrialException("pipeline " + intToString(number) + ", line " + intToString(line) +
                          ": " + message);
}

std::string defined(const std::string &what, const std::map<int, int> &numbers) {
    if (numbers.empty())
        return "no " + what + "s";
    std::string result = numbers.size() == 1 ? what + " " : what + "s ";
    bool first = true;
    for (const auto &entry : numbers) {
        if (!first)
            result += ", ";
        result += intToString(entry.first);
        first = false;
    }
    return result;
}

// --------------------------------------------------------------------------
// Reading
// --------------------------------------------------------------------------

/** Shared header validation. Map values are one-based source lines for diagnostics. */
struct SectionIndex {
    std::vector<std::string> lines;
    std::map<int, int> trials;
    std::map<int, int> pipelines;

    explicit SectionIndex(const std::string &text) : lines(splitLines(text)) {
        for (std::size_t i = 0; i < lines.size(); i++) {
            std::string header = trim(lines[i]);
            if (header.rfind("[", 0) != 0)
                continue;
            index(header, static_cast<int>(i) + 1, "trial", trials);
            index(header, static_cast<int>(i) + 1, "pipeline", pipelines);
        }
        for (const auto &entry : pipelines) {
            auto both = trials.find(entry.first);
            if (both != trials.end())
                throw TrialException("the number " + intToString(entry.first) +
                                     " is used by both [trial " + intToString(entry.first) +
                                     "] at line " + intToString(both->second) +
                                     " and [pipeline " + intToString(entry.first) +
                                     "] at line " + intToString(entry.second) +
                                     "; trials and pipelines share one numbering");
        }
    }

    static void index(const std::string &header, int line, const char *kind,
                      std::map<int, int> &sections) {
        std::string digits;
        if (!matchHeader(header, kind, digits)) {
            // Java's header.substring(1).stripLeading() starting with the kind.
            std::size_t from = 1;
            while (from < header.size() && isSpace(header[from]))
                from++;
            if (lower(header.substr(from)).rfind(kind, 0) == 0)
                throw TrialException("line " + intToString(line) + ": expected [" + kind +
                                     " <number>], found " + header);
            return;
        }
        int number = 0;
        if (!parseInt(digits, number))
            throw TrialException("line " + intToString(line) + ": " + kind +
                                 " number is out of range: " + digits);
        auto earlier = sections.find(number);
        if (earlier != sections.end())
            throw TrialException("[" + std::string(kind) + " " + intToString(number) +
                                 "] is defined twice, at lines " +
                                 intToString(earlier->second) + " and " + intToString(line));
        sections[number] = line;
    }
};

/** Values given one per field, kept with their line until the field count is known. */
struct PerField {
    std::vector<double> values;
    int line = 0;
};

enum class SelectionKind { All, AllExcept, List };

struct Selection {
    SelectionKind kind;
    std::vector<int> surfaces;
};

/** One aspheric term: the conic constant when index is -1, else a coefficient index. */
struct Term {
    int index;
    std::optional<double> scale;
};

struct AsphericRow {
    int surface;
    std::vector<Term> terms;
    int line;
};

struct ParaxialGoalRow {
    int id;
    double target;
    double weight;
};

struct MtfRows {
    int line = 0;
    std::optional<PerField> sagittal, tangential, weights, sagittalWeights, tangentialWeights;

    explicit MtfRows(int line_) : line(line_) {}
};

/** A LinkedHashMap of per-field rows: insertion order, looked up by frequency. */
class RowsByFrequency {
public:
    PerField *find(int frequency) {
        for (auto &entry : entries)
            if (entry.first == frequency)
                return &entry.second;
        return nullptr;
    }

    const PerField *find(int frequency) const {
        for (const auto &entry : entries)
            if (entry.first == frequency)
                return &entry.second;
        return nullptr;
    }

    void put(int frequency, const PerField &values) {
        if (PerField *existing = find(frequency))
            *existing = values;
        else
            entries.emplace_back(frequency, values);
    }

    bool empty() const { return entries.empty(); }
    const std::vector<std::pair<int, PerField>> &all() const { return entries; }

private:
    std::vector<std::pair<int, PerField>> entries;
};

} // namespace

/** One trial's settings as read, before they are applied to a builder. */
class OptimizationTrial::Reader {
public:
    static Reader parse(const std::string &text, int number);

    int trialNumber() const { return number; }

    /** Apply prescription-dependent checks at each stage, retaining source diagnostics. */
    OptimizationBuilder apply(spec::Prescription *prescription,
                              const OptimizationConfiguration &settings) const;

    /**
     * Resolve shorthand and omitted values once into the same settings the builder uses.
     * The Java's configuration(), renamed because the reader already has a configuration
     * member, the configuration number.
     */
    OptimizationConfiguration toConfiguration() const;

    bool weighted = true;
    bool dLineOnly = false;

private:
    Reader(int number, std::vector<std::string> radii)
        : number(number), radii(std::move(radii)) {}

    void read(int line, const std::string &raw);
    void solver(int line, const std::vector<std::string> &w);
    void vary(int line, const std::vector<std::string> &w);
    Selection selection(int line, const std::vector<std::string> &w, bool curvature);
    std::vector<int> surfaces(int line, const std::vector<std::string> &w, std::size_t from,
                              bool curvature);
    int surface(int line, const std::string &value);
    bool isStop(int surface) const;
    void aspherics(int line, const std::vector<std::string> &w);
    void constrain(int line, const std::vector<std::string> &w);
    void goal(int line, const std::vector<std::string> &w);
    void contrastGoal(int line, const std::vector<std::string> &w);
    void mtfGoal(int line, const std::vector<std::string> &w);
    void spotSampling(int line, const std::vector<std::string> &w);
    void paraxialGoal(int line, const std::vector<std::string> &w);

    /** Checks that need the whole trial: per-field counts, and rows that depend on others. */
    void finish();
    void checkContrastFrequency(int frequency, int line);
    void checkPerField(const std::optional<PerField> &values, const char *what);
    static int firstLine(const std::optional<PerField> &a,
                         const std::optional<PerField> &b = std::nullopt);

    /** One flag per field from 'all', 'all except <field values>' or yes/no per field. */
    std::vector<bool> balanceFlags() const;
    std::vector<double> weightsFor(const PerField *specific,
                                   const std::optional<PerField> &general) const;

    TrialException error(int line, const std::string &message) const;
    void once(int line, const std::string &key);
    void count(int line, const std::vector<std::string> &w, std::size_t expected,
               const std::string &form);
    bool yesNo(int line, const std::vector<std::string> &w, std::size_t index);
    bool yesNoValue(int line, const std::string &value) const;
    double number_(int line, const std::string &value, const std::string &what) const;
    double nonNegative(int line, const std::string &value, const std::string &what);
    double positive(int line, const std::string &value, const std::string &what);
    std::vector<double> numbers(int line, const std::vector<std::string> &w, std::size_t from);
    int integer(int line, const std::string &value, const std::string &what);
    int nonNegativeInt(int line, const std::string &value, const std::string &what);
    int positiveInt(int line, const std::string &value, const std::string &what);
    std::vector<int> positiveInts(int line, const std::vector<std::string> &w,
                                  std::size_t from, const std::string &what);
    std::vector<double> fieldList(int line, const std::vector<std::string> &w);
    Decimal decimal(int line, const std::string &value);

    int number;
    /** Radius column of each [lens data] row: surface numbers are positions in it. */
    std::vector<std::string> radii;
    std::set<std::string> seen;

    std::optional<std::string> description;
    std::optional<std::string> outdir;
    int configuration = 0;
    std::optional<std::vector<double>> fields;
    std::optional<std::vector<int>> frequencies;
    std::optional<spec::VigType> vignetting;
    bool freezeVignetting = false;
    std::optional<bool> checkSpotApertures;
    std::optional<double> solverFtol;
    std::optional<double> solverXtol;
    std::optional<double> solverGtol;
    std::optional<int> solverMaxEvaluations;

    std::optional<Selection> curvatures;
    std::optional<Selection> thicknesses;
    bool existingAspherics = false;
    std::vector<AsphericRow> asphericRows;

    std::optional<double> curvatureConstraint;
    int curvatureConstraintLine = 0;
    std::optional<double> thicknessConstraint;
    std::optional<double> edgeConstraint;

    std::optional<std::vector<int>> contrastFrequencies;
    std::optional<PerField> contrastSagittal;
    std::optional<PerField> contrastTangential;
    RowsByFrequency contrastSagittalFor;
    RowsByFrequency contrastTangentialFor;
    std::optional<std::vector<std::string>> balanceFields;
    int balanceLine = 0;
    double balanceWeight = OptimizationBuilder::NOMINAL_BALANCE_WEIGHT;
    std::optional<std::vector<int>> contrastSampling;
    std::optional<bool> calibrateContrast;
    std::optional<bool> exitPupilAiming;
    std::optional<bool> centerContrast;
    int contrastSettingsLine = 0;

    /** A LinkedHashMap in the Java: MTF goals keep the order their rows were given. */
    std::vector<std::pair<int, MtfRows>> mtf;
    MtfRows &mtfRowsFor(int frequency, int line);

    std::optional<PerField> spotRms;
    std::optional<PerField> spotRmsWeights;
    std::optional<PerField> spotMaxRadius;
    std::optional<PerField> spotMaxRadiusWeights;
    std::optional<PerField> spotDeviation;
    std::optional<PerField> spotDeviationX;
    std::optional<PerField> spotDeviationY;
    std::optional<std::vector<int>> gaussianSampling;
    double gaussianInnerRadius = 0.0;
    std::optional<int> hexapolarRays;
    int hexapolarLine = 0;

    bool rayAberrations = false;
    std::vector<ParaxialGoalRow> paraxialGoals;
};

OptimizationTrial::Reader OptimizationTrial::Reader::parse(const std::string &text, int number) {
    SectionIndex index(text);
    const std::vector<std::string> &lines = index.lines;

    // The lens data as the prescription reader sees it.
    std::vector<std::string> radii;
    std::string section;
    for (const std::string &raw : lines) {
        std::vector<std::string> w = OptimizationTrial::splitTabs(raw);
        if (w.empty() || w[0].rfind("#", 0) == 0)
            continue;
        if (w[0].rfind("[", 0) == 0) {
            section = w[0];
            continue;
        }
        if (section == "[lens data]" && w.size() >= 2)
            radii.push_back(w[1]);
    }

    const auto &headers = index.trials;
    if (headers.find(number) == headers.end())
        throw TrialException("there is no [trial " + intToString(number) +
                             "] in this prescription; " +
                             (headers.empty() ? "it defines no trials"
                                              : "it defines " + defined("trial", headers)));

    Reader reader(number, radii);
    // The stored one-based header line is the zero-based first body line.
    for (std::size_t i = static_cast<std::size_t>(headers.at(number));
         i < lines.size() && trim(lines[i]).rfind("[", 0) != 0; i++)
        reader.read(static_cast<int>(i) + 1, lines[i]);
    reader.finish();
    return reader;
}

void OptimizationTrial::Reader::read(int line, const std::string &raw) {
    std::size_t hash = raw.find('#');
    std::string text = trim(hash != std::string::npos ? raw.substr(0, hash) : raw);
    if (text.empty())
        return;
    std::vector<std::string> w = words(text);
    std::string keyword = lower(w[0]);
    if (keyword == "description") {
        once(line, "description");
        description = trim(text.substr(w[0].size()));
    } else if (keyword == "outdir") {
        once(line, "outdir");
        if (w.size() < 2)
            throw error(line, "expected 'outdir <directory>'");
        outdir = trim(text.substr(w[0].size()));
    } else if (keyword == "configuration") {
        once(line, "configuration");
        count(line, w, 2, "configuration <number>");
        configuration = nonNegativeInt(line, w[1], "configuration");
    } else if (keyword == "fields") {
        once(line, "fields");
        fields = fieldList(line, w);
    } else if (keyword == "frequencies") {
        once(line, "frequencies");
        frequencies = positiveInts(line, w, 1, "frequencies");
    } else if (keyword == "weighted") {
        once(line, "weighted");
        weighted = yesNo(line, w, 1);
    } else if (keyword == "d-line-only") {
        once(line, "d-line-only");
        dLineOnly = yesNo(line, w, 1);
    } else if (keyword == "check-spot-apertures") {
        once(line, "check-spot-apertures");
        checkSpotApertures = yesNo(line, w, 1);
    } else if (keyword == "vignetting") {
        once(line, "vignetting");
        if (w.size() < 2 || w.size() > 3)
            throw error(line, "expected 'vignetting <type> [frozen]'");
        try {
            vignetting = util::Args::parse_vig_type(w[1]);
        } catch (const IllegalArgumentException &) {
            throw error(line, "unknown vignetting type '" + w[1] +
                                  "', expected one of: " + util::Args::vig_type_names());
        }
        if (w.size() == 3) {
            if (lower(w[2]) != "frozen")
                throw error(line, "expected 'frozen' after the vignetting type, found '" +
                                      w[2] + "'");
            freezeVignetting = true;
        }
    } else if (keyword == "solver")
        solver(line, w);
    else if (keyword == "vary")
        vary(line, w);
    else if (keyword == "constrain")
        constrain(line, w);
    else if (keyword == "goal")
        goal(line, w);
    else
        throw error(line, "unknown keyword '" + w[0] + "'");
}

/**
 * The lmder stopping tolerances. Left alone, the solver keeps its own defaults;
 * see the note in LMDerSolver on why xtol is off by default.
 */
void OptimizationTrial::Reader::solver(int line, const std::vector<std::string> &w) {
    if (w.size() != 3)
        throw error(line, "expected 'solver ftol|xtol|gtol|max-evaluations <value>'");
    std::string what = lower(w[1]);
    if (what == "ftol") {
        once(line, "solver ftol");
        solverFtol = nonNegative(line, w[2], "ftol");
    } else if (what == "xtol") {
        once(line, "solver xtol");
        solverXtol = nonNegative(line, w[2], "xtol");
    } else if (what == "gtol") {
        once(line, "solver gtol");
        solverGtol = nonNegative(line, w[2], "gtol");
    } else if (what == "max-evaluations") {
        once(line, "solver max-evaluations");
        solverMaxEvaluations = positiveInt(line, w[2], "max-evaluations");
    } else
        throw error(line, "unknown solver setting '" + w[1] +
                              "'; expected ftol, xtol, gtol or max-evaluations");
}

void OptimizationTrial::Reader::vary(int line, const std::vector<std::string> &w) {
    if (w.size() < 3)
        throw error(line, "expected 'vary curvatures|thicknesses|aspherics ...'");
    std::string what = lower(w[1]);
    if (what == "curvatures") {
        once(line, "vary curvatures");
        curvatures = selection(line, w, true);
    } else if (what == "thicknesses") {
        once(line, "vary thicknesses");
        thicknesses = selection(line, w, false);
    } else if (what == "aspherics")
        aspherics(line, w);
    else
        throw error(line, "cannot vary '" + w[1] +
                              "'; expected curvatures, thicknesses or aspherics");
}

Selection OptimizationTrial::Reader::selection(int line, const std::vector<std::string> &w, bool curvature) {
    if (lower(w[2]) == "all") {
        if (w.size() == 3)
            return Selection{SelectionKind::All, {}};
        if (w.size() > 4 && lower(w[3]) == "except")
            return Selection{SelectionKind::AllExcept, surfaces(line, w, 4, false)};
        throw error(line, "expected 'all' or 'all except <surfaces>'");
    }
    return Selection{SelectionKind::List, surfaces(line, w, 2, curvature)};
}

std::vector<int> OptimizationTrial::Reader::surfaces(int line, const std::vector<std::string> &w,
                                  std::size_t from, bool curvature) {
    std::vector<int> result;
    std::set<int> unique;
    for (std::size_t i = from; i < w.size(); i++) {
        int s = surface(line, w[i]);
        if (!unique.insert(s).second)
            throw error(line, "surface " + intToString(s) + " is listed twice");
        if (curvature && isStop(s))
            throw error(line, "surface " + intToString(s) +
                                  " is a stop; it has no curvature to vary");
        result.push_back(s);
    }
    return result;
}

int OptimizationTrial::Reader::surface(int line, const std::string &value) {
    int s = 0;
    if (!parseInt(value, s))
        throw error(line, "expected a surface number, found '" + value +
                              "'; surfaces are numbered by their position in [lens data], "
                              "from 0");
    if (s < 0 || s >= static_cast<int>(radii.size()))
        throw error(line, "there is no surface " + value +
                              "; [lens data] has surfaces 0 to " +
                              intToString(static_cast<int>(radii.size()) - 1));
    return s;
}

bool OptimizationTrial::Reader::isStop(int surface) const {
    const std::string &radius = radii[static_cast<std::size_t>(surface)];
    return radius == "AS" || radius == "FS";
}

void OptimizationTrial::Reader::aspherics(int line, const std::vector<std::string> &w) {
    if (w.size() == 3 && lower(w[2]) == "existing") {
        once(line, "vary aspherics existing");
        existingAspherics = true;
        return;
    }
    if (w.size() < 4)
        throw error(line,
                    "expected 'vary aspherics existing' or 'vary aspherics <surface> <terms>'");
    int s = surface(line, w[2]);
    once(line, "vary aspherics " + intToString(s));
    std::vector<Term> terms;
    std::set<int> indices;
    for (std::size_t i = 3; i < w.size(); i++) {
        const std::string &token = w[i];
        std::size_t colon = token.find(':');
        std::string name = colon != std::string::npos ? token.substr(0, colon) : token;
        std::optional<double> scale;
        if (colon != std::string::npos)
            scale = positive(line, token.substr(colon + 1), "scale of " + name);
        int index = 0;
        if (lower(name) == "k") {
            if (scale.has_value())
                throw error(line, "K takes no scale");
            index = -1;
        } else {
            if (!parseInt(name, index))
                throw error(line, "expected K or a coefficient index, found '" + name + "'");
            if (index < 0)
                throw error(line, "coefficient indices start at 0, found " + name);
        }
        if (!indices.insert(index).second)
            throw error(line, (index < 0 ? std::string("K")
                                         : "coefficient " + intToString(index)) +
                                  " is listed twice");
        terms.push_back(Term{index, scale});
    }
    asphericRows.push_back(AsphericRow{s, terms, line});
}

void OptimizationTrial::Reader::constrain(int line, const std::vector<std::string> &w) {
    if (w.size() < 2 || w.size() > 3)
        throw error(line, "expected 'constrain curvatures|thicknesses|edges [weight]'");
    double weight = w.size() == 3 ? nonNegative(line, w[2], "constraint weight") : 1.0;
    std::string what = lower(w[1]);
    if (what == "curvatures") {
        once(line, "constrain curvatures");
        curvatureConstraint = weight;
        curvatureConstraintLine = line;
    } else if (what == "thicknesses") {
        once(line, "constrain thicknesses");
        thicknessConstraint = weight;
    } else if (what == "edges") {
        once(line, "constrain edges");
        edgeConstraint = weight;
    } else
        throw error(line, "cannot constrain '" + w[1] +
                              "'; expected curvatures, thicknesses or edges");
}

void OptimizationTrial::Reader::goal(int line, const std::vector<std::string> &w) {
    if (w.size() < 3)
        throw error(line, "expected 'goal <type> ...'");
    std::string what = lower(w[1]);
    if (what == "contrast")
        contrastGoal(line, w);
    else if (what == "mtf")
        mtfGoal(line, w);
    else if (what == "spot-rms") {
        if (lower(w[2]) == "weights") {
            once(line, "goal spot-rms weights");
            spotRmsWeights = PerField{numbers(line, w, 3), line};
        } else {
            once(line, "goal spot-rms");
            spotRms = PerField{numbers(line, w, 2), line};
        }
    } else if (what == "spot-max-radius") {
        if (lower(w[2]) == "weights") {
            once(line, "goal spot-max-radius weights");
            spotMaxRadiusWeights = PerField{numbers(line, w, 3), line};
        } else {
            once(line, "goal spot-max-radius");
            spotMaxRadius = PerField{numbers(line, w, 2), line};
        }
    } else if (what == "spot-deviation") {
        std::string direction = lower(w[2]);
        if (direction == "x") {
            once(line, "goal spot-deviation x");
            spotDeviationX = PerField{numbers(line, w, 3), line};
        } else if (direction == "y") {
            once(line, "goal spot-deviation y");
            spotDeviationY = PerField{numbers(line, w, 3), line};
        } else {
            once(line, "goal spot-deviation");
            spotDeviation = PerField{numbers(line, w, 2), line};
        }
    } else if (what == "spot")
        spotSampling(line, w);
    else if (what == "ray-aberrations") {
        once(line, "goal ray-aberrations");
        rayAberrations = yesNo(line, w, 2);
    } else if (what == "paraxial")
        paraxialGoal(line, w);
    else
        throw error(line, "unknown goal '" + w[1] +
                              "'; expected contrast, mtf, spot-rms, spot-max-radius, "
                              "spot-deviation, spot, ray-aberrations or paraxial");
}

void OptimizationTrial::Reader::contrastGoal(int line, const std::vector<std::string> &w) {
    std::string what = lower(w[2]);
    if (what == "sag") {
        once(line, "goal contrast sag");
        contrastSagittal = PerField{numbers(line, w, 3), line};
    } else if (what == "tan") {
        once(line, "goal contrast tan");
        contrastTangential = PerField{numbers(line, w, 3), line};
    } else if (what == "balance") {
        once(line, "goal contrast balance");
        std::vector<std::string> tokens(w.begin() + 3, w.end());
        std::size_t n = tokens.size();
        if (n >= 2 && lower(tokens[n - 2]) == "weight") {
            balanceWeight = nonNegative(line, tokens[n - 1], "balance weight");
            tokens.resize(n - 2);
        }
        if (tokens.empty())
            throw error(line, "expected the fields to balance: all, all except <fields>, or "
                              "yes/no per field");
        balanceFields = tokens;
        balanceLine = line;
    } else if (what == "sampling") {
        once(line, "goal contrast sampling");
        count(line, w, 5, "goal contrast sampling <rings> <spokes>");
        contrastSampling = std::vector<int>{positiveInt(line, w[3], "rings"),
                                            positiveInt(line, w[4], "spokes")};
        if ((*contrastSampling)[1] < 3)
            throw error(line, "contrast sampling requires at least 1 ring and 3 spokes");
        contrastSettingsLine = line;
    } else if (what == "calibrate") {
        once(line, "goal contrast calibrate");
        calibrateContrast = yesNo(line, w, 3);
        contrastSettingsLine = line;
    } else if (what == "exit-pupil-aiming") {
        once(line, "goal contrast exit-pupil-aiming");
        exitPupilAiming = yesNo(line, w, 3);
        contrastSettingsLine = line;
    } else if (what == "centering") {
        once(line, "goal contrast centering");
        centerContrast = yesNo(line, w, 3);
        contrastSettingsLine = line;
    } else if (w.size() >= 4 && (lower(w[3]) == "sag" || lower(w[3]) == "tan")) {
        int frequency = positiveInt(line, w[2], "contrast frequency");
        bool sagittal = lower(w[3]) == "sag";
        once(line, "goal contrast " + intToString(frequency) + " " + lower(w[3]));
        PerField values{numbers(line, w, 4), line};
        (sagittal ? contrastSagittalFor : contrastTangentialFor).put(frequency, values);
    } else {
        once(line, "goal contrast");
        contrastFrequencies = positiveInts(line, w, 2, "contrast frequencies");
    }
}

MtfRows &OptimizationTrial::Reader::mtfRowsFor(int frequency, int line) {
    for (auto &entry : mtf)
        if (entry.first == frequency)
            return entry.second;
    mtf.emplace_back(frequency, MtfRows(line));
    return mtf.back().second;
}

void OptimizationTrial::Reader::mtfGoal(int line, const std::vector<std::string> &w) {
    if (w.size() < 5)
        throw error(line, "expected 'goal mtf <frequency> sag|tan [weights] <values>' or "
                          "'goal mtf <frequency> weights <values>'");
    int frequency = positiveInt(line, w[2], "MTF frequency");
    MtfRows &rows = mtfRowsFor(frequency, line);
    std::string what = lower(w[3]);
    if (what == "weights") {
        once(line, "goal mtf " + intToString(frequency) + " weights");
        rows.weights = PerField{numbers(line, w, 4), line};
    } else if (what == "sag" || what == "tan") {
        bool sagittal = what == "sag";
        if (lower(w[4]) == "weights") {
            once(line, "goal mtf " + intToString(frequency) + " " + what + " weights");
            PerField values{numbers(line, w, 5), line};
            if (sagittal)
                rows.sagittalWeights = values;
            else
                rows.tangentialWeights = values;
        } else {
            once(line, "goal mtf " + intToString(frequency) + " " + what);
            PerField values{numbers(line, w, 4), line};
            if (sagittal)
                rows.sagittal = values;
            else
                rows.tangential = values;
        }
    } else
        throw error(line, "expected sag, tan or weights after the MTF frequency, found '" +
                              w[3] + "'");
}

void OptimizationTrial::Reader::spotSampling(int line, const std::vector<std::string> &w) {
    if (w.size() < 4 || lower(w[2]) != "sampling")
        throw error(line, "expected 'goal spot sampling gaussian <rings> <spokes> [<inner "
                          "radius>]' or 'goal spot sampling hexapolar <rings>'");
    std::string what = lower(w[3]);
    if (what == "gaussian") {
        once(line, "goal spot sampling gaussian");
        if (w.size() != 6 && w.size() != 7)
            throw error(line, "expected 'goal spot sampling gaussian <rings> <spokes> "
                              "[<inner radius>]'");
        gaussianSampling = std::vector<int>{positiveInt(line, w[4], "rings"),
                                            positiveInt(line, w[5], "spokes")};
        gaussianInnerRadius = w.size() == 7 ? nonNegative(line, w[6], "inner radius") : 0.0;
    } else if (what == "hexapolar") {
        once(line, "goal spot sampling hexapolar");
        count(line, w, 5, "goal spot sampling hexapolar <rings>");
        hexapolarRays = positiveInt(line, w[4], "rings");
        hexapolarLine = line;
    } else
        throw error(line,
                    "unknown spot sampling '" + w[3] + "'; expected gaussian or hexapolar");
}

void OptimizationTrial::Reader::paraxialGoal(int line, const std::vector<std::string> &w) {
    if (w.size() != 4 && w.size() != 6)
        throw error(line, "expected 'goal paraxial <quantity> <target> [weight <w>]'");
    std::string quantity = lower(w[2]);
    std::optional<int> id;
    for (const auto &entry : PARAXIAL_QUANTITIES)
        if (quantity == entry.first)
            id = entry.second;
    if (!id.has_value())
        throw error(line, "unknown paraxial quantity '" + w[2] +
                              "'; expected one of: " + paraxialQuantityNames());
    once(line, "goal paraxial " + quantity);
    double target = number_(line, w[3], "paraxial target");
    double weight = 1.0;
    if (w.size() == 6) {
        if (lower(w[4]) != "weight")
            throw error(line, "expected 'weight' after the target, found '" + w[4] + "'");
        weight = nonNegative(line, w[5], "paraxial weight");
    }
    paraxialGoals.push_back(ParaxialGoalRow{*id, target, weight});
}

void OptimizationTrial::Reader::finish() {
    if (!fields.has_value())
        throw TrialException("trial " + intToString(number) + ": 'fields' is required");
    if (!frequencies.has_value())
        throw TrialException("trial " + intToString(number) + ": 'frequencies' is required");

    checkPerField(contrastSagittal, "weights");
    checkPerField(contrastTangential, "weights");
    for (const auto &entry : contrastSagittalFor.all())
        checkPerField(entry.second, "weights");
    for (const auto &entry : contrastTangentialFor.all())
        checkPerField(entry.second, "weights");
    if (!contrastFrequencies.has_value()) {
        int line = firstLine(contrastSagittal, contrastTangential);
        if (line == 0 && !contrastSagittalFor.empty())
            line = contrastSagittalFor.all().front().second.line;
        if (line == 0 && !contrastTangentialFor.empty())
            line = contrastTangentialFor.all().front().second.line;
        if (line == 0 && balanceFields.has_value())
            line = balanceLine;
        if (line == 0)
            line = contrastSettingsLine;
        if (line != 0)
            throw error(line, "contrast settings need a 'goal contrast <frequencies>' line");
    } else {
        for (const auto &entry : contrastSagittalFor.all())
            checkContrastFrequency(entry.first, entry.second.line);
        for (const auto &entry : contrastTangentialFor.all())
            checkContrastFrequency(entry.first, entry.second.line);
    }

    for (auto &entry : mtf) {
        int frequency = entry.first;
        MtfRows &rows = entry.second;
        if (!rows.sagittal.has_value() || !rows.tangential.has_value())
            throw error(rows.line, "MTF goals at " + intToString(frequency) +
                                       " need both a sag and a tan row of targets");
        if (std::find(frequencies->begin(), frequencies->end(), frequency) ==
            frequencies->end())
            throw error(rows.line, "MTF goal frequency " + intToString(frequency) +
                                       " is not one of 'frequencies'");
        checkPerField(rows.sagittal, "targets");
        checkPerField(rows.tangential, "targets");
        for (const std::optional<PerField> *targets : {&rows.sagittal, &rows.tangential})
            for (double target : (*targets)->values)
                if (target < 0.0 || target > 100.0)
                    throw error((*targets)->line,
                                "MTF targets are percentages, between 0 and 100");
        checkPerField(rows.weights, "weights");
        checkPerField(rows.sagittalWeights, "weights");
        checkPerField(rows.tangentialWeights, "weights");
    }

    checkPerField(spotRms, "targets");
    checkPerField(spotRmsWeights, "weights");
    checkPerField(spotMaxRadius, "targets");
    checkPerField(spotMaxRadiusWeights, "weights");
    for (const std::optional<PerField> *targets : {&spotRms, &spotMaxRadius})
        if (targets->has_value())
            OptimizationValidation::range(
                (*targets)->values, std::numeric_limits<double>::infinity(), [&] {
                    return error((*targets)->line, "spot targets must be finite and non-negative");
                });
    if (spotRmsWeights.has_value() && !spotRms.has_value())
        throw error(spotRmsWeights->line,
                    "spot-rms weights need a 'goal spot-rms <targets>' line");
    if (spotMaxRadiusWeights.has_value() && !spotMaxRadius.has_value())
        throw error(spotMaxRadiusWeights->line,
                    "spot-max-radius weights need a 'goal spot-max-radius <targets>' line");

    checkPerField(spotDeviation, "weights");
    checkPerField(spotDeviationX, "weights");
    checkPerField(spotDeviationY, "weights");
    if (spotDeviation.has_value() &&
        (spotDeviationX.has_value() || spotDeviationY.has_value()))
        throw error(firstLine(spotDeviationX, spotDeviationY),
                    "give spot deviation weights either as one row or as x and y rows, not "
                    "both");
    if (spotDeviationX.has_value() != spotDeviationY.has_value())
        throw error(firstLine(spotDeviationX, spotDeviationY),
                    "spot deviation needs both an x and a y row");
    if ((spotDeviation.has_value() || spotDeviationX.has_value()) &&
        (hexapolarRays.has_value() || spotMaxRadius.has_value())) {
        int conflictLine =
            std::max(hexapolarLine, spotMaxRadius.has_value() ? spotMaxRadius->line : 0);
        for (const std::optional<PerField> *row : {&spotDeviation, &spotDeviationX, &spotDeviationY})
            if (row->has_value())
                conflictLine = std::max(conflictLine, (*row)->line);
        throw error(conflictLine, "spot deviation goals require Gaussian-quadrature spot sampling");
    }
}

void OptimizationTrial::Reader::checkContrastFrequency(int frequency, int line) {
    if (std::find(contrastFrequencies->begin(), contrastFrequencies->end(), frequency) ==
        contrastFrequencies->end())
        throw error(line, "contrast frequency " + intToString(frequency) +
                              " is not one of the 'goal contrast' frequencies");
}

void OptimizationTrial::Reader::checkPerField(const std::optional<PerField> &values, const char *what) {
    if (!values.has_value())
        return;
    if (values->values.size() != fields->size())
        throw error(values->line,
                    "expected " + intToString(static_cast<int>(fields->size())) + " " + what +
                        ", one per field, but found " +
                        intToString(static_cast<int>(values->values.size())));
    if (std::string(what) == "weights")
        for (double value : values->values)
            if (value < 0.0)
                throw error(values->line, "weights must not be negative");
}

int OptimizationTrial::Reader::firstLine(const std::optional<PerField> &a, const std::optional<PerField> &b) {
    if (a.has_value())
        return a->line;
    if (b.has_value())
        return b->line;
    return 0;
}

std::vector<bool> OptimizationTrial::Reader::balanceFlags() const {
    std::vector<bool> flags(fields->size(), false);
    std::string first = lower((*balanceFields)[0]);
    if (first == "all") {
        std::fill(flags.begin(), flags.end(), true);
        if (balanceFields->size() == 1)
            return flags;
        if (balanceFields->size() < 3 || lower((*balanceFields)[1]) != "except")
            throw error(balanceLine, "expected 'all' or 'all except <field values>'");
        for (std::size_t i = 2; i < balanceFields->size(); i++) {
            const std::string &value = (*balanceFields)[i];
            double field = number_(balanceLine, value, "field");
            int index = -1;
            for (std::size_t f = 0; f < fields->size(); f++)
                if ((*fields)[f] == field)
                    index = static_cast<int>(f);
            if (index < 0)
                throw error(balanceLine, "there is no field " + value + " in 'fields'");
            flags[static_cast<std::size_t>(index)] = false;
        }
        return flags;
    }
    if (balanceFields->size() != fields->size())
        throw error(balanceLine, "expected 'all', 'all except <field values>' or " +
                                     intToString(static_cast<int>(fields->size())) +
                                     " yes/no values, one per field");
    for (std::size_t i = 0; i < fields->size(); i++)
        flags[i] = yesNoValue(balanceLine, (*balanceFields)[i]);
    return flags;
}

OptimizationBuilder OptimizationTrial::Reader::apply(spec::Prescription *prescription,
                                  const OptimizationConfiguration &settings) const {
    const auto &surfaceList = prescription->_surface_list;
    if (surfaceList.size() != radii.size())
        throw IllegalArgumentException(
            "the prescription has " + intToString(static_cast<int>(surfaceList.size())) +
            " surfaces but the trial's [lens data] has " +
            intToString(static_cast<int>(radii.size())) +
            "; build the prescription from the same text as the trial");
    if (curvatureConstraint.has_value() && curvatures.has_value() &&
        curvatures->kind == SelectionKind::List) {
        for (int s : curvatures->surfaces)
            if (surfaceList[static_cast<std::size_t>(s)]._radius == 0.0)
                throw error(curvatureConstraintLine,
                            "cannot constrain curvature of flat surface " + intToString(s) +
                                "; a fractional curvature constraint needs a non-zero "
                                "starting curvature; remove this surface from 'vary "
                                "curvatures' or omit 'constrain curvatures'");
    }
    // Add explicit terms through the public API to keep its prescription-dependent
    // validation (asphere kind and coefficient scaling) and source-line errors.
    OptimizationConfiguration stage = settings.copy();
    stage.asphericTerms.clear();
    OptimizationBuilder builder(prescription, stage);
    for (const AsphericRow &row : asphericRows) {
        for (const Term &term : row.terms) {
            try {
                if (term.index < 0)
                    builder.varyConic(row.surface);
                else if (term.scale.has_value())
                    builder.varyAsphericCoefficient(row.surface, term.index, *term.scale);
                else
                    builder.varyAsphericCoefficient(row.surface, term.index);
            } catch (const IllegalArgumentException &e) {
                throw error(row.line, e.getMessage());
            }
        }
    }
    return builder;
}

OptimizationConfiguration OptimizationTrial::Reader::toConfiguration() const {
    OptimizationConfiguration c;
    c.description = description;
    c.outdir = outdir;
    c.fields = *fields;
    c.mtfFrequencies = *frequencies;
    c.scenario = configuration;
    c.weighted = weighted;
    c.dLineOnly = dLineOnly;
    if (vignetting.has_value())
        c.vigType = *vignetting;
    c.freezeVignetting = freezeVignetting;
    c.solverTolerances = SolverTolerances(
        solverFtol.value_or(SolverTolerances::defaultFtol()),
        solverXtol.value_or(SolverTolerances::DEFAULT_XTOL),
        solverGtol.value_or(SolverTolerances::defaultGtol()),
        solverMaxEvaluations.value_or(SolverTolerances::FROM_VARIABLE_COUNT));
    if (checkSpotApertures.has_value())
        c.checkSpotApertures = *checkSpotApertures;
    if (curvatures.has_value()) {
        c.allCurvatureSurfaces = curvatures->kind != SelectionKind::List;
        if (curvatures->kind == SelectionKind::List)
            c.curvatureSurfaces = curvatures->surfaces;
        if (curvatures->kind == SelectionKind::AllExcept)
            c.curvatureExclusions = curvatures->surfaces;
    }
    if (thicknesses.has_value()) {
        c.allThicknessSurfaces = thicknesses->kind != SelectionKind::List;
        if (thicknesses->kind == SelectionKind::List)
            c.thicknessSurfaces = thicknesses->surfaces;
        if (thicknesses->kind == SelectionKind::AllExcept)
            c.thicknessExclusions = thicknesses->surfaces;
    }
    c.includeExistingAspherics = existingAspherics;
    for (const AsphericRow &row : asphericRows)
        for (const Term &term : row.terms)
            c.asphericTerms.push_back(
                OptimizationBuilder::AsphericTerm{row.surface, term.index, term.scale});
    c.curvatureConstraintWeight = curvatureConstraint;
    c.thicknessConstraintWeight = thicknessConstraint;
    c.edgeThicknessConstraintWeight = edgeConstraint;
    if (contrastFrequencies.has_value()) {
        for (int frequency : *contrastFrequencies)
            c.contrastGoals.push_back(OptimizationBuilder::contrast(
                frequency, weightsFor(contrastSagittalFor.find(frequency), contrastSagittal),
                weightsFor(contrastTangentialFor.find(frequency), contrastTangential)));
        if (contrastSampling.has_value()) {
            c.contrastRings = (*contrastSampling)[0];
            c.contrastSpokes = (*contrastSampling)[1];
        }
        if (calibrateContrast.has_value())
            c.calibrateContrastFrequency = *calibrateContrast;
        if (exitPupilAiming.has_value())
            c.aimContrastAtExitPupil = *exitPupilAiming;
        if (centerContrast.has_value())
            c.centerContrastResiduals = *centerContrast;
        if (balanceFields.has_value()) {
            c.contrastBalanceFields = balanceFlags();
            c.contrastBalanceWeight = balanceWeight;
        }
    }
    for (const auto &entry : mtf) {
        const MtfRows &rows = entry.second;
        const std::optional<PerField> &both = rows.weights;
        const std::optional<PerField> &sagittalWeights =
            rows.sagittalWeights.has_value() ? rows.sagittalWeights : both;
        const std::optional<PerField> &tangentialWeights =
            rows.tangentialWeights.has_value() ? rows.tangentialWeights : both;
        std::vector<double> ones(fields->size(), 1.0);
        c.mtfGoals.push_back(OptimizationBuilder::mtf(
            entry.first, rows.sagittal->values, rows.tangential->values,
            sagittalWeights.has_value() ? sagittalWeights->values : ones,
            tangentialWeights.has_value() ? tangentialWeights->values : ones));
    }
    if (spotRms.has_value())
        c.spotRmsGoals = OptimizationBuilder::SpotGoals{
            spotRms->values, spotRmsWeights.has_value()
                                 ? spotRmsWeights->values
                                 : std::vector<double>(spotRms->values.size(), 1.0)};
    if (spotMaxRadius.has_value())
        c.spotMaxRadiusGoals = OptimizationBuilder::SpotGoals{
            spotMaxRadius->values, spotMaxRadiusWeights.has_value()
                                       ? spotMaxRadiusWeights->values
                                       : std::vector<double>(spotMaxRadius->values.size(), 1.0)};
    if (spotDeviation.has_value() || spotDeviationX.has_value()) {
        c.addSpotDeviationGoals = true;
        c.spotDeviationXWeights = (spotDeviation.has_value() ? spotDeviation : spotDeviationX)->values;
        c.spotDeviationYWeights = (spotDeviation.has_value() ? spotDeviation : spotDeviationY)->values;
    }
    if (gaussianSampling.has_value())
        c.gaussianSampling((*gaussianSampling)[0], (*gaussianSampling)[1], gaussianInnerRadius);
    if (hexapolarRays.has_value()) {
        c.useHexapolarSpotPattern = true;
        c.hexapolarSpotRays = *hexapolarRays;
    }
    c.addRayAberrationGoals = rayAberrations;
    for (const ParaxialGoalRow &goal : paraxialGoals)
        c.paraxialGoals.push_back(
            OptimizationBuilder::ParaxialGoal{goal.id, goal.target, goal.weight});
    return c;
}

std::vector<double> OptimizationTrial::Reader::weightsFor(const PerField *specific,
                                       const std::optional<PerField> &general) const {
    if (specific != nullptr)
        return specific->values;
    if (general.has_value())
        return general->values;
    return std::vector<double>(fields->size(), 1.0);
}

// --------------------------------------------------------------------------
// Helpers
// --------------------------------------------------------------------------

TrialException OptimizationTrial::Reader::error(int line, const std::string &message) const {
    return TrialException("trial " + intToString(number) + ", line " + intToString(line) +
                          ": " + message);
}

void OptimizationTrial::Reader::once(int line, const std::string &key) {
    if (!seen.insert(key).second)
        throw error(line, "'" + key + "' is given more than once");
}

void OptimizationTrial::Reader::count(int line, const std::vector<std::string> &w, std::size_t expected,
                   const std::string &form) {
    if (w.size() != expected)
        throw error(line, "expected '" + form + "'");
}

bool OptimizationTrial::Reader::yesNo(int line, const std::vector<std::string> &w, std::size_t index) {
    if (w.size() != index + 1)
        throw error(line, "expected yes or no after '" + join(w, 0, index) + "'");
    return yesNoValue(line, w[index]);
}

bool OptimizationTrial::Reader::yesNoValue(int line, const std::string &value) const {
    std::string text = lower(value);
    if (text == "yes" || text == "true" || text == "on")
        return true;
    if (text == "no" || text == "false" || text == "off")
        return false;
    throw error(line, "expected yes or no, found '" + value + "'");
}

double OptimizationTrial::Reader::number_(int line, const std::string &value, const std::string &what) const {
    double d = 0.0;
    if (!parseDouble(value, d))
        throw error(line, "expected a number for the " + what + ", found '" + value + "'");
    if (!std::isfinite(d))
        throw error(line, "the " + what + " must be finite");
    return d;
}

double OptimizationTrial::Reader::nonNegative(int line, const std::string &value, const std::string &what) {
    double d = number_(line, value, what);
    if (d < 0.0)
        throw error(line, "the " + what + " must not be negative");
    return d;
}

double OptimizationTrial::Reader::positive(int line, const std::string &value, const std::string &what) {
    double d = number_(line, value, what);
    if (d <= 0.0)
        throw error(line, "the " + what + " must be positive");
    return d;
}

std::vector<double> OptimizationTrial::Reader::numbers(int line, const std::vector<std::string> &w,
                                    std::size_t from) {
    if (from >= w.size())
        throw error(line, "expected values after '" + join(w, 0, from) + "'");
    std::vector<double> values;
    for (std::size_t i = from; i < w.size(); i++)
        values.push_back(number_(line, w[i], "value"));
    return values;
}

int OptimizationTrial::Reader::integer(int line, const std::string &value, const std::string &what) {
    int i = 0;
    if (!parseInt(value, i))
        throw error(line, "expected a whole number for the " + what + ", found '" + value +
                              "'");
    return i;
}

int OptimizationTrial::Reader::nonNegativeInt(int line, const std::string &value, const std::string &what) {
    int i = integer(line, value, what);
    if (i < 0)
        throw error(line, "the " + what + " must not be negative");
    return i;
}

int OptimizationTrial::Reader::positiveInt(int line, const std::string &value, const std::string &what) {
    int i = integer(line, value, what);
    if (i <= 0)
        throw error(line, "the " + what + " must be positive");
    return i;
}

std::vector<int> OptimizationTrial::Reader::positiveInts(int line, const std::vector<std::string> &w,
                                      std::size_t from, const std::string &what) {
    if (from >= w.size())
        throw error(line, "expected " + what + " after '" + join(w, 0, from) + "'");
    // Java's what.replaceAll("s$", ""): the singular for the per-value message.
    std::string singular = what;
    if (!singular.empty() && singular.back() == 's')
        singular.pop_back();
    std::vector<int> values;
    std::set<int> unique;
    for (std::size_t i = from; i < w.size(); i++) {
        int value = positiveInt(line, w[i], singular);
        if (!unique.insert(value).second)
            throw error(line, w[i] + " is listed twice");
        values.push_back(value);
    }
    return values;
}

std::vector<double> OptimizationTrial::Reader::fieldList(int line, const std::vector<std::string> &w) {
    std::vector<double> values;
    if (w.size() == 6 && lower(w[2]) == "to" && lower(w[4]) == "step") {
        Decimal from = decimal(line, w[1]);
        Decimal to = decimal(line, w[3]);
        Decimal step = decimal(line, w[5]);
        if (step.signum() <= 0)
            throw error(line, "the field step must be positive");
        for (Decimal value = from; value.compareTo(to) <= 0; value = value.add(step))
            values.push_back(value.toDouble());
    } else {
        if (w.size() < 2)
            throw error(line, "expected field values after 'fields'");
        for (std::size_t i = 1; i < w.size(); i++)
            values.push_back(number_(line, w[i], "field"));
    }
    for (double value : values)
        if (value < 0.0 || value > 1.0)
            throw error(line, "fields are relative heights, between 0 and 1");
    if (values.empty() || values[0] != 0.0)
        throw error(line, "the first field must be 0");
    return values;
}

Decimal OptimizationTrial::Reader::decimal(int line, const std::string &value) {
    Decimal result;
    if (!Decimal::parse(value, result))
        throw error(line, "expected a number, found '" + value + "'");
    return result;
}

// ---------------------------------------------------------------------------
// OptimizationTrial
// ---------------------------------------------------------------------------

OptimizationTrial::Trial OptimizationTrial::read(const std::string &text, int number,
                                                 bool useGlassTypes) {
    return parse(text, number).createBuilder(text, useGlassTypes);
}

OptimizationTrial::TrialDefinition OptimizationTrial::parse(const std::string &text,
                                                            int number) {
    auto settings = std::make_shared<const Reader>(Reader::parse(text, number));
    auto configuration =
        std::make_shared<const OptimizationConfiguration>(settings->toConfiguration());
    return TrialDefinition(std::move(settings), std::move(configuration));
}

OptimizationTrial::TrialDefinition::TrialDefinition(
    std::shared_ptr<const Reader> settings_,
    std::shared_ptr<const OptimizationConfiguration> configuration_)
    : settings(std::move(settings_)), configuration(std::move(configuration_)) {}

int OptimizationTrial::TrialDefinition::number() const {
    return settings->trialNumber();
}

OptimizationTrial::Trial OptimizationTrial::TrialDefinition::createBuilder(
    const std::string &prescriptionText, bool useGlassTypes) const {
    importers::OpticalBenchDataImporter::LensSpecifications specs;
    specs.parse_buffer(prescriptionText);
    auto prescription = std::make_unique<spec::Prescription>(spec::Prescription::build_prescription(
        specs, useGlassTypes, configuration->weighted, configuration->dLineOnly));
    OptimizationBuilder builder = settings->apply(prescription.get(), *configuration);
    return Trial{std::move(prescription), std::move(builder)};
}

std::string OptimizationTrial::TrialDefinition::toTrial() const {
    return configuration->toTrial(number());
}

std::string OptimizationTrial::TrialDefinition::toTrial(spec::Prescription *prescription) const {
    return settings->apply(prescription, *configuration).toTrial(number());
}

std::optional<OptimizationPipeline> OptimizationTrial::readPipeline(const std::string &text,
                                                                    int number) {
    SectionIndex index(text);
    const auto &lines = index.lines;
    const auto &trials = index.trials;
    const auto &pipelines = index.pipelines;
    auto startEntry = pipelines.find(number);
    if (startEntry == pipelines.end()) {
        if (trials.find(number) != trials.end())
            return std::nullopt;
        throw TrialException("there is no [trial " + intToString(number) + "] or [pipeline " +
                             intToString(number) + "] in this prescription; it defines " +
                             defined("trial", trials) + " and " +
                             defined("pipeline", pipelines));
    }
    int start = startEntry->second - 1;

    std::optional<std::string> description;
    std::optional<std::string> outdir;
    std::optional<std::vector<int>> stages;
    int stagesLine = 0;
    std::set<std::string> seen;
    for (std::size_t i = static_cast<std::size_t>(start) + 1;
         i < lines.size() && trim(lines[i]).rfind("[", 0) != 0; i++) {
        std::size_t hash = lines[i].find('#');
        std::string text_ =
            trim(hash != std::string::npos ? lines[i].substr(0, hash) : lines[i]);
        if (text_.empty())
            continue;
        int line = static_cast<int>(i) + 1;
        std::vector<std::string> w = words(text_);
        std::string keyword = lower(w[0]);
        if (!seen.insert(keyword).second)
            throw pipelineError(number, line, "'" + keyword + "' is given more than once");
        if (keyword == "description")
            description = trim(text_.substr(w[0].size()));
        else if (keyword == "outdir") {
            if (w.size() < 2)
                throw pipelineError(number, line, "expected 'outdir <directory>'");
            outdir = trim(text_.substr(w[0].size()));
        } else if (keyword == "trials") {
            if (w.size() < 2)
                throw pipelineError(number, line,
                                    "expected 'trials <number> <number> ...'");
            stages = std::vector<int>();
            stagesLine = line;
            for (std::size_t t = 1; t < w.size(); t++) {
                int stage = 0;
                if (!parseInt(w[t], stage))
                    throw pipelineError(number, line,
                                        "expected a trial number, found '" + w[t] + "'");
                if (pipelines.find(stage) != pipelines.end())
                    throw pipelineError(number, line,
                                        "stage " + intToString(stage) +
                                            " is a pipeline; a pipeline runs trials, not "
                                            "other pipelines");
                if (trials.find(stage) == trials.end())
                    throw pipelineError(number, line,
                                        "there is no [trial " + intToString(stage) +
                                            "] in this prescription");
                stages->push_back(stage);
            }
        } else
            throw pipelineError(number, line,
                                "unknown keyword '" + w[0] +
                                    "'; a pipeline takes description, outdir and trials");
    }
    if (!stages.has_value())
        throw TrialException("pipeline " + intToString(number) +
                             ": 'trials' is required, naming the trials to run in order");
    if (stages->empty())
        throw pipelineError(number, stagesLine, "a pipeline needs at least one trial");
    return OptimizationPipeline(number, description, outdir, *stages);
}

std::string OptimizationTrial::describe(const Var &variable) {
    if (const auto *radius = dynamic_cast<const VarRadius *>(&variable))
        return "surface " + intToString(radius->_surface_id) + " radius";
    if (const auto *thickness = dynamic_cast<const VarThickness *>(&variable))
        return "surface " + intToString(thickness->_surface_id) + " thickness";
    if (const auto *conic = dynamic_cast<const VarAsphK *>(&variable))
        return "surface " + intToString(conic->_surface_id) + " K";
    if (const auto *coefficient = dynamic_cast<const VarAsphCoeff *>(&variable))
        return "surface " + intToString(coefficient->_surface_id) + " coefficient " +
               intToString(coefficient->_index);
    return variable.toString();
}

std::string OptimizationTrial::paraxialName(int paraxId) {
    for (const auto &entry : PARAXIAL_QUANTITIES)
        if (entry.second == paraxId)
            return entry.first;
    throw IllegalArgumentException(std::string("paraxial quantity ") +
                                   ParaxHelper::Names[paraxId] + " has no name in a trial");
}

std::string OptimizationTrial::format(double value) {
    if (value == std::rint(value) && std::fabs(value) < 1e7)
        return std::to_string(static_cast<long long>(value));
    return doubleToString(value);
}

std::string OptimizationTrial::format(const std::vector<double> &values) {
    std::string result;
    for (double value : values) {
        if (!result.empty())
            result += " ";
        result += format(value);
    }
    return result;
}

std::string OptimizationTrial::format(const std::vector<int> &values) {
    std::string result;
    for (int value : values) {
        if (!result.empty())
            result += " ";
        result += intToString(value);
    }
    return result;
}

void OptimizationTrial::line(std::string &sb, const std::string &key,
                             const std::string &values) {
    sb += key;
    int padding = std::max(1, 22 - static_cast<int>(key.size()));
    sb.append(static_cast<std::size_t>(padding), ' ');
    sb += values;
    sb += "\n";
}

std::string OptimizationTrial::kebab(const std::string &name) {
    std::string sb;
    for (std::size_t i = 0; i < name.size(); i++) {
        char c = name[i];
        if (i > 0 && std::isupper(static_cast<unsigned char>(c)))
            sb.push_back('-');
        sb.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    }
    return sb;
}

std::vector<std::string> OptimizationTrial::splitTabs(std::string line) {
    std::vector<std::string> result;
    while (!line.empty()) {
        std::size_t pos = line.find('\t');
        if (pos == std::string::npos) {
            result.push_back(line);
            break;
        }
        result.push_back(line.substr(0, pos));
        line = line.substr(pos + 1);
    }
    return result;
}

} // namespace redukti::optim
