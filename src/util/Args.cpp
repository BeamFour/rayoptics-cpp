// C++ port of org.redukti.util.Args and org.redukti.util.Helper
#include "redukti/util/Args.h"

#include "redukti/Exceptions.h"

#include <cctype>
#include <climits>
#include <cerrno>
#include <cstdlib>
#include <filesystem>
#include <fstream>

namespace redukti::util {

using spec::VigType;

namespace {

/** The VigType constants in declaration order, with their Java enum names. */
struct VigTypeName {
    VigType value;
    const char *name;
};

const VigTypeName kVigTypes[] = {
    {VigType::None, "None"},
    {VigType::Paraxial, "Paraxial"},
    {VigType::SetVig, "SetVig"},
    {VigType::SetPupil, "SetPupil"},
    {VigType::SetStopAperture, "SetStopAperture"},
    {VigType::SetApertures, "SetApertures"},
    {VigType::SetFnum, "SetFnum"},
};

std::string strip(const std::string &s, char c) {
    std::string out;
    for (char ch : s)
        if (ch != c)
            out.push_back(ch);
    return out;
}

std::string lower(const std::string &s) {
    std::string out;
    for (char c : s)
        out.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    return out;
}

bool equalsIgnoreCase(const std::string &a, const std::string &b) {
    return lower(a) == lower(b);
}

std::string trim(const std::string &s) {
    std::size_t b = s.find_first_not_of(" \t\n\r\f\v");
    if (b == std::string::npos)
        return "";
    std::size_t e = s.find_last_not_of(" \t\n\r\f\v");
    return s.substr(b, e - b + 1);
}

std::string to_kebab_case(const std::string &name) {
    std::string sb;
    for (std::size_t i = 0; i < name.size(); i++) {
        char c = name[i];
        if (i > 0 && std::isupper(static_cast<unsigned char>(c)))
            sb.push_back('-');
        sb.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    }
    return sb;
}

} // namespace

Args Args::parseArguments(const std::vector<std::string> &args) {
    Args arguments;
    for (std::size_t i = 0; i < args.size(); i++) {
        const std::string &arg1 = args[i];
        std::optional<std::string> arg2;
        if (i + 1 < args.size())
            arg2 = args[i + 1];
        if (arg1 == "--specfile") {
            arguments.specfile = arg2;
            i++;
        } else if (arg1 == "-o") {
            arguments.outputFile = arg2;
            i++;
        } else if (arg1 == "--scenario") {
            arguments.scenario = std::atoi(arg2.value().c_str());
            i++;
        } else if (arg1 == "--outdir") {
            arguments.outdir = arg2;
            i++;
        } else if (arg1 == "--dont-use-glass-types") {
            arguments.use_glass_types = false;
        } else if (arg1 == "--force") {
            arguments.force = true;
        } else if (arg1 == "--only-d-line") {
            arguments.only_d_line = true;
        } else if (arg1 == "--output-ray-aberration-plots") {
            arguments.do_ray_aberrations = true;
        } else if (arg1 == "--output-wavelength-mtfs") {
            arguments.do_mono_chrome_mtfs = true;
        } else if (arg1 == "--mtf") {
            arguments.mtf_freqs = parse_mtf_freqs(arg2);
            i++;
        } else if (arg1 == "--use-spot-pattern") {
            arguments.spot_pattern = parse_spot_pattern(arg2);
            i++;
        } else if (arg1 == "--spot-grid-size") {
            arguments.spot_grid_size = parse_positive_int(arg1, arg2);
            i++;
        } else if (arg1 == "--auto-size-spot-diagrams") {
            arguments.auto_size_spots = true;
        } else if (arg1 == "--assign-glass-types") {
            arguments.assign_glass_types = true;
        } else if (arg1 == "--index-line") {
            arguments.index_line = parse_index_line(arg2);
            i++;
        } else if (arg1 == "--update-specfile") {
            arguments.update_specfile = true;
        } else if (arg1 == "--optimize") {
            arguments.optimize = true;
        } else if (arg1 == "--optimize-goal") {
            arguments.optimize_goal = parse_optimize_goal(arg2);
            i++;
        } else if (arg1 == "--vig-type") {
            arguments.vig_type = parse_vig_type(arg2);
            i++;
        } else if (arg1 == "--real-ray-aiming") {
            arguments.real_ray_aiming = true;
        } else if (arg1 == "--paraxial-ray-aiming") {
            arguments.real_ray_aiming = false;
        } else if (arg1 == "--generate-java") {
            arguments.generate_java = true;
        } else if (arg1 == "--legacy-notebook") {
            arguments.legacy_notebook = true;
        } else if (arg1 == "--reference") {
            arguments.reference_file = arg2;
            i++;
        }
    }
    return arguments;
}

rayoptics::seq::Glass::IndexLine Args::index_line_value() const {
    return index_line == "e" ? rayoptics::seq::Glass::IndexLine::E
                             : rayoptics::seq::Glass::IndexLine::D;
}

std::string Args::parse_index_line(const std::optional<std::string> &value) {
    if (!value.has_value())
        throw IllegalArgumentException("--index-line requires a value, one of: d, e");
    std::string normalized = lower(trim(*value));
    if (normalized == "d" || normalized == "e")
        return normalized;
    throw IllegalArgumentException("Unrecognized --index-line '" + *value +
                                   "', expected one of: d, e");
}

std::string Args::parse_optimize_goal(const std::optional<std::string> &value) {
    if (!value.has_value())
        throw IllegalArgumentException(
            "--optimize-goal requires a value, one of: contrast, mtf");
    std::string normalized = lower(trim(*value));
    if (normalized == "contrast" || normalized == "mtf")
        return normalized;
    throw IllegalArgumentException("Unrecognized --optimize-goal '" + *value +
                                   "', expected one of: contrast, mtf");
}

int Args::parse_positive_int(const std::string &option,
                             const std::optional<std::string> &value) {
    if (!value.has_value())
        throw IllegalArgumentException(option + " requires a positive integer");
    // Integer.parseInt: an optional sign, then digits, nothing else -- no
    // surrounding whitespace, and out-of-range values rejected rather than
    // clamped. Anything it would reject is reported the same way.
    const std::string &text = *value;
    bool wellFormed = !text.empty();
    std::size_t digitsFrom = (!text.empty() && (text[0] == '+' || text[0] == '-')) ? 1 : 0;
    if (digitsFrom == text.size())
        wellFormed = false;
    for (std::size_t k = digitsFrom; wellFormed && k < text.size(); k++)
        if (!std::isdigit(static_cast<unsigned char>(text[k])))
            wellFormed = false;
    long long parsed = 0;
    if (wellFormed) {
        errno = 0;
        parsed = std::strtoll(text.c_str(), nullptr, 10);
        if (errno == ERANGE || parsed < INT_MIN || parsed > INT_MAX)
            wellFormed = false;
    }
    if (!wellFormed)
        throw IllegalArgumentException(option + " requires a positive integer, got '" +
                                       text + "'");
    if (parsed < 2)
        throw IllegalArgumentException(option + " must be at least 2, got " +
                                       std::to_string(parsed));
    return static_cast<int>(parsed);
}

std::vector<int> Args::parse_mtf_freqs(const std::optional<std::string> &value) {
    if (!value.has_value())
        throw IllegalArgumentException(
            "--mtf requires a comma separated list of frequencies in cycles/mm, e.g. "
            "10,30,50");
    std::vector<std::string> parts;
    std::size_t start = 0;
    while (true) {
        std::size_t c = value->find(',', start);
        if (c == std::string::npos) {
            parts.push_back(value->substr(start));
            break;
        }
        parts.push_back(value->substr(start, c - start));
        start = c + 1;
    }
    std::vector<int> freqs;
    for (const auto &raw : parts) {
        std::string part = trim(raw);
        const char *s = part.c_str();
        char *end = nullptr;
        long v = std::strtol(s, &end, 10);
        if (end == s || *end != '\0')
            throw IllegalArgumentException("Unrecognized --mtf frequency '" + part +
                                           "' in '" + *value +
                                           "', expected a comma separated list such as "
                                           "10,30,50");
        int freq = static_cast<int>(v);
        if (freq <= 0)
            throw IllegalArgumentException("--mtf frequency must be positive, got " +
                                           std::to_string(freq));
        freqs.push_back(freq);
    }
    if (freqs.empty())
        throw IllegalArgumentException("--mtf requires at least one frequency, e.g. "
                                       "10,30,50");
    return freqs;
}

VigType Args::parse_vig_type(const std::optional<std::string> &value) {
    if (!value.has_value())
        throw IllegalArgumentException("--vig-type requires a value, one of: " +
                                       vig_type_names());
    std::string normalized = strip(strip(*value, '-'), '_');
    for (const auto &vt : kVigTypes) {
        if (equalsIgnoreCase(vt.name, normalized))
            return vt.value;
    }
    throw IllegalArgumentException("Unrecognized --vig-type '" + *value +
                                   "', expected one of: " + vig_type_names());
}

std::string Args::vig_type_names() {
    std::string sb;
    for (const auto &vt : kVigTypes) {
        if (!sb.empty())
            sb += ", ";
        sb += to_kebab_case(vt.name);
    }
    return sb;
}

int Args::parse_spot_pattern(const std::optional<std::string> &value) {
    if (!value.has_value())
        throw IllegalArgumentException("--use-spot-pattern requires a value, one of: " +
                                       spot_pattern_names());
    std::string normalized = lower(strip(strip(*value, '-'), '_'));
    if (normalized == "hex" || normalized == "hexapolar")
        return 1; // PATTERN_HEXAPOLAR
    if (normalized == "grid")
        return 3; // PATTERN_GRID
    if (normalized == "gq" || normalized == "gauss" || normalized == "gaussian" ||
        normalized == "gaussianquadrature")
        return 2; // PATTERN_GAUSS_QUADRATURE
    throw IllegalArgumentException("Unrecognized --use-spot-pattern '" + *value +
                                   "', expected one of: " + spot_pattern_names());
}

std::string Args::spot_pattern_names() { return "hex, grid, gaussian"; }

// ---------------------------------------------------------------------------
// Helper
// ---------------------------------------------------------------------------

namespace fs = std::filesystem;

std::string Helper::getOutputPath(const Args &arguments) {
    if (!arguments.outputFile.has_value())
        throw IllegalArgumentException("Output file name not specified");
    if (arguments.outdir.has_value())
        return (fs::path(*arguments.outdir) / *arguments.outputFile).string();
    fs::path path = fs::absolute(fs::path(*arguments.specfile));
    return (path.parent_path() / *arguments.outputFile).string();
}

std::string Helper::getOutputPath(const Args &arguments, const std::string &extension) {
    if (!arguments.outputFile.has_value())
        throw IllegalArgumentException("Output file name not specified");
    if (arguments.outdir.has_value())
        return (fs::path(*arguments.outdir) / (*arguments.outputFile + extension))
            .string();
    fs::path path = fs::absolute(fs::path(*arguments.specfile));
    return (path.parent_path() / (*arguments.outputFile + extension)).string();
}

std::string Helper::getOutputFileWithPath(const std::string &specfile,
                                          const std::string &outputFile,
                                          const std::optional<std::string> &outdir) {
    if (outdir.has_value())
        return (fs::path(*outdir) / outputFile).string();
    fs::path path = fs::absolute(fs::path(specfile));
    return (path.parent_path() / outputFile).string();
}

std::string Helper::getFilename(const std::string &specfile) {
    return fs::absolute(fs::path(specfile)).filename().string();
}

std::string Helper::replaceExtension(const std::string &fileName_in,
                                     const std::string &extension) {
    std::string fileName = fileName_in;
    std::size_t dotIndex = fileName.rfind('.');
    if (dotIndex != std::string::npos && dotIndex > 0 && dotIndex < fileName.size() - 1)
        fileName = fileName.substr(0, dotIndex);
    return fileName + extension;
}

std::string Helper::getOutputPathChangeExt(const std::string &specfile,
                                           const std::string &extension) {
    fs::path path = fs::absolute(fs::path(specfile));
    std::string fileName = path.filename().string();
    std::size_t dotIndex = fileName.rfind('.');
    if (dotIndex != std::string::npos && dotIndex > 0 && dotIndex < fileName.size() - 1)
        fileName = fileName.substr(0, dotIndex);
    return (path.parent_path() / (fileName + extension)).string();
}

void Helper::createOutputFile(const std::string &outpath, const std::string &content) {
    // Binary so the "\n" the renderers emit is not translated to CRLF on
    // Windows: the committed reference files use bare newlines.
    std::ofstream out(outpath, std::ios::binary | std::ios::trunc);
    if (!out)
        throw RuntimeException("Failed to create file " + outpath);
    out.write(content.data(), static_cast<std::streamsize>(content.size()));
    if (!out)
        throw RuntimeException("Failed to write file " + outpath);
}

} // namespace redukti::util
