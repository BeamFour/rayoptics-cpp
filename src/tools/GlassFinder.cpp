// C++ port of org.redukti.tools.GlassFinder
#include "redukti/tools/GlassFinder.h"

#include "redukti/Exceptions.h"
#include "redukti/Text.h"

#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <vector>

namespace redukti::tools {

namespace {

using rayoptics::seq::Glass;

const char *const LENS_DATA_SECTION = "[lens data]";
const std::string CANDIDATE_PREFIX = "candidate=";

/** Java's String.trim(): strips characters <= ' ' from both ends. */
std::string javaTrim(const std::string &s) {
    std::size_t b = 0, e = s.size();
    while (b < e && static_cast<unsigned char>(s[b]) <= ' ')
        b++;
    while (e > b && static_cast<unsigned char>(s[e - 1]) <= ' ')
        e--;
    return s.substr(b, e - b);
}

/** Java's String.isBlank(). */
bool isBlankString(const std::string &s) {
    for (char c : s)
        if (!std::isspace(static_cast<unsigned char>(c)))
            return false;
    return true;
}

bool equalsIgnoreCase(const std::string &a, const std::string &b) {
    if (a.size() != b.size())
        return false;
    for (std::size_t i = 0; i < a.size(); i++)
        if (std::tolower(static_cast<unsigned char>(a[i])) !=
            std::tolower(static_cast<unsigned char>(b[i])))
            return false;
    return true;
}

bool startsWith(const std::string &s, const std::string &prefix) {
    return s.compare(0, prefix.size(), prefix) == 0;
}

/**
 * Java's input.split("\\r?\\n", -1): every line, keeping trailing empty
 * strings, with a carriage return dropped only where it sits right before a
 * line feed.
 */
std::vector<std::string> splitLines(const std::string &input) {
    std::vector<std::string> lines;
    std::size_t start = 0;
    while (true) {
        std::size_t nl = input.find('\n', start);
        if (nl == std::string::npos) {
            lines.push_back(input.substr(start));
            break;
        }
        std::size_t end = nl;
        if (end > start && input[end - 1] == '\r')
            end--;
        lines.push_back(input.substr(start, end - start));
        start = nl + 1;
    }
    return lines;
}

/** Java's line.split("\\t", -1): keeps empty fields, trailing ones included. */
std::vector<std::string> splitTabs(const std::string &line) {
    std::vector<std::string> fields;
    std::size_t start = 0;
    while (true) {
        std::size_t tab = line.find('\t', start);
        if (tab == std::string::npos) {
            fields.push_back(line.substr(start));
            break;
        }
        fields.push_back(line.substr(start, tab - start));
        start = tab + 1;
    }
    return fields;
}

std::string join(const std::vector<std::string> &parts, const std::string &separator) {
    std::string out;
    for (std::size_t i = 0; i < parts.size(); i++) {
        if (i > 0)
            out += separator;
        out += parts[i];
    }
    return out;
}

/**
 * Java's GlassFinder.parseDouble: null for a blank field or anything
 * Double.parseDouble rejects. Follows the same rules as the importer's number
 * parser -- surrounding whitespace, the spelled-out Infinity and NaN, one type
 * suffix -- but reports a rejection as empty rather than as 0.0, which is the
 * distinction this caller needs.
 */
std::optional<double> parseDouble(const std::string &value) {
    if (isBlankString(value))
        return std::nullopt;
    std::string body = javaTrim(value);
    bool neg = false;
    if (!body.empty() && (body[0] == '+' || body[0] == '-')) {
        neg = body[0] == '-';
        body.erase(0, 1);
    }
    if (body == "Infinity")
        return neg ? -std::numeric_limits<double>::infinity()
                   : std::numeric_limits<double>::infinity();
    if (body == "NaN")
        return std::numeric_limits<double>::quiet_NaN();
    for (char c : body)
        if (c == 'i' || c == 'I' || c == 'n' || c == 'N' || c == 'x' || c == 'X')
            return std::nullopt;
    if (!body.empty()) {
        char last = body.back();
        if (last == 'd' || last == 'D' || last == 'f' || last == 'F')
            body.pop_back();
    }
    if (body.empty())
        return std::nullopt;
    std::string full = (neg ? "-" : "") + body;
    char *end = nullptr;
    double v = std::strtod(full.c_str(), &end);
    if (end == full.c_str() || *end != '\0')
        return std::nullopt;
    return v;
}

std::string formatCandidate(const Glass::GlassMatch &match) {
    const auto &glass = match.glass;
    // Java string concatenation renders a null as "null".
    return "candidate=" + glass->catalog_name.value_or("null") + "/" +
           glass->label.value_or("null") +
           ",nd=" + formatF(glass->nd, 5) + ",vd=" + formatF(glass->vd, 2) +
           ",dnd=" + formatF(match.nd_difference, 5) +
           ",dvd=" + formatF(match.vd_difference, 2);
}

std::vector<std::string> appendCandidates(std::vector<std::string> fields,
                                          const std::vector<Glass::GlassMatch> &matches) {
    for (const auto &match : matches)
        fields.push_back(formatCandidate(match));
    return fields;
}

std::vector<std::string> removeCandidateFields(const std::vector<std::string> &fields) {
    std::vector<std::string> retained;
    for (const auto &field : fields)
        if (!startsWith(field, CANDIDATE_PREFIX))
            retained.push_back(field);
    return retained;
}

bool isEmpty(const std::vector<std::string> &fields, std::size_t index) {
    return index >= fields.size() || isBlankString(fields[index]);
}

std::vector<std::string> ensureLength(std::vector<std::string> fields, std::size_t length) {
    if (fields.size() < length)
        fields.resize(length, "");
    return fields;
}

std::string readFile(const std::string &path) {
    std::ifstream in(path, std::ios::binary);
    if (!in)
        throw IOException("Failed to read " + path);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

} // namespace

GlassFinder::EnrichmentResult GlassFinder::enrich(const std::string &input, bool force) {
    return enrich(input, force, Glass::IndexLine::D);
}

GlassFinder::EnrichmentResult GlassFinder::enrich(const std::string &input, bool force,
                                                  Glass::IndexLine indexLine) {
    std::string newline = input.find("\r\n") != std::string::npos ? "\r\n" : "\n";
    bool endsWithNewline = !input.empty() && input.back() == '\n';
    std::vector<std::string> lines = splitLines(input);
    std::vector<std::string> output;
    output.reserve(lines.size());
    bool inLensData = false;
    int selected = 0;
    int ambiguous = 0;
    int unmatched = 0;

    std::size_t lineCount = endsWithNewline ? lines.size() - 1 : lines.size();
    for (std::size_t i = 0; i < lineCount; i++) {
        const std::string &line = lines[i];
        std::string trimmed = javaTrim(line);
        if (startsWith(trimmed, "[") && !trimmed.empty() && trimmed.back() == ']')
            inLensData = equalsIgnoreCase(trimmed, LENS_DATA_SECTION);

        if (!inLensData || trimmed.empty()) {
            output.push_back(line);
            continue;
        }

        std::vector<std::string> fields = splitTabs(line);
        if (fields.size() < 6 || fields[1] == "AS" || fields[1] == "FS" ||
            fields[1] == "CG") {
            output.push_back(line);
            continue;
        }

        auto nd = parseDouble(fields[3]);
        auto vd = parseDouble(fields[5]);
        if (!nd.has_value() || !vd.has_value() || *nd == 0.0 || *vd == 0.0) {
            output.push_back(line);
            continue;
        }

        if (!isEmpty(fields, 6) && !isEmpty(fields, 7) && !force &&
            Glass::get_catalog_name(fields[7]).has_value()) {
            output.push_back(line);
            continue;
        }

        auto matches = Glass::find_glasses(*nd, *vd, indexLine);
        if (matches.empty()) {
            unmatched++;
            output.push_back(line);
            continue;
        }

        std::vector<Glass::GlassMatch> exactMatches;
        for (const auto &match : matches)
            if (match.exact)
                exactMatches.push_back(match);
        fields = removeCandidateFields(fields);
        if (!exactMatches.empty() || matches.size() == 1) {
            const auto &glass = (!exactMatches.empty() ? exactMatches : matches)[0].glass;
            fields = ensureLength(fields, 8);
            // String.join renders a null element as "null".
            fields[6] = glass->label.value_or("null");
            fields[7] = glass->catalog_name.value_or("null");
            if (indexLine == Glass::IndexLine::E) {
                // The columns held ne and vd; restate them at the d line so that
                // anything reading this file without knowing the original
                // convention still gets the right medium.
                fields[3] = formatF(glass->nd, 5);
                fields[5] = formatF(glass->vd, 2);
            }
            selected++;
            if (matches.size() > 1)
                fields = appendCandidates(fields, matches);
        } else {
            fields = ensureLength(fields, 8);
            fields = appendCandidates(fields, matches);
            ambiguous++;
        }
        output.push_back(join(fields, "\t"));
    }

    std::string text = join(output, newline) + (endsWithNewline ? newline : "");
    return EnrichmentResult{text, selected, ambiguous, unmatched};
}

void GlassFinder::run(const util::Args &arguments) {
    std::string input = *arguments.specfile;
    std::string output = arguments.outputFile.has_value()
                             ? *arguments.outputFile
                             : util::Helper::getOutputFileWithPath(*arguments.specfile,
                                                                   "specs.txt", std::nullopt);

    EnrichmentResult result =
        enrich(readFile(input), arguments.force, arguments.index_line_value());
    util::Helper::createOutputFile(output, result.text);
    std::cout << "Selected " << result.selected << " glass types; " << result.ambiguous
              << " ambiguous; " << result.unmatched << " unmatched" << std::endl;
}

void GlassFinder::usage() {
    std::cerr << "Usage: GlassFinder --specfile input.txt -o output.txt [--index-line d|e] "
                 "[--force]"
              << std::endl;
    std::cerr << "       --index-line e when the prescription quotes the refractive index "
                 "at the e line"
              << std::endl;
}

} // namespace redukti::tools
