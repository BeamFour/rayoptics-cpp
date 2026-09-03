// C++ port of org.redukti.importers.obench.OpticalBenchDataImporter
#include "redukti/importers/OpticalBenchDataImporter.h"

#include "redukti/Exceptions.h"
#include "redukti/mathlib/M.h"
#include "redukti/rayoptics/seq/Glass.h"

#include <cctype>
#include <cstdlib>
#include <fstream>
#include <map>

namespace redukti::importers {

using rayoptics::seq::Glass;

// ---------------------------------------------------------------------------
// Java number parsing
// ---------------------------------------------------------------------------

double parse_double(const std::string &s) {
    if (s.empty())
        return 0.0;
    // Java allows surrounding whitespace and nothing else around the number.
    std::size_t b = s.find_first_not_of(" \t\n\r\f\v");
    if (b == std::string::npos)
        return 0.0;
    std::size_t e = s.find_last_not_of(" \t\n\r\f\v");
    std::string t = s.substr(b, e - b + 1);
    if (t.empty())
        return 0.0;

    // Java spells these words out; strtod would accept "inf"/"infinity"
    // case-insensitively, which Java does not.
    std::string body = t;
    bool neg = false;
    if (body[0] == '+' || body[0] == '-') {
        neg = body[0] == '-';
        body = body.substr(1);
    }
    if (body == "Infinity")
        return neg ? -std::numeric_limits<double>::infinity()
                   : std::numeric_limits<double>::infinity();
    if (body == "NaN")
        return std::numeric_limits<double>::quiet_NaN();
    // Reject the spellings strtod would take but Java rejects.
    for (char c : body) {
        if (c == 'i' || c == 'I' || c == 'n' || c == 'N' || c == 'x' || c == 'X')
            return 0.0;
    }
    // Java permits one trailing type suffix.
    if (!body.empty()) {
        char last = body.back();
        if (last == 'd' || last == 'D' || last == 'f' || last == 'F')
            body.pop_back();
    }
    if (body.empty())
        return 0.0;
    std::string full = (neg ? "-" : "") + body;
    const char *start = full.c_str();
    char *end = nullptr;
    double v = std::strtod(start, &end);
    // Anything left over means Java would have thrown.
    if (end == start || *end != '\0')
        return 0.0;
    return v;
}

int parse_integer(const std::string &s, int defaultValue) {
    if (s.empty())
        return defaultValue;
    const char *start = s.c_str();
    char *end = nullptr;
    long v = std::strtol(start, &end, 10);
    if (end == start || *end != '\0')
        return defaultValue;
    return static_cast<int>(v);
}

// ---------------------------------------------------------------------------
// VarSet
// ---------------------------------------------------------------------------

OpticalBenchDataImporter::Variable *OpticalBenchDataImporter::VarSet::add_variable(
    const std::string &name) {
    variables_.push_back(Variable(name));
    return &variables_.back();
}

OpticalBenchDataImporter::Variable *OpticalBenchDataImporter::VarSet::find_variable(
    const std::string &name) {
    for (auto &v : variables_) {
        if (name == v.name())
            return &v;
    }
    return nullptr;
}

const OpticalBenchDataImporter::Variable *
OpticalBenchDataImporter::VarSet::find_variable(const std::string &name) const {
    for (const auto &v : variables_) {
        if (name == v.name())
            return &v;
    }
    return nullptr;
}

std::string OpticalBenchDataImporter::VarSet::get_value(const std::string &name) const {
    const Variable *variable = find_variable(name);
    if (variable != nullptr)
        return variable->get_value(0);
    return "";
}

// ---------------------------------------------------------------------------
// AsphericalData
// ---------------------------------------------------------------------------

std::vector<double> OpticalBenchDataImporter::AsphericalData::get_coeffs() const {
    int a = 0;
    if (get_asphere_type() == AsphereType::Odd)
        a = 2;
    else if (get_asphere_type() == AsphereType::Even)
        a = 1;
    std::vector<double> coeffs(static_cast<std::size_t>(
                                   static_cast<int>(_data.size()) - 2 + a),
                               0.0);
    for (int i = 2; i < static_cast<int>(_data.size()); i++, a++)
        coeffs[static_cast<std::size_t>(a)] = data(i);
    return coeffs;
}

// ---------------------------------------------------------------------------
// LensSurface
// ---------------------------------------------------------------------------

double OpticalBenchDataImporter::LensSurface::get_thickness(int scenario) const {
    if (scenario < static_cast<int>(_thickness_by_scenario.size()))
        return _thickness_by_scenario[static_cast<std::size_t>(scenario)];
    // Java asserts there is exactly one and falls back to it.
    return _thickness_by_scenario[0];
}

double OpticalBenchDataImporter::LensSurface::get_diameter(int scenario) const {
    if (scenario < static_cast<int>(_diameter_by_scenario.size()))
        return _diameter_by_scenario[static_cast<std::size_t>(scenario)];
    return _diameter_by_scenario[0];
}

// ---------------------------------------------------------------------------
// Sections
// ---------------------------------------------------------------------------

namespace {

enum class Section {
    DESCRIPTIVE_DATA,
    CONSTANTS,
    VARIABLE_DISTANCES,
    LENS_DATA,
    ASPHERICAL_DATA,
    PATENT_INFO,
    REPORT_DATA,
};

struct SectionMapping {
    const char *name;
    Section section;
};

const SectionMapping g_SectionMappings[] = {
    {"[descriptive data]", Section::DESCRIPTIVE_DATA},
    {"[constants]", Section::CONSTANTS},
    {"[variable distances]", Section::VARIABLE_DISTANCES},
    {"[lens data]", Section::LENS_DATA},
    {"[aspherical data]", Section::ASPHERICAL_DATA},
    {"[patent info]", Section::PATENT_INFO},
    {"[report data]", Section::REPORT_DATA},
};

/** Null when the header is not one this importer knows. */
const Section *find_section(const std::string &name) {
    for (const auto &m : g_SectionMappings) {
        if (name == m.name)
            return &m.section;
    }
    return nullptr;
}

} // namespace

// ---------------------------------------------------------------------------
// LensSpecifications
// ---------------------------------------------------------------------------

std::vector<std::string> OpticalBenchDataImporter::LensSpecifications::splitLine(
    const std::string &line_in) {
    std::vector<std::string> words;
    std::string line = line_in;
    while (line.length() > 0) {
        std::size_t pos = line.find('\t');
        if (pos == std::string::npos) {
            words.push_back(line);
            break;
        } else if (pos == 0) {
            words.push_back("");
            line = line.substr(1);
        } else {
            words.push_back(line.substr(0, pos));
            line = line.substr(pos + 1);
        }
    }
    return words;
}

bool OpticalBenchDataImporter::LensSpecifications::parse_file(
    const std::string &file_name) {
    std::ifstream in(file_name, std::ios::binary);
    if (!in)
        throw RuntimeException("Cannot read " + file_name);
    std::string content((std::istreambuf_iterator<char>(in)),
                        std::istreambuf_iterator<char>());
    return parse_buffer(content);
}

bool OpticalBenchDataImporter::LensSpecifications::parse_buffer(
    const std::string &buffer) {
    std::vector<std::string> lines;
    std::size_t start = 0;
    while (start <= buffer.size()) {
        std::size_t nl = buffer.find('\n', start);
        std::string line;
        if (nl == std::string::npos) {
            line = buffer.substr(start);
            start = buffer.size() + 1;
        } else {
            line = buffer.substr(start, nl - start);
            start = nl + 1;
        }
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        lines.push_back(line);
    }
    // Java splits on \r?\n, which drops a trailing empty line only when the
    // buffer does not end in a separator; an empty trailing line is harmless
    // here because splitLine gives it no words.
    return parse_lines(lines);
}

bool OpticalBenchDataImporter::LensSpecifications::parse_lines(
    const std::vector<std::string> &lines) {
    const Section *current_section = nullptr;
    int surface_id = 1;
    std::map<std::string, int> surfaceIdMap;
    AsphereType asphere_type = AsphereType::Even;
    for (const std::string &line : lines) {
        std::vector<std::string> words = splitLine(line);
        if (words.empty())
            continue;
        if (!words[0].empty() && words[0][0] == '#')
            continue;
        if (!words[0].empty() && words[0][0] == '[') {
            current_section = find_section(words[0]);
            continue;
        }
        if (current_section == nullptr)
            continue;
        switch (*current_section) {
        case Section::DESCRIPTIVE_DATA:
            if (words.size() >= 2) {
                Variable *var = descriptive_data_.add_variable(words[0]);
                for (std::size_t i = 1; i < words.size(); i++)
                    var->add_value(words[i]);
            }
            break;
        case Section::CONSTANTS: {
            Variable *var = constants_.add_variable(words[0]);
            for (std::size_t i = 1; i < words.size(); i++)
                var->add_value(words[i]);
            break;
        }
        case Section::VARIABLE_DISTANCES:
            if (words.size() >= 2) {
                Variable *var = variables_.add_variable(words[0]);
                for (std::size_t i = 1; i < words.size(); i++)
                    var->add_value(words[i]);
            }
            break;
        case Section::PATENT_INFO:
            if (words.size() >= 2) {
                Variable *var = patent_info_.add_variable(words[0]);
                for (std::size_t i = 1; i < words.size(); i++)
                    var->add_value(words[i]);
            }
            break;
        case Section::REPORT_DATA:
            if (words.size() >= 2) {
                Variable *var = report_data_.add_variable(words[0]);
                for (std::size_t i = 1; i < words.size(); i++)
                    var->add_value(words[i]);
            }
            break;
        case Section::LENS_DATA: {
            if (words.size() < 2)
                break;
            int id = surface_id++;
            surfaceIdMap[words[0]] = id;
            LensSurface surface_data(id);
            SurfaceType type = SurfaceType::surface;
            if (words[1] == "AS") {
                type = SurfaceType::aperture_stop;
                surface_data.set_radius(0.0);
            } else if (words[1] == "FS") {
                type = SurfaceType::field_stop;
                surface_data.set_radius(0.0);
            } else if (words[1] == "CG") {
                surface_data.set_radius(0.0);
                surface_data.set_is_cover_glass(true);
            } else {
                if (words[1] == "Infinity")
                    surface_data.set_radius(0.0);
                else
                    surface_data.set_radius(parse_double(words[1]));
            }
            surface_data.set_surface_type(type);
            if (words.size() >= 3 && !words[2].empty())
                parse_thickness(words[2], surface_data);
            if (words.size() >= 4 && !words[3].empty())
                surface_data.set_refractive_index(parse_double(words[3]));
            if (words.size() >= 5 && !words[4].empty())
                parse_diameter(words[4], type == SurfaceType::aperture_stop,
                               surface_data);
            if (words.size() >= 6 && !words[5].empty())
                surface_data.set_abbe_vd(parse_double(words[5]));
            if (words.size() >= 7 && !words[6].empty())
                surface_data.set_glass_name(words[6]);
            if (words.size() >= 8 && !words[7].empty())
                surface_data.set_catalog_name(Glass::get_catalog_name(words[7]));
            surfaces_.push_back(surface_data);
            break;
        }
        case Section::ASPHERICAL_DATA: {
            if (has_constant("AsphericalOddCount"))
                asphere_type = AsphereType::Odd;
            else if (has_constant("AsphericalA2"))
                asphere_type = AsphereType::EvenA2;
            else
                asphere_type = AsphereType::Even;
            const std::string &optBenchID = words[0];
            auto it = surfaceIdMap.find(optBenchID);
            if (it == surfaceIdMap.end())
                // Java indexes the map and NPEs on a miss; this says why.
                throw RuntimeException("Unknown surface " + optBenchID);
            int id = it->second;
            auto aspherical_data = std::make_unique<AsphericalData>(asphere_type, id);
            for (std::size_t i = 1; i < words.size(); i++)
                aspherical_data->add_data(parse_double(words[i]));
            aspherical_data_.push_back(std::move(aspherical_data));
            LensSurface *surface_builder = find_surface(id);
            if (surface_builder == nullptr)
                throw RuntimeException("Unknown surface " + optBenchID);
            surface_builder->set_aspherical_data(aspherical_data_.back().get());
            break;
        }
        }
    }
    return true;
}

OpticalBenchDataImporter::LensSurface *
OpticalBenchDataImporter::LensSpecifications::find_surface(int id) {
    for (auto &s : surfaces_) {
        if (s.get_id() == id)
            return &s;
    }
    return nullptr;
}

void OpticalBenchDataImporter::LensSpecifications::parse_thickness(
    const std::string &value, LensSurface &surface_builder) {
    if (value.length() == 0) {
        surface_builder.add_thickness(0.0);
        return;
    }
    if (std::isalpha(static_cast<unsigned char>(value[0]))) {
        Variable *var = find_variable(value);
        if (var != nullptr) {
            for (int i = 0; i < var->num_scenarios(); i++)
                surface_builder.add_thickness(parse_double(var->get_value(i)));
        } else {
            surface_builder.add_thickness(0.0);
        }
    } else {
        surface_builder.add_thickness(parse_double(value));
    }
}

void OpticalBenchDataImporter::LensSpecifications::parse_diameter(
    const std::string &value, bool isApertureStop, LensSurface &surface_builder) {
    double dValue = parse_double(value);
    if (!isApertureStop) {
        surface_builder.set_diameter(dValue);
    } else {
        Variable *var = find_variable("Aperture Diameter");
        if (var != nullptr) {
            for (int i = 0; i < var->num_scenarios(); i++)
                surface_builder.set_diameter(parse_double(var->get_value(i)));
        } else {
            surface_builder.set_diameter(dValue);
        }
    }
}

double OpticalBenchDataImporter::LensSpecifications::get_image_height() const {
    const Variable *var = find_variable("Image Height");
    if (var != nullptr)
        return var->get_value_as_double(0);
    return 43.2;
}

double OpticalBenchDataImporter::LensSpecifications::get_focal_length() const {
    const Variable *var = find_variable("Focal Length");
    if (var != nullptr)
        return var->get_value_as_double(0);
    throw IllegalArgumentException("Focal Length");
}

double OpticalBenchDataImporter::LensSpecifications::get_focal_length(
    int scenario) const {
    return find_variable("Focal Length")->get_value_as_double(scenario);
}

double OpticalBenchDataImporter::LensSpecifications::get_stop_diameter(
    int scenario) const {
    return find_variable("Aperture Diameter")->get_value_as_double(scenario);
}

double OpticalBenchDataImporter::LensSpecifications::get_angle_of_view_in_degrees(
    int scenario) const {
    return find_variable("Angle of View")->get_value_as_double(scenario);
}

double OpticalBenchDataImporter::LensSpecifications::get_f_number(int scenario) const {
    return find_variable("F-Number")->get_value_as_double(scenario);
}

double OpticalBenchDataImporter::LensSpecifications::get_half_angle_of_view_in_radians(
    int scenario) const {
    return mathlib::M::toRadians(
        find_variable("Angle of View")->get_value_as_double(scenario) / 2.0);
}

double OpticalBenchDataImporter::LensSpecifications::get_half_angle_of_view_in_degrees(
    int scenario) const {
    return find_variable("Angle of View")->get_value_as_double(scenario) / 2.0;
}

} // namespace redukti::importers
