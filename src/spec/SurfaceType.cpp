// C++ port of org.redukti.spec.SurfaceType
#include "redukti/spec/SurfaceType.h"

#include "redukti/Text.h"
#include "redukti/rayoptics/seq/Glass.h"

namespace redukti::spec {

using rayoptics::seq::Glass;

namespace {
/** Java's StringBuilder.append(double), which is Double.toString. */
std::string d(double v) { return doubleToString(v); }
/** Java's StringBuilder.append(int). */
std::string i(int v) { return std::to_string(v); }
} // namespace

double SurfaceType::get_diameter_by_scenario(int scenario) const {
    if (_diameter_by_scenario.has_value())
        return (*_diameter_by_scenario)[static_cast<std::size_t>(scenario)];
    return _diameter;
}

double SurfaceType::get_thickness_by_scenario(int scenario) const {
    if (_thickness_by_scenario.has_value())
        return (*_thickness_by_scenario)[static_cast<std::size_t>(scenario)];
    return _thickness;
}

std::string &SurfaceType::to_opt_bench_str(std::string &sb, bool is_last) const {
    sb += _id;
    sb += "\t";
    if (_is_aperture_stop)
        sb += "AS";
    else if (_is_field_stop)
        sb += "FS";
    else
        sb += d(_radius);
    sb += "\t";
    if (_thickness_by_scenario.has_value()) {
        if (is_last)
            sb += "Bf";
        else
            sb += "d" + _id;
    } else {
        sb += d(_thickness);
    }
    sb += "\t";
    std::shared_ptr<Glass> glass;
    double nd = _nd;
    double vd = _vd;
    if (nd != 0.0 && _glass_name.has_value()) {
        glass = Glass::glass_by_catalog_name(_catalog_name, *_glass_name);
        if (glass != nullptr) {
            nd = glass->nd;
            vd = glass->vd;
        }
    }
    if (nd != 0.0)
        sb += d(nd);
    sb += "\t";
    sb += d(_diameter);
    sb += "\t";
    if (nd != 0.0)
        sb += d(vd);
    sb += "\t";
    if (glass != nullptr) {
        sb += *_glass_name;
        sb += "\t";
        // Always present on a catalog glass, which is the only way to get here.
        sb += glass->catalog_name.value();
    }
    sb += "\n";
    return sb;
}

std::string &SurfaceType::aspherics_to_opt_bench_str(std::string &sb) const {
    if (_k == 0 && (!_coeffs.has_value() || _coeffs->empty()))
        return sb;
    sb += _id;
    sb += "\t";
    sb += d(_radius);
    sb += "\t";
    sb += d(_k);
    sb += "\t";
    int start = 0;
    if (_asph_type == ASPH_EVEN)
        start = 1;
    else if (_asph_type == ASPH_ODD)
        start = 2;
    for (int j = start; j < static_cast<int>(_coeffs->size()); j++) {
        sb += d((*_coeffs)[static_cast<std::size_t>(j)]);
        sb += "\t";
    }
    sb += "\n";
    return sb;
}

std::string &SurfaceType::aspheric_markdown_table_header(std::string &sb,
                                                         int max_coeffs) {
    sb += "## Aspherical Data";
    sb += "\n";
    sb += "| ID  | Type | k   |";
    for (int j = 0; j < max_coeffs; j++) {
        sb += " P" + i(j + 1) + " |";
    }
    sb += "\n";
    sb += "| --- | --- | --- |";
    for (int j = 0; j < max_coeffs; j++)
        sb += " --- |";
    sb += "\n";
    return sb;
}

std::string SurfaceType::asphere_type() const {
    switch (_asph_type) {
    case ASPH_EVEN:
    case ASPH_EVEN_A2:
        return "EVEN";
    case ASPH_ODD:
        return "ODD";
    }
    return "";
}

std::string &SurfaceType::asherics_to_markdown_table_row(std::string &sb,
                                                         int max_coeffs) const {
    sb += "| " + _id;
    sb += "| " + asphere_type();
    sb += " | " + d(_k);
    for (int j = 0; j < max_coeffs; j++) {
        if (_coeffs.has_value() && j < static_cast<int>(_coeffs->size()))
            sb += " | " + d((*_coeffs)[static_cast<std::size_t>(j)]);
        else
            sb += " | 0 ";
    }
    sb += " |\n";
    return sb;
}

std::string &SurfaceType::variables_to_markdown_table_row(std::string &sb) const {
    if (_diameter_by_scenario.has_value()) {
        sb += "| a" + _id + " |";
        for (double v : *_diameter_by_scenario)
            sb += " " + d(v) + " |";
        sb += "\n";
    }
    if (_thickness_by_scenario.has_value()) {
        sb += "| d" + _id + " |";
        for (double v : *_thickness_by_scenario)
            sb += " " + d(v) + " |";
        sb += "\n";
    }
    return sb;
}

std::string &SurfaceType::to_markdown_table_header(std::string &sb) {
    sb += "## Surface Data";
    sb += "\n";
    sb += "Note that where glass types are shown the refractive index and abbe number "
          "is as per assigned glass type\n\n";
    sb += "| ID  | Radius | Thickness | Diameter | nd  | vd  | Glass Make | Glass |\n";
    sb += "| --- | ---    | ---       | ---      | --- | --- | ---        | ---   |\n";
    return sb;
}

std::string &SurfaceType::to_markdown_table_row(std::string &sb) const {
    sb += "| " + _id + " | ";
    if (_is_aperture_stop)
        sb += "AS";
    else if (_is_field_stop)
        sb += "FS";
    else
        sb += d(_radius);
    sb += " | ";
    if (_thickness_by_scenario.has_value())
        sb += "d" + _id;
    else
        sb += d(_thickness);
    sb += " | ";
    if (_diameter_by_scenario.has_value())
        sb += "a" + _id;
    else
        sb += d(_diameter);
    sb += " | ";
    std::string glassMaker;
    double nd = _nd;
    double vd = _vd;
    if (nd != 0.0 && _glass_name.has_value()) {
        auto glass = Glass::glass_by_catalog_name(_catalog_name, *_glass_name);
        if (glass != nullptr) {
            nd = glass->nd;
            vd = glass->vd;
            glassMaker = glass->catalog_name.value();
        }
    }
    if (nd != 0.0)
        sb += d(nd);
    sb += " | ";
    if (nd != 0.0)
        sb += d(vd);
    sb += " | ";
    if (nd != 0.0 && _glass_name.has_value()) {
        sb += glassMaker;
        sb += " | ";
        sb += *_glass_name;
    }
    sb += " |\n";
    return sb;
}

std::string SurfaceType::toString() const {
    std::string sb;
    to_opt_bench_str(sb, false);
    return sb;
}

} // namespace redukti::spec
