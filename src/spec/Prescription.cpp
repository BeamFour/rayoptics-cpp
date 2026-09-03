// C++ port of org.redukti.spec.Prescription
#include "redukti/spec/Prescription.h"

#include "redukti/Exceptions.h"
#include "redukti/Text.h"
#include "redukti/mathlib/M.h"
#include "redukti/rayoptics/elem/surface/Surface.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/seq/Glass.h"
#include "redukti/rayoptics/seq/SequentialModel.h"

#include <cmath>

namespace redukti::spec {

using importers::OpticalBenchDataImporter;
using rayoptics::seq::Glass;

namespace {
std::string d(double v) { return doubleToString(v); }
std::string i(int v) { return std::to_string(v); }
} // namespace

Prescription::Prescription(double focal_length, double fno,
                           double angle_of_view_degrees, double diameter_image_circle,
                           bool d_line_only)
    : Prescription(focal_length, fno, angle_of_view_degrees, diameter_image_circle,
                   d_line_only ? std::vector<double>{587.5618}
                               : std::vector<double>{587.5618, 486.1327, 656.2725},
                   d_line_only ? std::vector<double>{1.0}
                               : std::vector<double>{1.0, 1.0, 1.0}) {}

Prescription::Prescription(double focal_length, double fno,
                           double angle_of_view_degrees, double diameter_image_circle,
                           std::vector<double> wvls, std::vector<double> wts)
    : _focal_length(focal_length), _fno(fno),
      _angle_of_view_in_degrees(angle_of_view_degrees),
      _diameter_image_circle(diameter_image_circle), _wvls(std::move(wvls)),
      _wts(std::move(wts)) {}

Prescription &Prescription::surf(double radius, double thickness, double diameter,
                                 double nd, double vd,
                                 const std::optional<std::string> &glass_name,
                                 const std::optional<std::string> &catalog_name) {
    _surface_list.push_back(SurfaceType(i(static_cast<int>(_surface_list.size()) + 1),
                                        false, radius, thickness, diameter, nd, vd,
                                        glass_name, catalog_name));
    return *this;
}

Prescription &Prescription::surf(double radius, double thickness, double diameter,
                                 double nd, double vd) {
    _surface_list.push_back(SurfaceType(i(static_cast<int>(_surface_list.size()) + 1),
                                        false, radius, thickness, diameter, nd, vd,
                                        std::nullopt, std::nullopt));
    return *this;
}

Prescription &Prescription::surf(double radius, double thickness, double diameter) {
    _surface_list.push_back(SurfaceType(i(static_cast<int>(_surface_list.size()) + 1),
                                        false, radius, thickness, diameter, 0, 0,
                                        std::nullopt, std::nullopt));
    return *this;
}

Prescription &Prescription::stop(double thickness, double diameter) {
    _surface_list.push_back(SurfaceType(i(static_cast<int>(_surface_list.size()) + 1),
                                        true, 0, thickness, diameter, 0, 0, std::nullopt,
                                        std::nullopt));
    return *this;
}

Prescription &Prescription::field_stop(double thickness, double diameter) {
    SurfaceType surface(i(static_cast<int>(_surface_list.size()) + 1), false, 0, thickness,
                        diameter, 0, 0, std::nullopt, std::nullopt);
    surface._is_field_stop = true;
    _surface_list.push_back(surface);
    return *this;
}

Prescription &Prescription::asph(int asph_type, double k,
                                 const std::vector<double> &coeffs) {
    if (asph_type == SurfaceType::ASPH_EVEN && coeffs[0] != 0.0)
        throw IllegalArgumentException("EVEN aspheres must have 0 as first coefficient");
    else if (asph_type == SurfaceType::ASPH_ODD && (coeffs[0] != 0.0 || coeffs[1] != 0.0))
        throw IllegalArgumentException(
            "ODD aspheres must have 0 as first and second coefficients");
    SurfaceType &lastSurface = _surface_list.back();
    lastSurface._asph_type = asph_type;
    lastSurface._k = k;
    lastSurface._coeffs = coeffs;
    return *this;
}

double Prescription::image_diameter_for_field(double field) const {
    return _diameter_image_circle * field;
}

double Prescription::full_angle_of_view_degrees(double field) const {
    if (field == 0.0)
        return 0.0;
    if (field > 0 && field <= 1.0) {
        auto radius = image_diameter_for_field(field) / 2.0;
        auto radians = std::atan(radius / _focal_length);
        return 2.0 * mathlib::M::toDegrees(radians);
    }
    throw IllegalArgumentException("Field must be between 0 and 1.");
}

double Prescription::get_half_angle_of_view_in_radians() const {
    return mathlib::M::toRadians(_angle_of_view_in_degrees / 2.0);
}

Prescription &Prescription::build() {
    _built = true;
    return *this;
}

Prescription &Prescription::import_surface(
    const OpticalBenchDataImporter::LensSurface &surface, int scenario,
    bool use_glass_types) {
    double thickness = surface.get_thickness(scenario);
    double radius = surface.get_radius();
    double refractive_index = surface.get_refractive_index();
    double abbe_vd = surface.get_abbe_vd();
    double diameter = surface.get_diameter(scenario);
    const auto &glass_name = surface.get_glass_name();
    const auto &catalog_name = surface.get_catalog_name();
    if (surface.get_surface_type() ==
        OpticalBenchDataImporter::SurfaceType::aperture_stop) {
        stop(thickness, diameter);
        return *this;
    } else if (surface.get_surface_type() ==
               OpticalBenchDataImporter::SurfaceType::field_stop) {
        field_stop(thickness, diameter);
        return *this;
    }
    if (use_glass_types && glass_name.has_value() &&
        Glass::glass_by_catalog_name(catalog_name, *glass_name) != nullptr) {
        surf(radius, thickness, diameter, refractive_index, abbe_vd, glass_name,
             catalog_name);
    } else if (refractive_index != 0.0) {
        surf(radius, thickness, diameter, refractive_index, abbe_vd);
    } else {
        surf(radius, thickness, diameter);
    }
    const OpticalBenchDataImporter::AsphericalData *aspherical_data =
        surface.get_aspherical_data();
    if (aspherical_data != nullptr) {
        int asph_type = SurfaceType::ASPH_EVEN;
        switch (aspherical_data->get_asphere_type()) {
        case OpticalBenchDataImporter::AsphereType::Even:
            asph_type = SurfaceType::ASPH_EVEN;
            break;
        case OpticalBenchDataImporter::AsphereType::EvenA2:
            asph_type = SurfaceType::ASPH_EVEN_A2;
            break;
        case OpticalBenchDataImporter::AsphereType::Odd:
            asph_type = SurfaceType::ASPH_ODD;
            break;
        }
        double k = aspherical_data->get_cc();
        std::vector<double> coeffs = aspherical_data->get_coeffs();
        asph(asph_type, k, coeffs);
    }
    return *this;
}

Prescription Prescription::build_prescription_d_line(const LensSpecifications &specs) {
    return build_prescription(specs, true, false, true);
}

Prescription Prescription::build_prescription_e_line(const LensSpecifications &specs) {
    return build_prescription(specs, true, std::vector<double>{546.074},
                              std::vector<double>{1.0}, 0);
}

Prescription Prescription::build_prescription(const LensSpecifications &specs,
                                              bool use_glass_types) {
    return build_prescription(specs, use_glass_types, false, false);
}

Prescription Prescription::build_prescription(const LensSpecifications &specs,
                                              bool use_glass_types,
                                              const std::vector<double> &wvls,
                                              const std::vector<double> &wts) {
    return build_prescription(specs, use_glass_types, wvls, wts, 0);
}

Prescription Prescription::build_prescription(const LensSpecifications &specs,
                                              bool use_glass_types, bool weighted,
                                              bool d_line) {
    std::vector<double> wvls;
    std::vector<double> wts;
    if (d_line) {
        wvls = {Glass::d};
        wts = {1.0};
    } else if (weighted) {
        wvls = {Glass::d, Glass::C, Glass::e, Glass::F, Glass::g};
        wts = {WT_d, WT_C, WT_e, WT_F, WT_g};
    } else {
        wvls = {Glass::d, Glass::F, Glass::C};
        wts = {1.0, 1.0, 1.0};
    }
    return build_prescription(specs, use_glass_types, wvls, wts, 0);
}

namespace {

const char *const FOCAL_LENGTH = "Focal Length";
const char *const F_NUMBER = "F-Number";
const char *const ANGLE_OF_VIEW = "Angle of View";

double require_positive_value(const Prescription::LensSpecifications &specs,
                              const char *name, const char *expected, int scenario) {
    const auto *variable = specs.find_variable(name);
    if (variable == nullptr)
        throw IllegalArgumentException(
            std::string("The prescription does not specify '") + name +
            "'; add it to the [variable distances] section as " + expected);
    if (scenario >= variable->num_values())
        throw IllegalArgumentException(
            std::string("The prescription specifies '") + name + "' for " +
            std::to_string(variable->num_values()) + " scenario(s), but scenario " +
            std::to_string(scenario) + " was requested");
    double value = variable->get_value_as_double(scenario);
    if (!(value > 0.0))
        throw IllegalArgumentException(
            std::string("The prescription specifies '") + name + "' as '" +
            variable->get_value(scenario) + "' for scenario " + std::to_string(scenario) +
            "; expected " + expected + ", which must be positive");
    return value;
}

} // namespace

Prescription Prescription::build_prescription(const LensSpecifications &specs,
                                              bool use_glass_types,
                                              const std::vector<double> &wvls,
                                              const std::vector<double> &wts,
                                              int default_scenario) {
    Prescription prescription(
        require_positive_value(specs, FOCAL_LENGTH, "the focal length in mm",
                               default_scenario),
        require_positive_value(specs, F_NUMBER, "the f-number", default_scenario),
        require_positive_value(specs, ANGLE_OF_VIEW,
                               "the full angle of view in degrees", default_scenario),
        specs.get_image_height(), wvls, wts);
    prescription._title = specs.get_descriptive_data().get_value("title");
    const auto &patent_info_n = specs.get_patent_info();
    if (patent_info_n.count() > 0) {
        prescription._patent_country = patent_info_n.get_value("country");
        prescription._patent_number = patent_info_n.get_value("number");
        prescription._patent_example = patent_info_n.get_value("example");
        prescription._application_year = patent_info_n.get_value("year applied");
        prescription._inventors = patent_info_n.get_value("inventors");
        prescription._original_assignee = patent_info_n.get_value("original assignee");
        prescription._current_assignee = patent_info_n.get_value("current assignee");
        prescription._patent_link = patent_info_n.get_value("link");
    } else if (specs.get_descriptive_data().find_variable("patent") != nullptr) {
        const auto *patentInfo = specs.get_descriptive_data().find_variable("patent");
        prescription._patent_country = patentInfo->get_value(0);
        prescription._patent_number = patentInfo->get_value(1);
        prescription._patent_example = patentInfo->get_value(2);
        prescription._application_year = patentInfo->get_value(3);
        prescription._inventors = patentInfo->get_value(4);
        prescription._current_assignee = patentInfo->get_value(5);
        prescription._original_assignee = patentInfo->get_value(5);
        prescription._patent_link = patentInfo->get_value(6);
    }
    const auto &report_data = specs.get_report_data();
    std::optional<std::string> lensName;
    if (report_data.count() > 0)
        lensName = report_data.get_value("lens name");
    if (!lensName.has_value()) {
        const auto *variable = specs.get_descriptive_data().find_variable("lens name");
        if (variable != nullptr)
            lensName = variable->get_value(0);
    }
    if (lensName.has_value())
        prescription._lens_name = lensName;
    prescription.add_configurations(specs);
    const auto &surfaces = specs.get_surfaces();
    for (std::size_t k = 0; k < surfaces.size(); k++) {
        prescription.import_surface(surfaces[k], default_scenario, use_glass_types);
        if (prescription._configurations.has_value())
            prescription.add_configuration_data(surfaces[k]);
    }
    return prescription.build();
}

void Prescription::add_configuration_data(
    const OpticalBenchDataImporter::LensSurface &lensSurface) {
    if (!_configurations.has_value() || _configurations->empty())
        return;
    SurfaceType &lastSurface = _surface_list.back();
    const auto &thickness_by_scenario = lensSurface.get_thickness_by_scenario();
    if (thickness_by_scenario.size() > 1) {
        std::vector<double> thickness(_configurations->size(), 0.0);
        for (std::size_t k = 0; k < _configurations->size(); k++) {
            int scenario = (*_configurations)[k];
            thickness[k] = thickness_by_scenario[static_cast<std::size_t>(scenario)];
        }
        lastSurface.set_thickness_by_scenario(thickness);
    }
    if (lensSurface.is_aperture_stop()) {
        const auto &diameter_by_scenario = lensSurface.get_diameter_by_scenario();
        if (diameter_by_scenario.size() > 1) {
            std::vector<double> diameter(_configurations->size(), 0.0);
            for (std::size_t k = 0; k < diameter.size(); k++) {
                int scenario = (*_configurations)[k];
                diameter[k] = diameter_by_scenario[static_cast<std::size_t>(scenario)];
            }
            lastSurface.set_diameter_by_scenario(diameter);
        }
    }
}

Prescription &Prescription::add_configurations(const LensSpecifications &specs) {
    const auto *configurations = specs.get_report_data().find_variable("scenarios");
    const auto *configuration_names = specs.get_report_data().find_variable("names");
    if (configurations != nullptr && configurations->num_values() > 0 &&
        configuration_names != nullptr &&
        configuration_names->num_values() == configurations->num_values()) {
        std::vector<int> configs(static_cast<std::size_t>(configurations->num_values()));
        std::vector<std::string> names(
            static_cast<std::size_t>(configuration_names->num_values()));
        for (int k = 0; k < configurations->num_values(); k++) {
            configs[static_cast<std::size_t>(k)] =
                configurations->get_value_as_integer(k, 0);
            names[static_cast<std::size_t>(k)] = configuration_names->get_value(k);
        }
        _configurations = configs;
        _configuration_names = names;
        _focal_length_by_scenario.assign(configs.size(), 0.0);
        _f_number_by_scenario.assign(configs.size(), 0.0);
        _angle_of_views_by_scenario.assign(configs.size(), 0.0);
        for (std::size_t k = 0; k < configs.size(); k++) {
            int scenario = configs[k];
            _focal_length_by_scenario[k] =
                require_positive_value(specs, FOCAL_LENGTH, "the focal length in mm",
                                       scenario);
            _f_number_by_scenario[k] =
                require_positive_value(specs, F_NUMBER, "the f-number", scenario);
            _angle_of_views_by_scenario[k] = require_positive_value(
                specs, ANGLE_OF_VIEW, "the full angle of view in degrees", scenario);
        }
    }
    return *this;
}

namespace {

/**
 * Java's `Math.round(value * scale) / scale`. Math.round is floor(x + 0.5) --
 * it rounds a half *up*, not to even, so std::round is not the same function
 * for negative halves and llround is not either.
 */
double round_to(double value, int decimals) {
    double scale = std::pow(10, decimals);
    return static_cast<double>(
               static_cast<long long>(std::floor(value * scale + 0.5))) /
           scale;
}

} // namespace

int Prescription::update_apertures_from(rayoptics::optical::OpticalModel *opm, int config,
                                        const std::vector<int> *avoid_list,
                                        const std::vector<int> *include_list,
                                        int decimals) {
    auto &ifcs = opm->seq_model->ifcs;
    const int n = static_cast<int>(_surface_list.size());
    if (static_cast<int>(ifcs.size()) != n + 2)
        throw IllegalArgumentException(
            "Model has " + std::to_string(ifcs.size()) + " interfaces, expected " +
            std::to_string(n + 2) + " for this prescription's " + std::to_string(n) +
            " surfaces. The model must be built from this prescription.");
    std::vector<int> targets;
    if (avoid_list != nullptr) {
        for (int k = 0; k < n; k++) {
            if (std::find(avoid_list->begin(), avoid_list->end(), k) == avoid_list->end())
                targets.push_back(k);
        }
    } else if (include_list != nullptr) {
        targets = *include_list;
    } else {
        for (int k = 0; k < n; k++)
            targets.push_back(k);
    }
    int changed = 0;
    for (int k : targets) {
        if (k < 0 || k >= n)
            throw IllegalArgumentException("Surface index out of range: " +
                                           std::to_string(k));
        SurfaceType &surface = _surface_list[static_cast<std::size_t>(k)];
        double diameter =
            round_to(ifcs[static_cast<std::size_t>(k + 1)]->surface_od() * 2.0, decimals);
        bool modified = false;
        if (surface._diameter_by_scenario.has_value()) {
            if (config < 0 ||
                config >= static_cast<int>(surface._diameter_by_scenario->size()))
                throw IllegalArgumentException(
                    "Config " + std::to_string(config) + " out of range for surface " +
                    std::to_string(k) + ", which has " +
                    std::to_string(surface._diameter_by_scenario->size()) +
                    " configurations");
            if ((*surface._diameter_by_scenario)[static_cast<std::size_t>(config)] !=
                diameter) {
                (*surface._diameter_by_scenario)[static_cast<std::size_t>(config)] =
                    diameter;
                modified = true;
            }
        }
        if ((!surface._diameter_by_scenario.has_value() || config == 0) &&
            surface._diameter != diameter) {
            surface._diameter = diameter;
            modified = true;
        }
        if (modified)
            changed++;
    }
    return changed;
}

bool Prescription::has_odd_aspheric() const {
    for (const auto &s : _surface_list) {
        if (s.is_odd_asphere())
            return true;
    }
    return false;
}

bool Prescription::has_even_a2_aspheric() const {
    for (const auto &s : _surface_list) {
        if (s.is_even_a2_asphere())
            return true;
    }
    return false;
}

void Prescription::add_patent_section(std::string &sb) const {
    sb += "[patent info]\n";
    if (!_patent_country.empty())
        sb += "country\t" + _patent_country + "\n";
    if (_patent_number.has_value() && !_patent_number->empty())
        sb += "number\t" + *_patent_number + "\n";
    if (!_patent_example.empty())
        sb += "example\t" + _patent_example + "\n";
    if (!_application_year.empty())
        sb += "year applied\t" + _application_year + "\n";
    if (!_inventors.empty())
        sb += "inventors\t" + _inventors + "\n";
    if (!_current_assignee.empty())
        sb += "current assignee\t" + _current_assignee + "\n";
    if (!_original_assignee.empty())
        sb += "original assignee\t" + _original_assignee + "\n";
    if (!_patent_link.empty())
        sb += "link\t" + _patent_link + "\n";
}

void Prescription::add_report_section(std::string &sb) const {
    sb += "[report data]\n";
    if (_lens_name.has_value() && !_lens_name->empty())
        sb += "lens name\t" + *_lens_name + "\n";
    if (_configurations.has_value()) {
        sb += "scenarios";
        for (std::size_t k = 0; k < _configurations->size(); k++)
            sb += "\t" + i(static_cast<int>(k));
        sb += "\n";
        sb += "names";
        for (const auto &nm : *_configuration_names)
            sb += "\t" + nm;
        sb += "\n";
    }
}

std::string &Prescription::to_opt_bench_str(std::string &sb) const {
    sb += "[descriptive data]\n";
    sb += "title\t" + _title + "\n";
    sb += "[constants]\n";
    if (has_odd_aspheric())
        sb += "AsphericalOddCount\t1\n";
    else if (has_even_a2_aspheric())
        sb += "AsphericalA2\n";
    sb += "[variable distances]\n";
    sb += "Focal Length";
    if (!_configurations.has_value())
        sb += "\t" + d(_focal_length);
    else
        for (double v : _focal_length_by_scenario)
            sb += "\t" + d(v);
    sb += "\n";
    sb += "Angle of View";
    if (!_configurations.has_value())
        sb += "\t" + d(_angle_of_view_in_degrees);
    else
        for (double v : _angle_of_views_by_scenario)
            sb += "\t" + d(v);
    sb += "\n";
    sb += "F-Number";
    if (!_configurations.has_value())
        sb += "\t" + d(_fno);
    else
        for (double v : _f_number_by_scenario)
            sb += "\t" + d(v);
    sb += "\n";
    sb += "Image Height";
    if (!_configurations.has_value())
        sb += "\t" + d(_diameter_image_circle);
    else
        for (std::size_t k = 0; k < _configurations->size(); k++)
            sb += "\t" + d(_diameter_image_circle);
    sb += "\n";
    sb += "Magnification";
    if (!_configurations.has_value())
        sb += "\t0";
    else
        for (std::size_t k = 0; k < _configurations->size(); k++)
            sb += "\t0";
    sb += "\n";
    const SurfaceType &last_surf = _surface_list.back();
    sb += "Bf";
    if (!_configurations.has_value()) {
        sb += "\t" + d(last_surf._thickness);
    } else {
        if (!last_surf._thickness_by_scenario.has_value()) {
            for (std::size_t k = 0; k < _configurations->size(); k++)
                sb += "\t" + d(last_surf._thickness);
        } else {
            for (double v : *last_surf._thickness_by_scenario)
                sb += "\t" + d(v);
        }
    }
    sb += "\n";
    for (std::size_t k = 0; k < _surface_list.size(); k++) {
        const SurfaceType &surf_k = _surface_list[k];
        if (k < _surface_list.size() - 1 && surf_k._thickness_by_scenario.has_value()) {
            sb += "d" + surf_k._id;
            for (double v : *surf_k._thickness_by_scenario)
                sb += "\t" + d(v);
            sb += "\n";
        }
        if (surf_k._diameter_by_scenario.has_value() && surf_k._is_aperture_stop) {
            sb += "Aperture Diameter";
            for (double v : *surf_k._diameter_by_scenario)
                sb += "\t" + d(v);
            sb += "\n";
        }
    }
    sb += "[lens data]\n";
    for (std::size_t k = 0; k < _surface_list.size(); k++) {
        bool is_last = k == _surface_list.size() - 1;
        _surface_list[k].to_opt_bench_str(sb, is_last);
    }
    sb += "[aspherical data]\n";
    for (const auto &surface : _surface_list)
        surface.aspherics_to_opt_bench_str(sb);
    sb += "[notes]\n";
    sb += "Generated by Beam42\n";
    add_patent_section(sb);
    add_report_section(sb);
    return sb;
}

std::string Prescription::organization() const {
    if (!_original_assignee.empty())
        return _original_assignee;
    return _current_assignee;
}

std::string &Prescription::to_markdown_str(std::string &sb) const {
    if (_lens_name.has_value())
        sb += "# " + *_lens_name + "\n";
    if (_patent_number.has_value()) {
        sb += "## Patent Information\n";
        sb += "| Country | Patent Number | Example | Year of Application | Inventors | "
              "Organisation | Link |\n";
        sb += "| ---     | ---           | ---     | ---                 | ---       | "
              "---          | ---  |\n";
        sb += "|" + _patent_country + " | " + *_patent_number + " | " + _patent_example +
              " | " + _application_year + " | " + _inventors + " | " + organization() +
              " | " + "[link](" + _patent_link + ") |\n";
    }
    bool sawAsph = false;
    SurfaceType::to_markdown_table_header(sb);
    for (const auto &surface : _surface_list) {
        surface.to_markdown_table_row(sb);
        if (surface.is_aspheric())
            sawAsph = true;
    }
    if (sawAsph) {
        int max_coeffs = 0;
        for (const auto &surface : _surface_list) {
            if (surface._coeffs.has_value() &&
                static_cast<int>(surface._coeffs->size()) > max_coeffs)
                max_coeffs = static_cast<int>(surface._coeffs->size());
        }
        SurfaceType::aspheric_markdown_table_header(sb, max_coeffs);
        for (const auto &surface : _surface_list) {
            if (surface.is_aspheric())
                surface.asherics_to_markdown_table_row(sb, max_coeffs);
        }
    }
    if (get_num_configurations() > 1) {
        sb += "## Variables\n";
        sb += "| Variable |";
        for (int k = 0; k < get_num_configurations(); k++)
            sb += " " + (*_configuration_names)[static_cast<std::size_t>(k)] + " |";
        sb += "\n";
        sb += "| --- |";
        for (int k = 0; k < get_num_configurations(); k++)
            sb += " --- |";
        sb += "\n";
        sb += "| Focal length |";
        for (int k = 0; k < get_num_configurations(); k++)
            sb += d(_focal_length_by_scenario[static_cast<std::size_t>(k)]) + " |";
        sb += "\n";
        sb += "| F-Number |";
        for (int k = 0; k < get_num_configurations(); k++)
            sb += d(_f_number_by_scenario[static_cast<std::size_t>(k)]) + " |";
        sb += "\n";
        sb += "| Angle of View |";
        for (int k = 0; k < get_num_configurations(); k++)
            sb += d(_angle_of_views_by_scenario[static_cast<std::size_t>(k)]) + " |";
        sb += "\n";
        for (const auto &surface : _surface_list)
            surface.variables_to_markdown_table_row(sb);
    }
    return sb;
}

std::vector<std::pair<double, double>> Prescription::get_wvl_wts() const {
    std::vector<std::pair<double, double>> out;
    for (std::size_t k = 0; k < _wvls.size(); k++)
        out.push_back({_wvls[k], _wts[k]});
    return out;
}

std::string Prescription::toString() const {
    std::string sb;
    to_opt_bench_str(sb);
    return sb;
}

} // namespace redukti::spec
