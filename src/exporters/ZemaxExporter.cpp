// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.exporters.ZemaxExporter
#include "redukti/exporters/ZemaxExporter.h"

#include "redukti/Exceptions.h"
#include "redukti/Text.h"

namespace redukti::exporters {

using spec::Prescription;
using spec::SurfaceType;

namespace {

/** Java StringBuilder.append(double) is Double.toString. */
std::string d(double v) { return doubleToString(v); }
/** Java StringBuilder.append(int). */
std::string i(int v) { return std::to_string(v); }

} // namespace

std::string ZemaxExporter::generate(const Prescription &prescription,
                                    bool d_line_only) const {
    std::string sb;
    outputHeading(prescription, d_line_only, sb);
    output_object(prescription, sb);
    output_surfaces(prescription, sb);
    output_image_plane(prescription, sb);
    output_configurations(prescription, sb);
    return sb;
}

void ZemaxExporter::outputHeading(const Prescription &prescription, bool d_line_only,
                                  std::string &sb) {
    sb += "VERS 161019 507 33785\n";
    sb += "MODE SEQ\n";
    sb += "NAME " + prescription.get_title() + "\n";
    sb += "PFIL 0 0 0\nLANG 0\nUNIT MM X W X CM MR CPMM\n";
    if (prescription.get_num_configurations() == 0)
        sb += "FNUM " + d(prescription.get_f_number()) + "\n";
    else
        sb += "FLOA\n";
    double img_ht = prescription._diameter_image_circle / 2.0;
    // The Java uses a text block here, which strips the common indentation, so
    // the emitted lines have no leading whitespace.
    sb += "ENVD 20 1 0\n";
    sb += "GFAC 0 0\n";
    sb += "GCAT SCHOTT HOYA OHARA NIKON-HIKARI NIKON HIKARI SUMITA CDGM CORNING "
          "LACROIX\n";
    sb += "RAIM 0 2 1 1 0 1 0 0 0\n";
    sb += "SDMA 0 1 0\n";
    if (!prescription._vignetting.has_value()) {
        sb += "FTYP 3 0 9 3 0 0 0 9\n";
        output_rop_lines(sb);
        sb += "XFLN 0 0 0 0 0 0 0 0 0 0 0 0\n";
        sb += "YFLN 0 2.16 4.33 6.5 8.65 10.8 12.98 15.14 " + formatF(img_ht, 6) +
              " 0 0 0 0 0 0 0 0\n";
        sb += "FWGN 10 9 9 8 8 7 7 7 3 1 1 1\n";
        sb += "VDXN 0 0 0 0 0 0 0 0 0 0 0 0\n";
        sb += "VDYN 0 0 0 0 0 0 0 0 0 0 0 0\n";
        sb += "VCXN 0 0 0 0 0 0 0 0 0 0 0 0\n";
        sb += "VCYN 0 0 0 0 0 0 0 0 0 0 0 0\n";
        sb += "VANN 0 0 0 0 0 0 0 0 0 0 0 0\n";
    } else
        output_vignetted_fields(prescription, sb);
    if (d_line_only) {
        sb += "WAVM 1 0.5875618 1\n";
        sb += "WAVM 2 0.550 0\n";
        sb += "WAVM 3 0.550 0\n";
    } else {
        sb += "WAVM 1 0.4861327 1\n";
        sb += "WAVM 2 0.5875618 1\n";
        sb += "WAVM 3 0.6562725 1\n";
    }
    for (int w = 4; w <= 24; w++)
        sb += "WAVM " + i(w) + " 0.550 0\n";
    sb += "PWAV 2\n";
    sb += "POLS 1 0 1 0 0 1 0\n";
    sb += "GSTD 0 100 100 100 100 100 100 0 1 1 0 0 1 1 1 1 1 1\n";
    sb += "NSCD 100 500 0 1.0E-3 5 1.0E-6 0 0 0 0 0 0 1000000 0 2\n";
}

void ZemaxExporter::output_rop_lines(std::string &sb) {
    sb += "ROPD 2\n";
    sb += "HYPR 0\n";
    sb += "PICB 1\n";
}

void ZemaxExporter::output_vignetted_fields(const Prescription &prescription,
                                            std::string &sb) {
    auto fields = vignetting_by_config(prescription, 0);
    std::size_t n = fields.size();
    double half_angle = half_angle_of_view(prescription, 0);
    sb += "FTYP 0 0 " + i(static_cast<int>(n)) + " 3 0 0 0 " + i(static_cast<int>(n)) + "\n";
    output_rop_lines(sb);
    std::vector<double> x(n), y(n), weights(n), vdx(n), vdy(n), vcx(n), vcy(n);
    for (std::size_t k = 0; k < n; k++) {
        const auto &v = fields[k];
        y[k] = v.field * half_angle;
        weights[k] = 1;
        vdx[k] = decenter(v.vux, v.vlx);
        vdy[k] = decenter(v.vuy, v.vly);
        vcx[k] = compression(v.vux, v.vlx);
        vcy[k] = compression(v.vuy, v.vly);
    }
    field_line(sb, "XFLN", x);
    field_line(sb, "YFLN", y);
    field_line(sb, "FWGN", weights);
    field_line(sb, "VDXN", vdx);
    field_line(sb, "VDYN", vdy);
    field_line(sb, "VCXN", vcx);
    field_line(sb, "VCYN", vcy);
    field_line(sb, "VANN", {});
}

void ZemaxExporter::field_line(std::string &sb, const std::string &key,
                               const std::vector<double> &values) {
    sb += key;
    // Java's `i < values.length ? values[i] : 0` promotes the 0 to a double, so the
    // padding is written as 0.0.
    for (std::size_t k = 0; k < static_cast<std::size_t>(MAX_FIELDS); k++)
        sb += " " + d(k < values.size() ? values[k] : 0.0);
    sb += "\n";
}

double ZemaxExporter::decenter(double upper, double lower) { return (lower - upper) / 2.0; }

double ZemaxExporter::compression(double upper, double lower) { return (upper + lower) / 2.0; }

double ZemaxExporter::half_angle_of_view(const Prescription &prescription, int config) {
    double angle_of_view =
        config == 0 ? prescription._angle_of_view_in_degrees
                    : prescription._angle_of_views_by_scenario[static_cast<std::size_t>(config)];
    return angle_of_view / 2.0;
}

std::vector<Prescription::FieldVignetting> ZemaxExporter::vignetting_by_config(
    const Prescription &prescription, int config) {
    std::vector<Prescription::FieldVignetting> rows;
    std::vector<double> first;
    for (const auto &v : *prescription._vignetting) {
        if (v.config == config)
            rows.push_back(v);
        if (v.config == 0)
            first.push_back(v.field);
    }
    if (rows.empty() || rows.size() > static_cast<std::size_t>(MAX_FIELDS))
        throw IllegalArgumentException("Configuration " + i(config) + " has " +
                                       i(static_cast<int>(rows.size())) +
                                       " vignetted fields; Zemax needs 1 to " + i(MAX_FIELDS));
    if (rows.size() != first.size())
        throw IllegalArgumentException("Configuration " + i(config) +
                                       " has different fields from the first configuration");
    for (std::size_t k = 0; k < rows.size(); k++)
        if (rows[k].field != first[k])
            throw IllegalArgumentException("Configuration " + i(config) +
                                           " has different fields from the first configuration");
    return rows;
}

void ZemaxExporter::output_object(const Prescription &prescription, std::string &sb) {
    (void)prescription;
    sb += "SURF 0\n";
    sb += "  TYPE STANDARD\n";
    sb += "  CURV 0.0 0 0 0 0 \"\"\n";
    sb += "  HIDE 0 0 0 0 0 0 0 0 0 0\n";
    sb += "  MIRR 2 0\n";
    sb += "  DISZ INFINITY\n";
    sb += "  DIAM 0 1 0 0 1 \"\"\n";
    sb += "  POPS 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0\n";
}

void ZemaxExporter::output_surfaces(const Prescription &prescription, std::string &sb) {
    const auto &surfaces = prescription.get_surfaces();
    for (std::size_t k = 0; k < surfaces.size(); k++) {
        const SurfaceType &s = surfaces[k];
        double thickness = 0.0;
        double diameter = s.get_diameter();
        diameter /= 2.0;
        thickness += s.get_thickness();
        sb += "SURF " + i(static_cast<int>(k) + 1) + "\n";
        if (s.is_aperture_stop())
            sb += "  STOP\n";
        if (s.is_aspheric()) {
            if (s.is_radial_asphere())
                sb += "  TYPE XOSPHERE\n";
            else
                sb += "  TYPE EVENASPH\n";
        } else {
            sb += "  TYPE STANDARD\n";
        }
        double curvature =
            s.get_radius_of_curvature() == 0.0 ? 0 : 1.0 / s.get_radius_of_curvature();
        sb += "  CURV " + d(curvature) + " 0 0 0 0\n";
        sb += "  HIDE 0 0 0 0 0 0 0 0 0 0\n";
        sb += "  MIRR 2 0\n";
        if (s.is_aspheric()) {
            std::vector<double> aspherics = s.get_aspheric_coeffs();
            if (s.is_radial_asphere()) {
                int count = static_cast<int>(aspherics.size());
                if (count > 240)
                    throw IllegalArgumentException(
                        "XOSPHERE supports at most 240 radial terms");
                // Extra Data: term count, normalization radius, then A1, A2, ... .
                // A unit normalization radius preserves our dimensional coefficients.
                sb += "  XDAT 1 " + i(count) + "\n";
                sb += "  XDAT 2 1.0\n";
                for (int a = 0; a < count; a++)
                    sb += "  XDAT " + i(a + 3) + " " +
                          d(aspherics[static_cast<std::size_t>(a)]) + "\n";
            } else {
                for (std::size_t a = 1; a <= aspherics.size(); a++) {
                    sb += "  PARM " + i(static_cast<int>(a)) + " ";
                    sb += d(aspherics[a - 1]) + "\n";
                }
            }
        }
        sb += "  DISZ " + d(thickness) + "\n";
        if (s.is_aspheric()) {
            sb += "  CONI " + d(s.get_cc()) + "\n";
        }
        if (s.get_refractive_index() != 0.0) {
            sb += "  GLAS ";
            if (s.get_glass_name().has_value())
                sb += *s.get_glass_name() + " 0 0 ";
            else
                sb += "___BLANK 1 0 ";
            sb += d(s.get_refractive_index()) + " " + d(s.get_abbe_vd()) +
                  " 0 0 0 0 0 0\n";
        }
        sb += "  DIAM " + d(diameter) + " 1 0 0 1 \"\"\n";
        if (s.is_field_stop())
            sb += "  CLAP 0 " + d(diameter) + " 0\n";
        sb += "  POPS 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0\n";
    }
}

void ZemaxExporter::output_image_plane(const Prescription &prescription,
                                       std::string &sb) {
    int sid = static_cast<int>(prescription.get_surfaces().size()) + 1;
    sb += "SURF " + i(sid) + "\n";
    double img_ht = prescription._diameter_image_circle / 2.0;
    sb += "  TYPE STANDARD\n";
    sb += "  CURV 0.0 0 0 0 0 \"\"\n";
    sb += "  HIDE 0 0 0 0 0 0 0 0 0 0\n";
    sb += "  MIRR 2 0\n";
    sb += "  DISZ 0\n";
    sb += "  DIAM " + formatF(img_ht, 6) + " 1 0 0 1 \"\"\n";
    sb += "  POPS 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0\n";
    sb += "TOL TOFF   0   0              0              0   0 0 0 0\n";
}

void ZemaxExporter::output_configurations(const Prescription &prescription,
                                          std::string &sb) {
    if (prescription.get_num_configurations() <= 1)
        return;
    sb += "MNUM " + i(prescription.get_num_configurations()) + " 1\n";
    const auto &surfaces = prescription.get_surfaces();
    for (std::size_t k = 0; k < surfaces.size(); k++) {
        const SurfaceType &surface = surfaces[k];
        if (surface._diameter_by_scenario.has_value()) {
            const auto &v = *surface._diameter_by_scenario;
            for (std::size_t j = 0; j < v.size(); j++)
                sb += "SDIA    " + i(static_cast<int>(k) + 1) + "   " +
                      i(static_cast<int>(j) + 1) + " " + d(v[j] / 2.0) +
                      "  0 0 0 1 1 1 0 0\n";
        }
        if (surface._thickness_by_scenario.has_value()) {
            const auto &v = *surface._thickness_by_scenario;
            for (std::size_t j = 0; j < v.size(); j++)
                sb += "THIC    " + i(static_cast<int>(k) + 1) + "   " +
                      i(static_cast<int>(j) + 1) + " " + d(v[j]) +
                      "  0 0 0 1 1 1 0 0\n";
        }
    }
    if (prescription._vignetting.has_value())
        output_vignetting_configurations(prescription, sb);
}

void ZemaxExporter::output_vignetting_configurations(const Prescription &prescription,
                                                     std::string &sb) {
    int num_configs = prescription.get_num_configurations();
    std::vector<std::vector<Prescription::FieldVignetting>> by_config;
    for (int config = 0; config < num_configs; config++)
        by_config.push_back(vignetting_by_config(prescription, config));
    std::size_t num_fields = by_config[0].size();
    for (std::size_t fi = 0; fi < num_fields; fi++) {
        for (const char *op : {"YFIE", "FVDX", "FVDY", "FVCX", "FVCY"}) {
            for (int config = 0; config < num_configs; config++) {
                const auto &v = by_config[static_cast<std::size_t>(config)][fi];
                std::string name = op;
                double value;
                if (name == "YFIE")
                    value = v.field * half_angle_of_view(prescription, config);
                else if (name == "FVDX")
                    value = decenter(v.vux, v.vlx);
                else if (name == "FVDY")
                    value = decenter(v.vuy, v.vly);
                else if (name == "FVCX")
                    value = compression(v.vux, v.vlx);
                else
                    value = compression(v.vuy, v.vly);
                sb += name + "    " + i(static_cast<int>(fi) + 1) + "   " + i(config + 1) + " " +
                      d(value) + "  0 0 0 1 1 1 0 0\n";
            }
        }
    }
}

} // namespace redukti::exporters
