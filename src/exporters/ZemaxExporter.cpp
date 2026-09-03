// C++ port of org.redukti.exporters.ZemaxExporter
#include "redukti/exporters/ZemaxExporter.h"

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
    sb += "FTYP 3 0 9 3 0 0 0 9\n";
    sb += "ROPD 2\n";
    sb += "HYPR 0\n";
    sb += "PICB 1\n";
    sb += "XFLN 0 0 0 0 0 0 0 0 0 0 0 0\n";
    sb += "YFLN 0 2.16 4.33 6.5 8.65 10.8 12.98 15.14 " + formatF(img_ht, 6) +
          " 0 0 0 0 0 0 0 0\n";
    sb += "FWGN 10 9 9 8 8 7 7 7 3 1 1 1\n";
    sb += "VDXN 0 0 0 0 0 0 0 0 0 0 0 0\n";
    sb += "VDYN 0 0 0 0 0 0 0 0 0 0 0 0\n";
    sb += "VCXN 0 0 0 0 0 0 0 0 0 0 0 0\n";
    sb += "VCYN 0 0 0 0 0 0 0 0 0 0 0 0\n";
    sb += "VANN 0 0 0 0 0 0 0 0 0 0 0 0\n";
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
            if (s.is_odd_asphere())
                sb += "  TYPE ODDASPHE\n";
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
            for (std::size_t a = 1; a <= aspherics.size(); a++) {
                sb += "  PARM " + i(static_cast<int>(a)) + " ";
                sb += d(aspherics[a - 1]) + "\n";
            }
        }
        sb += "  DISZ " + d(thickness) + "\n";
        if (s.is_aspheric()) {
            double k_conic = s.is_odd_asphere() ? s.get_cc() + 1 : s.get_cc();
            sb += "  CONI " + d(k_conic) + "\n";
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
}

} // namespace redukti::exporters
