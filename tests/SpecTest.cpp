// End-to-end check of the obench importer and the spec package against the JVM.
//
// Reads three real prescriptions from Examples/ -- a zoom with aspheres and
// multiple configurations, a prime with aspheres and named glasses, and a
// purely spherical lens -- prints the parse back out in both output formats,
// builds an optical model under three vignetting modes, and finally pulls the
// computed apertures back into the prescription. Compared line by line with
// output dumped from JDK 25.
//
// This is the first test that starts from a file rather than a hand-built
// model, so it pins the whole path an application actually takes.
#include "SpecExpected.h"
#include "TestHarness.h"

#include "redukti/Text.h"
#include "redukti/importers/OpticalBenchDataImporter.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/parax/ParaxTypes.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/spec/Prescription.h"

#include <string>
#include <vector>

namespace {

using redukti::doubleToString;
using redukti::importers::OpticalBenchDataImporter;
using redukti::spec::Prescription;
using redukti::spec::RayOpticsModelBuilder;
using redukti::spec::VigType;

/** Set by CMake to <repo>/Examples/jfotoptix/. */
const char *const ROOT = REDUKTI_EXAMPLES_DIR;

const char *const FILES[] = {
    "canon-ef11-24mm-f4L/US20150146085_Example01P.txt",
    "canon-ef35mm-f1.4L/JP1999-211978_Example01P.txt",
    "angenieux-180mm-f2.3/US004726669_Example01P.txt",
};

std::string d(double v) { return doubleToString(v); }

/** Java prints a null String reference as "null". */
std::string opt(const std::optional<std::string> &v) {
    return v.has_value() ? *v : std::string("null");
}

std::string b(bool v) { return v ? "true" : "false"; }

const char *surface_type_name(OpticalBenchDataImporter::SurfaceType t) {
    switch (t) {
    case OpticalBenchDataImporter::SurfaceType::surface: return "surface";
    case OpticalBenchDataImporter::SurfaceType::aperture_stop: return "aperture_stop";
    case OpticalBenchDataImporter::SurfaceType::field_stop: return "field_stop";
    }
    return "?";
}

const char *vig_type_name(VigType v) {
    switch (v) {
    case VigType::None: return "None";
    case VigType::Paraxial: return "Paraxial";
    case VigType::SetVig: return "SetVig";
    case VigType::SetPupil: return "SetPupil";
    case VigType::SetStopAperture: return "SetStopAperture";
    case VigType::SetApertures: return "SetApertures";
    case VigType::SetFnum: return "SetFnum";
    }
    return "?";
}

std::vector<std::string> lines(const std::string &text) {
    std::vector<std::string> out;
    std::size_t start = 0;
    while (start < text.size()) {
        auto nl = text.find('\n', start);
        if (nl == std::string::npos) {
            out.push_back(text.substr(start));
            break;
        }
        out.push_back(text.substr(start, nl - start));
        start = nl + 1;
    }
    while (!out.empty() && out.back().empty())
        out.pop_back();
    return out;
}

void report(std::string &sb, const std::string &path) {
    sb += "--- FILE " + path + " ---\n";
    OpticalBenchDataImporter::LensSpecifications specs;
    specs.parse_file(std::string(ROOT) + path);

    const auto &surfaces = specs.get_surfaces();
    sb += "surfaces=" + std::to_string(surfaces.size()) + "\n";
    for (std::size_t i = 0; i < surfaces.size(); i++) {
        const auto &s = surfaces[i];
        sb += "  surf " + std::to_string(i) +
              " type=" + surface_type_name(s.get_surface_type()) +
              " r=" + d(s.get_radius()) + " t0=" + d(s.get_thickness(0)) +
              " dia0=" + d(s.get_diameter(0)) + " nd=" + d(s.get_refractive_index()) +
              " vd=" + d(s.get_abbe_vd()) + " glass=" + opt(s.get_glass_name()) +
              " cat=" + opt(s.get_catalog_name()) +
              " asph=" + b(s.get_aspherical_data() != nullptr) + "\n";
    }
    sb += "image_height=" + d(specs.get_image_height()) + "\n";
    sb += "focal_length=" + d(specs.get_focal_length()) + "\n";
    sb += "aov0=" + d(specs.get_angle_of_view_in_degrees(0)) + "\n";
    sb += "fno0=" + d(specs.get_f_number(0)) + "\n";
    sb += "half_aov_rad0=" + d(specs.get_half_angle_of_view_in_radians(0)) + "\n";

    Prescription p = Prescription::build_prescription(specs, true, false, false);
    sb += "num_configurations=" + std::to_string(p.get_num_configurations()) + "\n";
    sb += "title=" + p.get_title() + "\n";
    sb += "--- OPT BENCH ---\n";
    p.to_opt_bench_str(sb);
    sb += "--- MARKDOWN ---\n";
    p.to_markdown_str(sb);

    const VigType vts[] = {VigType::None, VigType::Paraxial, VigType::SetVig};
    for (VigType vt : vts) {
        auto model = RayOpticsModelBuilder(p).build_optical_model(
            true, std::vector<double>{0.0, 0.707, 1.0}, false, vt, false, 0);
        const auto &fod = model->optical_spec->parax_data->fod;
        sb += std::string("--- MODEL ") + vig_type_name(vt) + " ---\n";
        sb += "efl=" + d(fod.efl) + "\n";
        sb += "bfl=" + d(fod.bfl) + "\n";
        sb += "ffl=" + d(fod.ffl) + "\n";
        sb += "fno=" + d(fod.fno) + "\n";
        sb += "img_ht=" + d(fod.img_ht) + "\n";
        sb += "enp_dist=" + d(fod.enp_dist) + "\n";
        sb += "enp_radius=" + d(fod.enp_radius) + "\n";
        sb += "exp_dist=" + d(fod.exp_dist) + "\n";
        sb += "exp_radius=" + d(fod.exp_radius) + "\n";
        sb += "opt_inv=" + d(fod.opt_inv) + "\n";
        sb += "obj_ang=" + d(fod.obj_ang) + "\n";
        sb += "num_ifcs=" + std::to_string(model->seq_model->ifcs.size()) + "\n";
        const auto &fields = model->optical_spec->fov->fields;
        for (std::size_t fi = 0; fi < fields.size(); fi++) {
            const auto &f = *fields[fi];
            sb += "  fld " + std::to_string(fi) + " y=" + d(f.y) + " vlx=" + d(f.vlx) +
                  " vly=" + d(f.vly) + " vux=" + d(f.vux) + " vuy=" + d(f.vuy) + "\n";
        }
    }
    auto model = RayOpticsModelBuilder(p).build_optical_model(
        true, std::vector<double>{0.0, 0.707, 1.0}, true, VigType::SetVig, false, 0);
    int changed = p.update_apertures_from(model.get(), 0);
    sb += "apertures_changed=" + std::to_string(changed) + "\n";
    for (const auto &s : p.get_surfaces())
        sb += "  dia " + s._id + "=" + d(s.get_diameter()) + "\n";
}

} // namespace

TEST(spec_importer_and_model_builder_match_jvm) {
    std::string sb;
    for (const char *f : FILES)
        report(sb, f);

    auto actual = lines(sb);
    const std::size_t expected_n = sizeof(EXPECTED_SPEC_LINES) / sizeof(EXPECTED_SPEC_LINES[0]);
    CHECK_EQ(actual.size(), expected_n);
    for (std::size_t i = 0; i < actual.size() && i < expected_n; i++)
        CHECK_STR_EQ(actual[i], std::string(EXPECTED_SPEC_LINES[i]));
}
