// C++ port of org.redukti.tools.LensTool2
#include "redukti/tools/LensTool2.h"

#include "redukti/Text.h"
#include "redukti/exporters/ZemaxExporter.h"
#include "redukti/mathlib/M.h"
#include "redukti/plotter/Plotter.h"
#include "redukti/rayoptics/analysis/MTF.h"
#include "redukti/rayoptics/analysis/RayAberrations.h"
#include "redukti/rayoptics/layout/Layout2D.h"
#include "redukti/rayoptics/raytr/RayTrace.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <iostream>

namespace redukti::tools {

using importers::OpticalBenchDataImporter;
using plotter::GeoMTFByFieldPlot;
using plotter::GeoMTFPlot;
using plotter::RayAberrationPlot;
using plotter::SpotDiagram;
using spec::Prescription;
using spec::RayOpticsModelBuilder;
using spec::VigType;
using util::Args;
using util::Helper;
using namespace redukti::rayoptics;

namespace {

/** Java's static DecimalFormat decimalFormat = M.decimal_format(). */
const DecimalFormat &decimalFormat() {
    static const DecimalFormat instance = mathlib::M::decimal_format();
    return instance;
}

std::string d(double v) { return doubleToString(v); }

} // namespace

LensTool2::LensSpecifications LensTool2::getSpecsFromFile(const std::string &specfile) {
    LensSpecifications specs;
    specs.parse_file(specfile);
    return specs;
}

Prescription LensTool2::createPrescription(const LensSpecifications &specs,
                                           bool use_glass_types, bool d_line) {
    return Prescription::build_prescription(specs, use_glass_types, false, d_line);
}

Prescription LensTool2::createPrescription(const LensSpecifications &specs,
                                           bool use_glass_types, bool weighted,
                                           bool d_line) {
    return Prescription::build_prescription(specs, use_glass_types, weighted, d_line);
}

std::unique_ptr<optical::OpticalModel> LensTool2::createSystem(
    const Prescription &prescription, bool fov_angle, VigType vig_type,
    bool use_wideangle_aiming, const std::vector<double> &fields, int config) {
    return RayOpticsModelBuilder(prescription)
        .build_optical_model(fov_angle, fields, false, vig_type, use_wideangle_aiming,
                             config);
}

void LensTool2::outputSpotAnalysis(
    const analysis::SpotAnalysisResult::SpotResultsForField &result,
    const std::optional<std::string> &output_file, std::optional<double> radius) {
    if (output_file.has_value())
        Helper::createOutputFile(*output_file, SpotDiagram(result).plot(radius));
    else
        std::cout << SpotDiagram(result).plot(radius) << std::endl;
}

std::string &LensTool2::fodToMarkdown(const parax::FirstOrderData &fod, std::string &sb) {
    const DecimalFormat &f = decimalFormat();
    sb += "| parameter | value |\n";
    sb += "| ---       | ---   |\n";
    // obj_dist and red are appended raw in the Java, not through the format.
    sb += "| effective_focal_length |" + f.format(fod.efl) +
          "\n| back_focal_length | " + f.format(fod.bfl) +
          "\n| optical_invariant | " + f.format(fod.opt_inv) +
          "\n| object_distance | " + d(fod.obj_dist) +
          "\n| image_distance | " + f.format(fod.img_dist) +
          "\n| power | " + f.format(fod.power) +
          "\n| pp1_H | " + f.format(fod.pp1) +
          "\n| ppk_H' | " + f.format(fod.ppk) +
          "\n| ffl_F | " + f.format(fod.ffl) +
          "\n| fno | " + f.format(fod.fno) +
          "\n| enp_dist_P | " + f.format(fod.enp_dist) +
          "\n| enp_radius | " + f.format(fod.enp_radius) +
          "\n| exp_dist_P' | " + f.format(fod.exp_dist) +
          "\n| exp_radius | " + f.format(fod.exp_radius) +
          "\n| m | " + f.format(fod.m) +
          "\n| red | " + d(fod.red) +
          "\n| n_obj | " + f.format(fod.n_obj) +
          "\n| n_img | " + f.format(fod.n_img) +
          "\n| img_ht | " + f.format(fod.img_ht) +
          "\n| obj_ang | " + f.format(fod.obj_ang) +
          "\n| obj_na | " + f.format(fod.obj_na) +
          "\n| img_na | " + f.format(fod.img_na) + "|\n";
    return sb;
}

std::string &LensTool2::spotResultsMarkdownTable(
    const analysis::SpotAnalysisResult &spotAnalysisResult, std::string &sb) {
    const DecimalFormat &f = decimalFormat();
    sb += "| Field | Spot Mean Radius | Spot Max Radius |\n";
    sb += "| ---   | ---              | ---             |\n";
    for (const auto &result : spotAnalysisResult.spot_results) {
        sb += " | " + result.fld->toString();
        sb += " | " + f.format(result.get_mean_radius());
        sb += " | " + f.format(result.get_max_radius());
        sb += "|\n";
    }
    return sb;
}

std::string LensTool2::startREADME(const LensSpecifications &specs) {
    Prescription prescription = Prescription::build_prescription(specs, true);
    std::string sb;
    prescription.to_markdown_str(sb);
    return sb;
}

std::string &LensTool2::addConfigLabelToREADME(std::string &sb,
                                               const std::optional<std::string> &label) {
    if (label.has_value())
        sb += "# " + *label + "\n";
    return sb;
}

std::string &LensTool2::addLayoutsToREADME(std::string &sb,
                                           const std::string &scenario_filesuffix) {
    sb += "## Layouts\n";
    sb += "![Layout Elements](./layoutonly" + scenario_filesuffix + ".svg)\n";
    sb += "![Layout](./layout" + scenario_filesuffix + ".svg)\n";
    return sb;
}

std::string &LensTool2::addSpotDiagramsToREADME(std::string &sb,
                                                const std::string &scenario_filesuffix) {
    sb += "## Spot Diagrams\n";
    sb += "![Spot Diagram Field 0.0](./spot" + scenario_filesuffix + ".svg)\n";
    sb += "![Spot Diagram Field 0.7](./spot-semi-skew" + scenario_filesuffix + ".svg)\n";
    sb += "![Spot Diagram Field 1.0](./spot-skew" + scenario_filesuffix + ".svg)\n";
    return sb;
}

void LensTool2::addFodToREADME(std::string &sb, const parax::FirstOrderData &fod) {
    sb += "## Paraxial Parameters\n";
    fodToMarkdown(fod, sb);
}

void LensTool2::addSpotReportToREADME(
    std::string &sb, const analysis::SpotAnalysisResult &spotAnalysisResult) {
    sb += "## Spot Analysis\n";
    spotResultsMarkdownTable(spotAnalysisResult, sb);
}

void LensTool2::addMTFsToREADME(std::string &sb, const std::string &scenario_filesuffix,
                                const std::vector<int> &mtf_freqs) {
    std::string freq_legend = GeoMTFByFieldPlot::freq_legend(mtf_freqs);
    sb += "## Polychromatic Geometric MTF\n";
    sb += "![Polychromatic Geometrical MTF](./mtf" + scenario_filesuffix + ".svg)\n";
    sb += "* " + freq_legend + " cycles/mm\n";
    sb += "* Solid lines represent sagittal, dashed lines tangential\n";
    sb += "* To generate above, MTFs for wavelengths 587.5618(d), 486.1327(F), "
          "656.2725(C) were calculated across 10 fields, and then averaged\n";
    sb += "## Polychromatic Geometric MTF (Weighted)\n";
    sb += "![Polychromatic Geometrical MTF Weighted](./mtf-w" + scenario_filesuffix +
          ".svg)\n";
    sb += "* " + freq_legend + " cycles/mm\n";
    sb += "* Solid lines represent sagittal, dashed lines tangential\n";
    sb += "* To generate above, MTFs for wavelengths 587.5618(d) wt(1.0), 656.2725(C) "
          "wt(0.475), 546.074(e) wt(0.98), 486.1327(F) wt(0.49), 435.8343(g) wt(0.15) "
          "were calculated across 10 fields, and then combined using weighted average\n";
}

void LensTool2::createREADME(std::string &sb, const std::string &specFile,
                             const std::string &output_file,
                             const std::string &generated_on) {
    std::string filename = Helper::getFilename(specFile);
    std::string zmxFilename = Helper::replaceExtension(filename, ".zmx");
    sb += "## Resources\n";
    sb += "* [OpticalBench Compatible Data File, tab delimited](./prescription.txt)\n";
    sb += "* [Zemax file](./" + zmxFilename + ")\n\n";
    sb += "Report / Zemax file generated using "
          "[Beam42](https://github.com/BeamFour/Beam42) on " +
          generated_on + "\n";
    Helper::createOutputFile(output_file, sb);
}

std::string LensTool2::suffixed_name(const std::string &baseName,
                                     const std::string &suffix, const std::string &ext) {
    return baseName + suffix + ext;
}

analysis::SpotAnalysisResult LensTool2::generateSpotDiagrams(
    optical::OpticalModel *opm, const Args &arguments, bool standardSize,
    const std::string &filename_suffix) {
    auto spotAnalysis = analysis::SpotAnalysis::eval(opm, analysis::SpotOptions());
    Helper::createOutputFile(
        Helper::getOutputFileWithPath(*arguments.specfile,
                                      suffixed_name("spot-report", filename_suffix,
                                                    ".txt"),
                                      arguments.outdir),
        spotAnalysis.toString());
    for (std::size_t i = 0; i < spotAnalysis.spot_results.size(); i++) {
        const auto &spotFld = spotAnalysis.spot_results[i];
        std::optional<std::string> filename;
        if (spotFld.fld->y == 0.0)
            filename = suffixed_name("spot", filename_suffix, ".svg");
        else if (spotFld.fld->y == 0.7)
            filename = suffixed_name("spot-semi-skew", filename_suffix, ".svg");
        else if (spotFld.fld->y == 1.0)
            filename = suffixed_name("spot-skew", filename_suffix, ".svg");
        if (!filename.has_value())
            continue;
        auto outfile = Helper::getOutputFileWithPath(*arguments.specfile, *filename,
                                                     arguments.outdir);
        outputSpotAnalysis(spotFld, outfile,
                           standardSize ? std::optional<double>(600.)
                                        : std::optional<double>());
    }
    return spotAnalysis;
}

void LensTool2::generateMTFs(optical::OpticalModel *opm, const Args &arguments,
                             const std::vector<double> &fields,
                             const std::vector<std::pair<double, double>> &wv_wts,
                             const std::string &outname,
                             const std::string &filename_suffix) {
    auto spotAnalysis = analysis::SpotAnalysis::eval(opm, analysis::SpotOptions());
    std::vector<analysis::PolyMTF> mtfs;
    for (std::size_t i = 0; i < spotAnalysis.spot_results.size(); i++) {
        const auto &spotFld = spotAnalysis.spot_results[i];
        auto cfg = spotFld.mtfHistogramConfig();
        std::optional<analysis::PolyMTF> polyMtfForField;
        for (const auto &intercepts : spotFld.intercepts) {
            std::string filename =
                suffixed_name("mtf-fld" + std::to_string(i) + "-" +
                                  std::to_string(static_cast<int>(intercepts.wvl)),
                              filename_suffix, ".svg");
            auto output_file = Helper::getOutputFileWithPath(*arguments.specfile, filename,
                                                             arguments.outdir);
            analysis::MonochromaticGeometricMTF mtf(intercepts, cfg);
            if (!polyMtfForField.has_value())
                polyMtfForField.emplace(mtf.mtf.fft_size, mtf.h2d.pixel_size);
            // Java: wv_wts.getOrDefault(intercepts.wvl, 0.0) -- exact key match.
            double wt = 0.0;
            for (const auto &kv : wv_wts) {
                if (kv.first == intercepts.wvl) {
                    wt = kv.second;
                    break;
                }
            }
            if (wt != 0.0)
                polyMtfForField->add(mtf.mtf, wt);
            if (arguments.do_mono_chrome_mtfs)
                Helper::createOutputFile(output_file,
                                         GeoMTFPlot(*spotFld.fld, mtf).plot());
        }
        if (polyMtfForField.has_value()) {
            polyMtfForField->compute();
            mtfs.push_back(std::move(*polyMtfForField));
        }
    }
    const std::vector<int> &freqs = arguments.mtf_freqs;
    std::vector<analysis::MTFResultByFreq> mtfResults;
    for (int freq : freqs)
        mtfResults.push_back(analysis::MTFResultByFreq(mtfs, freq));
    auto mtffile = Helper::getOutputFileWithPath(
        *arguments.specfile, suffixed_name(outname, filename_suffix, ".svg"),
        arguments.outdir);
    GeoMTFByFieldPlot plot(mtfResults, fields);
    Helper::createOutputFile(mtffile, plot.plot());
    auto mtfdata = Helper::getOutputFileWithPath(
        *arguments.specfile, suffixed_name(outname, filename_suffix, ".csv"),
        arguments.outdir);
    Helper::createOutputFile(mtfdata, plot.toString());
}

void LensTool2::generateRayAberrationPlots(optical::OpticalModel *opm,
                                           const Args &arguments,
                                           const std::string &filname_suffix) {
    raytr::TraceOptions traceOptions;
    auto rayAber = analysis::TransverseRayAberrationAnalysis::eval(opm, 21, false,
                                                                   traceOptions);
    for (const auto &fan_results : rayAber.results) {
        std::string filename =
            suffixed_name("rayabbr-fld" + std::to_string(fan_results.fi) + "-" +
                              (fan_results.xy == 1 ? "tan" : "sag"),
                          filname_suffix, ".svg");
        auto output_file = Helper::getOutputFileWithPath(*arguments.specfile, filename,
                                                         arguments.outdir);
        Helper::createOutputFile(output_file,
                                 RayAberrationPlot(rayAber).plot(fan_results, 0));
    }
    raytr::TraceOptions traceOptions2;
    auto opdAber = analysis::WavefrontAberrationAnalysis::eval(opm, 21, false,
                                                                traceOptions2);
    for (const auto &fan_results : opdAber.results) {
        std::string filename =
            suffixed_name("opdabbr-fld" + std::to_string(fan_results.fi) + "-" +
                              (fan_results.xy == 1 ? "tan" : "sag"),
                          filname_suffix, ".svg");
        auto output_file = Helper::getOutputFileWithPath(*arguments.specfile, filename,
                                                         arguments.outdir);
        Helper::createOutputFile(output_file,
                                 RayAberrationPlot(opdAber).plot(fan_results, 0));
    }
}

std::unique_ptr<optical::OpticalModel> LensTool2::createLayoutSystem(
    const Prescription &prescription, int config, VigType vigType,
    bool useWideAngleAiming) {
    std::vector<double> layoutFields{0.0, 1.0};
    return createSystem(prescription, true, vigType, useWideAngleAiming, layoutFields,
                        config);
}

void LensTool2::doLayoutDiagrams(const Prescription &prescription, const Args &arguments,
                                 int config, const std::string &filename_suffix) {
    // For very wide angle lenses, blindly spraying rays doesn't work very well,
    // so the layout system is aimed through the pupil.
    auto opm = createLayoutSystem(prescription, config, VigType::SetPupil, true);
    layout::Layout2D lay;
    auto output = Helper::getOutputFileWithPath(
        *arguments.specfile, suffixed_name("layout-fan", filename_suffix, ".svg"),
        arguments.outdir);
    layout::LayoutOptions fanOpts;
    fanOpts.drawReferenceRays_(false).fanRayCount_(9).clipRays_(true).useTraceFan_(true);
    Helper::createOutputFile(output, lay.renderSvg(opm.get(), 1000, 500, &fanOpts));

    output = Helper::getOutputFileWithPath(
        *arguments.specfile, suffixed_name("layoutonly", filename_suffix, ".svg"),
        arguments.outdir);
    layout::LayoutOptions elementsOpts;
    elementsOpts.drawReferenceRays_(false);
    Helper::createOutputFile(output, lay.renderSvg(opm.get(), 1000, 500, &elementsOpts));

    layout::LayoutOptions defaultOpts;
    std::string reference = lay.renderSvg(opm.get(), 1000, 500, &defaultOpts);
    output = Helper::getOutputFileWithPath(
        *arguments.specfile, suffixed_name("layout", filename_suffix, ".svg"),
        arguments.outdir);
    Helper::createOutputFile(output, reference);
}

std::string LensTool2::today() {
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::tm tm{};
#ifdef _WIN32
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif
    char buf[16];
    std::snprintf(buf, sizeof(buf), "%04d-%02d-%02d", tm.tm_year + 1900, tm.tm_mon + 1,
                  tm.tm_mday);
    return buf;
}

void LensTool2::run(const Args &arguments, const std::string &generated_on) {
    const std::vector<double> fields{0.0, 0.1, 0.2, 0.3, 0.4, 0.5,
                                     0.6, 0.7, 0.8, 0.9, 1.0};
    VigType vigType = VigType::SetPupil;
    LensSpecifications specs = getSpecsFromFile(*arguments.specfile);
    auto prescription =
        createPrescription(specs, arguments.use_glass_types, arguments.only_d_line);
    std::string prescription_output;
    prescription.to_opt_bench_str(prescription_output);
    Helper::createOutputFile(Helper::getOutputFileWithPath(*arguments.specfile,
                                                           "prescription.txt",
                                                           arguments.outdir),
                             prescription_output);
    exporters::ZemaxExporter zemaxExporter;
    Helper::createOutputFile(
        Helper::getOutputPathChangeExt(*arguments.specfile, ".zmx"),
        zemaxExporter.generate(prescription, arguments.only_d_line));
    std::string SB = startREADME(specs);
    const int configs = std::max(prescription.get_num_configurations(), 1);
    for (int config = 0; config < configs; config++) {
        if (prescription.get_num_configurations() > 0)
            addConfigLabelToREADME(
                SB, (*prescription._configuration_names)[static_cast<std::size_t>(config)]);
        std::string scenario_filesuffix =
            prescription.get_num_configurations() > 0 ? ("-" + std::to_string(config))
                                                      : "";
        auto opm = createSystem(prescription, true, vigType, true, fields, config);
        auto sm = opm->seq_model.get();
        auto osp = opm->optical_spec.get();
        const auto &fod = opm->optical_spec->parax_data->fod;
        std::string surfaces;
        sm->list_surfaces(surfaces);
        std::cout << surfaces << std::endl;
        std::string gaps;
        sm->list_gaps(gaps);
        std::cout << gaps << std::endl;
        std::string specText;
        osp->list_str(specText);
        std::cout << specText << std::endl;
        Helper::createOutputFile(
            Helper::getOutputFileWithPath(*arguments.specfile,
                                          suffixed_name("vig", scenario_filesuffix,
                                                        ".txt"),
                                          arguments.outdir),
            specText);
        Helper::createOutputFile(
            Helper::getOutputFileWithPath(*arguments.specfile,
                                          suffixed_name("paraxial", scenario_filesuffix,
                                                        ".txt"),
                                          arguments.outdir),
            fod.toString());
        doLayoutDiagrams(prescription, arguments, config, scenario_filesuffix);
        auto spotAnalysis = generateSpotDiagrams(opm.get(), arguments,
                                                 !arguments.auto_size_spots,
                                                 scenario_filesuffix);
        addLayoutsToREADME(SB, scenario_filesuffix);
        addSpotDiagramsToREADME(SB, scenario_filesuffix);
        addFodToREADME(SB, fod);
        addSpotReportToREADME(SB, spotAnalysis);
        addMTFsToREADME(SB, scenario_filesuffix, arguments.mtf_freqs);
        generateMTFs(opm.get(), arguments, fields, prescription.get_wvl_wts(), "mtf",
                     scenario_filesuffix);
        if (arguments.do_ray_aberrations)
            generateRayAberrationPlots(opm.get(), arguments, scenario_filesuffix);
        auto prescriptionForWeightedMTF =
            createPrescription(specs, arguments.use_glass_types, true,
                               arguments.only_d_line);
        auto opm2 = createSystem(prescriptionForWeightedMTF, true, vigType, true, fields,
                                 config);
        generateMTFs(opm2.get(), arguments, fields,
                     prescriptionForWeightedMTF.get_wvl_wts(), "mtf-w",
                     scenario_filesuffix);
    }
    createREADME(SB, *arguments.specfile,
                 Helper::getOutputFileWithPath(*arguments.specfile, "README.md",
                                               arguments.outdir),
                 generated_on);
}

} // namespace redukti::tools
