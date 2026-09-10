// C++ port of org.redukti.tools.LensTool2
//
// Generates the report artifacts that live under Examples/: prescription.txt,
// the .zmx, the layout and spot SVGs, the MTF SVGs and CSVs, vig/paraxial
// dumps and README.md.
#ifndef REDUKTI_TOOLS_LENSTOOL2_H
#define REDUKTI_TOOLS_LENSTOOL2_H

#include "redukti/importers/OpticalBenchDataImporter.h"
#include "redukti/rayoptics/analysis/SpotAnalysis.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/parax/ParaxTypes.h"
#include "redukti/spec/Prescription.h"
#include "redukti/util/Args.h"

#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace redukti::tools {

class LensTool2 {
public:
    using LensSpecifications = importers::OpticalBenchDataImporter::LensSpecifications;

    static LensSpecifications getSpecsFromFile(const std::string &specfile);

    static spec::Prescription createPrescription(const LensSpecifications &specs,
                                                 bool use_glass_types, bool d_line);
    static spec::Prescription createPrescription(const LensSpecifications &specs,
                                                 bool use_glass_types, bool weighted,
                                                 bool d_line);
    static spec::Prescription createPrescription(const LensSpecifications &specs,
                                                 bool use_glass_types,
                                                 const std::vector<double> &wvls,
                                                 const std::vector<double> &wts);

    static std::unique_ptr<rayoptics::optical::OpticalModel> createSystem(
        const spec::Prescription &prescription, bool fov_angle, spec::VigType vig_type,
        bool use_wideangle_aiming, const std::vector<double> &fields, int config);

    static void outputSpotAnalysis(
        const rayoptics::analysis::SpotAnalysisResult::SpotResultsForField &result,
        const std::optional<std::string> &output_file, std::optional<double> radius);

    static std::string &fodToMarkdown(const rayoptics::parax::FirstOrderData &fod,
                                      std::string &sb);
    static std::string &spotResultsMarkdownTable(
        const rayoptics::analysis::SpotAnalysisResult &spotAnalysisResult,
        std::string &sb);

    static std::string startREADME(const LensSpecifications &specs);
    static std::string startREADME(const spec::Prescription &prescription);

    /**
     * Use the final geometry, including optimized airspaces, with the weighted
     * spectrum. Package-private in the Java; public here so the test can reach
     * it.
     */
    static spec::Prescription createWeightedPrescription(
        const spec::Prescription &prescription, bool dLineOnly);
    static std::string &addConfigLabelToREADME(std::string &sb,
                                               const std::optional<std::string> &label);
    static std::string &addLayoutsToREADME(std::string &sb,
                                           const std::string &scenario_filesuffix);
    static std::string &addSpotDiagramsToREADME(std::string &sb,
                                                const std::string &scenario_filesuffix);
    static void addFodToREADME(std::string &sb,
                               const rayoptics::parax::FirstOrderData &fod);
    static void addSpotReportToREADME(
        std::string &sb, const rayoptics::analysis::SpotAnalysisResult &spotAnalysisResult);
    static void addMTFsToREADME(std::string &sb, const std::string &scenario_filesuffix,
                                const std::vector<int> &mtf_freqs);

    /**
     * The Java stamps LocalDate.now() into the README, which would make the
     * output differ every day and defeat comparing against the committed
     * reports. The date is a parameter here; run() passes today's date, and the
     * test passes the date the reference README carries.
     */
    static void createREADME(std::string &sb, const std::string &specFile,
                             const std::string &output_file,
                             const std::string &generated_on);

    static void doLayoutDiagrams(const spec::Prescription &prescription,
                                 const util::Args &arguments, int config,
                                 const std::string &filename_suffix);

    /** The body of the Java `main`, minus argv parsing and the usage banner. */
    static void run(const util::Args &arguments, const std::string &generated_on);

    /** Today's date as Java LocalDate.now() renders it, i.e. ISO yyyy-MM-dd. */
    static std::string today();

private:
    /**
     * Loads the prescription, optionally running the glass type matcher over it
     * first. The matched prescription is used for this run only; it replaces
     * the input file just when --update-specfile asks for that.
     */
    static LensSpecifications loadSpecs(const util::Args &arguments);

    static rayoptics::analysis::SpotAnalysisResult generateSpotDiagrams(
        rayoptics::optical::OpticalModel *opm, const util::Args &arguments,
        bool standardSize, const std::string &filename_suffix);

    static void generateMTFs(rayoptics::optical::OpticalModel *opm,
                             const util::Args &arguments,
                             const std::vector<double> &fields,
                             const std::vector<std::pair<double, double>> &wv_wts,
                             const std::string &outname,
                             const std::string &filename_suffix);

    static void generateRayAberrationPlots(rayoptics::optical::OpticalModel *opm,
                                           const util::Args &arguments,
                                           const std::string &filname_suffix);

    static rayoptics::analysis::SpotOptions spotOptions(const util::Args &arguments);

    static std::unique_ptr<rayoptics::optical::OpticalModel> createLayoutSystem(
        const spec::Prescription &prescription, int config, spec::VigType vigType,
        bool useWideAngleAiming);

    /**
     * Runs the routine airspace optimization. A prime gets its back focus
     * varied; a zoom gets the other variable airspaces varied, one
     * configuration at a time, because a variable cannot yet be shared across
     * configurations.
     */
    static void runDefaultOptimizations(spec::Prescription &prescription,
                                        const util::Args &arguments,
                                        spec::VigType vigType);

    static std::string suffixed_name(const std::string &baseName,
                                     const std::string &suffix, const std::string &ext);
};

} // namespace redukti::tools

#endif // REDUKTI_TOOLS_LENSTOOL2_H
