// C++ port of org.redukti.spec.Prescription, VigType and RayOpticsModelBuilder
#ifndef REDUKTI_SPEC_PRESCRIPTION_H
#define REDUKTI_SPEC_PRESCRIPTION_H

#include "redukti/importers/OpticalBenchDataImporter.h"
#include "redukti/spec/SurfaceType.h"

#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace redukti::rayoptics::optical {
class OpticalModel;
}

namespace redukti::rayoptics::seq {
class SequentialModel;
}

namespace redukti::spec {

/**
 * Which aperture and vignetting calculation to run once a model is built.
 *
 * Note the direction each one runs in. #SetPupil derives the pupil spec
 * from the authored stop diameter, so it changes the f/# and leaves apertures
 * alone. #SetStopAperture and #SetFnum go the other way: they
 * hold the f/# and size the stop to satisfy it.
 */
enum class VigType {
    None,
    Paraxial,
    /** Vignetting factors from the existing apertures. */
    SetVig,
    /** Pupil spec derived from the existing stop diameter. Apertures unchanged. */
    SetPupil,
    /** Stop sized to satisfy the pupil spec, then vignetting recalculated. */
    SetStopAperture,
    /** Vignetting, then every clear aperture sized to pass the vignetted rays. */
    SetApertures,
    /**
     * Drive the whole aperture set from the defined f/#: size the stop to
     * satisfy it, then size every other clear aperture to pass the resulting
     * vignetted rays.
     *
     * For prescriptions that quote an exact f/# but whose apertures were
     * estimated off a drawing - the stop is trusted to the spec and everything
     * else is rebuilt from the rays.
     */
    SetFnum,
};

/**
 * A format for prescription that is easier to work with when trying to optimize.
 * Supports both single and multiple configurations - a multi configuration setup is
 * required for zoom lenses where thicknesses, focal lengths, diameters, fnumbers can
 * vary by zoom setting.
 * There is always a default scenario configuration.
 * During optimization the default scenario is used so for multi-config prescriptions
 * the default scenario must be set appropriately (to be tested).
 */
class Prescription {
public:
    /** Focal length of default scenario, in multi config this is defined by _focal_length_by_scenario */
    double _focal_length;
    /** F-number of default scenario, in multi config this is defined by _f_number_by_scenario */
    double _fno;
    /** The quoted angle of view - e.g. 47 degrees for 50mm; this is the full angle of view
     * This is defined by _angle_of_views_by_scenario in multi config
     */
    double _angle_of_view_in_degrees;
    /** For 35mm this is sqrt(36^2 + 24^2) = 43.27 */
    double _diameter_image_circle;
    /** wavelengths to use,
     * NOTE atm first wvl is made reference wvl
     */
    std::vector<double> _wvls;
    /** wavelength weights - mainly used when computing spot and MTFs */
    std::vector<double> _wts;

    /** when building we use a list */
    std::vector<SurfaceType> _surface_list;
    /** After construction this is the list of surfaces */
    /**
     * Java copies _surface_list into an array in build(); here the flag records
     * that build() ran, and _surface_list stays the single store so the
     * pointers callers take remain valid.
     */
    bool _built = false;

    /** Maps our config id to scenario number in the OpticalBench specs
     * The scenario in OBench corresponds to how it is defined in
     * source patent data.
     * Note that these scenario numbers get lost when we write new
     * prescription because the new prescription ends up only with the
     * configured scenarios.
     */
    /** Null until add_configurations finds a [report data] scenarios row. */
    std::optional<std::vector<int>> _configurations;
    /** Each config is given a name */
    std::optional<std::vector<std::string>> _configuration_names;
    /** Each config has its own angle of view */
    std::vector<double> _angle_of_views_by_scenario;
    /** Each config has its own focal length */
    std::vector<double> _focal_length_by_scenario;
    /** Each config has its own fnumber */
    std::vector<double> _f_number_by_scenario;
    // There are other values such as thickness and aperture
    // that vary for configurations but these are specified for
    // each surface

    /** Following are optional values for information only, used to generate
     * lens report.
     */
    std::string _title;
    std::optional<std::string> _lens_name;
    std::string _patent_country;
    std::optional<std::string> _patent_number;
    std::string _patent_example;
    std::string _application_year;
    std::string _inventors;
    std::string _original_assignee;
    std::string _current_assignee;
    std::string _patent_link;

    static constexpr double WT_d = 1.0;
    static constexpr double WT_C = 0.475;
    static constexpr double WT_e = 0.98;
    static constexpr double WT_F = 0.49;
    static constexpr double WT_g = 0.15;

    static constexpr int APERTURE_DECIMALS = 4;

    Prescription(double focal_length, double fno, double angle_of_view_degrees,
                 double diameter_image_circle, bool d_line_only);

    Prescription(double focal_length, double fno, double angle_of_view_degrees,
                 double diameter_image_circle, std::vector<double> wvls,
                 std::vector<double> wts);

    Prescription &surf(double radius, double thickness, double diameter, double nd,
                       double vd, const std::optional<std::string> &glass_name,
                       const std::optional<std::string> &catalog_name);
    Prescription &surf(double radius, double thickness, double diameter, double nd,
                       double vd);
    Prescription &surf(double radius, double thickness, double diameter);
    Prescription &stop(double thickness, double diameter);
    Prescription &field_stop(double thickness, double diameter);
    /** Sets up the last added surface as an asphere.
     * @param asph_type 1=EVEN,2=EVEN_A2,3=ODD
     * @param coeffs for EVEN first param is not used and should be 0, for ODD first 2 must be 0
     */
    Prescription &asph(int asph_type, double k, const std::vector<double> &coeffs);

    /**
     * Derive diameter for given field.
     * @param field fields are relative, 0 for axis to 1 for edge.
     */
    double image_diameter_for_field(double field) const;
    /**
     * Full angle of view in degrees for given field for default configuration
     * @param field fields are relative, 0 for axis to 1 for edge.
     */
    double full_angle_of_view_degrees(double field) const;
    /** Get angle of view for default configuration */
    double get_half_angle_in_degrees() const { return _angle_of_view_in_degrees / 2.0; }
    /** Get angle of view for default configuration */
    double get_half_angle_of_view_in_radians() const;

    Prescription &build();

    using LensSpecifications = importers::OpticalBenchDataImporter::LensSpecifications;

    static Prescription build_prescription_d_line(const LensSpecifications &specs);
    static Prescription build_prescription_e_line(const LensSpecifications &specs);
    static Prescription build_prescription(const LensSpecifications &specs,
                                           bool use_glass_types);
    static Prescription build_prescription(const LensSpecifications &specs,
                                           bool use_glass_types,
                                           const std::vector<double> &wvls,
                                           const std::vector<double> &wts);
    static Prescription build_prescription(const LensSpecifications &specs,
                                           bool use_glass_types, bool weighted,
                                           bool d_line);
    static Prescription build_prescription(const LensSpecifications &specs,
                                           bool use_glass_types,
                                           const std::vector<double> &wvls,
                                           const std::vector<double> &wts,
                                           int default_scenario);

    const std::string &get_title() const { return _title; }

    /**
     * Copy computed apertures from a traced model back into this prescription.
     *
     * Building a model never writes back, so a prescription keeps whatever
     * apertures it was authored with until this is called deliberately. Use it
     * after the model has been given real apertures, by
     * org.redukti.rayoptics.raytr.VigCalc#set_ape for every surface or
     * org.redukti.rayoptics.raytr.VigCalc#set_stop_aperture for the stop
     * alone.
     *
     * Diameters are taken from Interface.surface_od(), which is what the
     * layout and the Zemax exporter already use, and rounded so that rewriting a
     * prescription does not churn every diameter into full double precision.
     *
     * Surfaces here are indexed as in #_surfaces. The model carries an
     * object and an image interface either side of them, so prescription surface
     * i is model interface i + 1 - the lists taken by
     * org.redukti.rayoptics.raytr.VigCalc#set_clear_apertures are model
     * indices and are offset by one from these.
     *
     * @param opm          a model built from this prescription
     * @param config       configuration index the model was built for
     * @param avoid_list   prescription surfaces to leave alone, or null
     * @param include_list prescription surfaces to update, or null for all.
     *                     Ignored when avoid_list is given.
     * @param decimals     decimal places to round to
     * @return the number of surfaces whose diameter changed
     * @throws IllegalArgumentException if the model was not built from this
     *                                  prescription, or config is out of range
     */
    int update_apertures_from(rayoptics::optical::OpticalModel *opm, int config,
                              const std::vector<int> *avoid_list,
                              const std::vector<int> *include_list, int decimals);
    /** Update every surface for the given configuration. */
    int update_apertures_from(rayoptics::optical::OpticalModel *opm, int config) {
        return update_apertures_from(opm, config, nullptr, nullptr, APERTURE_DECIMALS);
    }

    int get_num_configurations() const {
        return _configurations.has_value() ? static_cast<int>(_configurations->size()) : 0;
    }

    double get_f_number() const { return _fno; }

    const std::vector<SurfaceType> &get_surfaces() const { return _surface_list; }
    std::vector<SurfaceType> &get_surfaces() { return _surface_list; }

    bool has_odd_aspheric() const;
    bool has_even_a2_aspheric() const;

    std::string &to_opt_bench_str(std::string &sb) const;
    std::string &to_markdown_str(std::string &sb) const;

    std::string organization() const;

    /** Java returns a LinkedHashMap; insertion order is the wavelength order. */
    std::vector<std::pair<double, double>> get_wvl_wts() const;

    std::string toString() const;

private:
    Prescription &import_surface(
        const importers::OpticalBenchDataImporter::LensSurface &surface, int scenario,
        bool use_glass_types);
    void add_configuration_data(
        const importers::OpticalBenchDataImporter::LensSurface &lensSurface);
    Prescription &add_configurations(const LensSpecifications &specs);

    //sb.append("angle of view = ").append(fullAngleOfViewDegrees(1.0)).append("\n");
    void add_patent_section(std::string &sb) const;
    void add_report_section(std::string &sb) const;
};

class RayOpticsModelBuilder {
public:
    Prescription _prescription;

    explicit RayOpticsModelBuilder(Prescription prescription)
        : _prescription(std::move(prescription)) {}

    std::unique_ptr<rayoptics::optical::OpticalModel> build_optical_model(
        bool fov_angle, const std::vector<double> &fields, bool do_apertures,
        VigType vig_type, bool use_wideangle_aiming, int config);

private:
    void add_surface(rayoptics::seq::SequentialModel *sm, const SurfaceType &s,
                     int config);
};

} // namespace redukti::spec

#endif // REDUKTI_SPEC_PRESCRIPTION_H
