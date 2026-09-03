// C++ port of org.redukti.spec.Prescription, VigType and RayOpticsModelBuilder
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
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

enum class VigType {
    None,
    Paraxial,
    SetVig,
    SetPupil,
    SetStopAperture,
    SetApertures,
    SetFnum,
};

class Prescription {
public:
    double _focal_length;
    double _fno;
    double _angle_of_view_in_degrees;
    double _diameter_image_circle;
    std::vector<double> _wvls;
    std::vector<double> _wts;

    std::vector<SurfaceType> _surface_list;
    /**
     * Java copies _surface_list into an array in build(); here the flag records
     * that build() ran, and _surface_list stays the single store so the
     * pointers callers take remain valid.
     */
    bool _built = false;

    /** Null until add_configurations finds a [report data] scenarios row. */
    std::optional<std::vector<int>> _configurations;
    std::optional<std::vector<std::string>> _configuration_names;
    std::vector<double> _angle_of_views_by_scenario;
    std::vector<double> _focal_length_by_scenario;
    std::vector<double> _f_number_by_scenario;

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
    Prescription &asph(int asph_type, double k, const std::vector<double> &coeffs);

    double image_diameter_for_field(double field) const;
    double full_angle_of_view_degrees(double field) const;
    double get_half_angle_in_degrees() const { return _angle_of_view_in_degrees / 2.0; }
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

    int update_apertures_from(rayoptics::optical::OpticalModel *opm, int config,
                              const std::vector<int> *avoid_list,
                              const std::vector<int> *include_list, int decimals);
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
