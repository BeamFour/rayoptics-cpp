// C++ port of org.redukti.spec.SurfaceType
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_SPEC_SURFACETYPE_H
#define REDUKTI_SPEC_SURFACETYPE_H

#include <optional>
#include <string>
#include <vector>

namespace redukti::spec {

/**
 * One prescription surface, plus the per-configuration overrides.
 *
 * `_thickness_by_scenario` and `_diameter_by_scenario` are null in the Java
 * when the surface does not vary across configurations, and the getters fall
 * back to the scalar; that distinction is kept with std::optional because the
 * fallback is what a single-configuration model relies on.
 */
class SurfaceType {
public:
    static constexpr int ASPH_EVEN = 1;
    static constexpr int ASPH_EVEN_A2 = 2;
    static constexpr int ASPH_ODD = 3;

    std::string _id;
    double _radius;
    double _thickness;
    double _diameter;
    bool _is_aperture_stop;
    bool _is_field_stop = false;
    double _nd;
    double _vd;
    /** Nullable in the Java. */
    std::optional<std::string> _catalog_name;
    std::optional<std::string> _glass_name;
    int _asph_type;
    double _k = 0.0;
    std::optional<std::vector<double>> _coeffs;
    std::optional<std::vector<double>> _thickness_by_scenario;
    std::optional<std::vector<double>> _diameter_by_scenario;

    SurfaceType(std::string id, bool is_aperture_stop, double radius, double thickness,
                double diameter, double nd, double vd,
                std::optional<std::string> glass_name,
                std::optional<std::string> catalog_name)
        : _id(std::move(id)), _radius(radius), _thickness(thickness), _diameter(diameter),
          _is_aperture_stop(is_aperture_stop), _nd(nd), _vd(vd),
          _catalog_name(std::move(catalog_name)), _glass_name(std::move(glass_name)),
          _asph_type(0) {}

    SurfaceType &set_thickness_by_scenario(const std::vector<double> &v) {
        _thickness_by_scenario = v;
        return *this;
    }

    SurfaceType &set_diameter_by_scenario(const std::vector<double> &v) {
        _diameter_by_scenario = v;
        return *this;
    }

    bool is_aperture_stop() const { return _is_aperture_stop; }
    bool is_field_stop() const { return _is_field_stop; }

    double get_diameter() const { return _diameter; }
    double get_diameter_by_scenario(int scenario) const;

    double get_thickness() const { return _thickness; }
    double get_thickness_by_scenario(int scenario) const;

    double get_radius_of_curvature() const { return _radius; }

    /** Empty when the surface is not aspheric. */
    std::vector<double> get_aspheric_coeffs() const {
        return _coeffs.has_value() ? *_coeffs : std::vector<double>{};
    }

    double get_cc() const { return _k; }
    double get_refractive_index() const { return _nd; }
    double get_abbe_vd() const { return _vd; }

    const std::optional<std::string> &get_glass_name() const { return _glass_name; }
    const std::optional<std::string> &get_catalog_name() const { return _catalog_name; }

    bool is_aspheric() const { return _asph_type != 0; }
    bool is_odd_asphere() const { return is_aspheric() && _asph_type == ASPH_ODD; }
    bool is_even_a2_asphere() const { return is_aspheric() && _asph_type == ASPH_EVEN_A2; }
    bool is_even_asphere() const { return is_aspheric() && _asph_type == ASPH_EVEN; }

    std::string &to_opt_bench_str(std::string &sb, bool is_last) const;
    std::string &aspherics_to_opt_bench_str(std::string &sb) const;

    static std::string &aspheric_markdown_table_header(std::string &sb, int max_coeffs);
    std::string &asherics_to_markdown_table_row(std::string &sb, int max_coeffs) const;
    std::string &variables_to_markdown_table_row(std::string &sb) const;

    static std::string &to_markdown_table_header(std::string &sb);
    std::string &to_markdown_table_row(std::string &sb) const;

    std::string toString() const;

private:
    std::string asphere_type() const;
};

} // namespace redukti::spec

#endif // REDUKTI_SPEC_SURFACETYPE_H
