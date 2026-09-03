// C++ port of org.redukti.importers.obench.OpticalBenchDataImporter
//
// Reads the tab-separated OpticalBench prescription format used by the files
// under Examples/.
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_IMPORTERS_OPTICALBENCHDATAIMPORTER_H
#define REDUKTI_IMPORTERS_OPTICALBENCHDATAIMPORTER_H

#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace redukti::importers {

/**
 * Java's `Double.parseDouble` wrapped to answer 0.0 on anything it rejects.
 *
 * Not strtod: Java accepts "Infinity" and "NaN" and a trailing d/D/f/F, and
 * rejects trailing garbage that strtod would silently stop at. An empty or
 * absent field is 0.0, which the format relies on.
 */
double parse_double(const std::string &s);

/** Java's `Integer.parseInt` wrapped to answer `defaultValue` on rejection. */
int parse_integer(const std::string &s, int defaultValue);

class OpticalBenchDataImporter {
public:
    /** One named row, which may carry a value per scenario. */
    class Variable {
    public:
        explicit Variable(std::string name) : _name(std::move(name)) {}

        const std::string &name() const { return _name; }

        void add_value(const std::string &value) { _values.push_back(value); }

        int num_scenarios() const { return static_cast<int>(_values.size()); }
        int num_values() const { return static_cast<int>(_values.size()); }

        const std::string &get_value(int scenario) const {
            return _values[static_cast<std::size_t>(scenario)];
        }

        double get_value_as_double(int scenario) const {
            return parse_double(get_value(scenario));
        }

        int get_value_as_integer(int scenario, int defaultValue) const {
            return parse_integer(get_value(scenario), defaultValue);
        }

    private:
        std::string _name;
        std::vector<std::string> _values;
    };

    class VarSet {
    public:
        Variable *add_variable(const std::string &name);

        /** Null when absent. */
        Variable *find_variable(const std::string &name);
        const Variable *find_variable(const std::string &name) const;

        /** The empty string when absent, as in the Java. */
        std::string get_value(const std::string &name) const;

        int count() const { return static_cast<int>(variables_.size()); }

    private:
        std::vector<Variable> variables_;
    };

    enum class AsphereType {
        Even,
        EvenA2,
        Odd,
    };

    class AsphericalData {
    public:
        AsphericalData(AsphereType asphere_type, int surface_number)
            : _asphere_type(asphere_type), _surface_number(surface_number) {}

        void add_data(double d) { _data.push_back(d); }

        int data_points() const { return static_cast<int>(_data.size()); }

        double data(int i) const {
            return i >= 0 && i < static_cast<int>(_data.size())
                       ? _data[static_cast<std::size_t>(i)]
                       : 0.0;
        }

        int get_surface_number() const { return _surface_number; }
        AsphereType get_asphere_type() const { return _asphere_type; }
        bool is_odd_asphere() const { return _asphere_type == AsphereType::Odd; }

        std::vector<double> get_coeffs() const;

        double get_cc() const { return data(1); }
        double get_r() const { return data(0); }

    private:
        AsphereType _asphere_type;
        int _surface_number;
        std::vector<double> _data;
    };

    enum class SurfaceType {
        surface,
        aperture_stop,
        field_stop,
    };

    class LensSurface {
    public:
        explicit LensSurface(int id) : _id(id) {}

        SurfaceType get_surface_type() const { return _surface_type; }
        void set_surface_type(SurfaceType surface_type) { _surface_type = surface_type; }

        bool is_aperture_stop() const {
            return _surface_type == SurfaceType::aperture_stop;
        }

        double get_radius() const { return _radius; }
        void set_radius(double radius) { _radius = radius; }

        double get_thickness(int scenario) const;
        void add_thickness(double thickness) { _thickness_by_scenario.push_back(thickness); }

        double get_diameter(int scenario) const;
        void set_diameter(double value) { _diameter_by_scenario.push_back(value); }

        double get_refractive_index() const { return _refractive_index; }
        void set_refractive_index(double v) { _refractive_index = v; }

        double get_abbe_vd() const { return _abbe_vd; }
        void set_abbe_vd(double v) { _abbe_vd = v; }

        /** Null when the surface is not aspheric; borrowed from the spec. */
        const AsphericalData *get_aspherical_data() const { return _aspherical_data; }
        void set_aspherical_data(const AsphericalData *d) { _aspherical_data = d; }

        int get_id() const { return _id; }

        bool is_cover_glass() const { return _is_cover_glass; }
        void set_is_cover_glass(bool v) { _is_cover_glass = v; }

        /** Nullable in the Java; absent means the field was not present. */
        void set_glass_name(const std::string &name) { _glass_name = name; }
        const std::optional<std::string> &get_glass_name() const { return _glass_name; }

        void set_catalog_name(const std::optional<std::string> &name) {
            _catalog_name = name;
        }
        const std::optional<std::string> &get_catalog_name() const {
            return _catalog_name;
        }

        const std::vector<double> &get_thickness_by_scenario() const {
            return _thickness_by_scenario;
        }
        const std::vector<double> &get_diameter_by_scenario() const {
            return _diameter_by_scenario;
        }

    private:
        int _id;
        SurfaceType _surface_type = SurfaceType::surface;
        double _radius = 0;
        std::vector<double> _thickness_by_scenario;
        std::vector<double> _diameter_by_scenario;
        double _refractive_index = 0;
        double _abbe_vd = 0;
        bool _is_cover_glass = false;
        const AsphericalData *_aspherical_data = nullptr;
        std::optional<std::string> _glass_name;
        std::optional<std::string> _catalog_name;
    };

    class LensSpecifications {
    public:
        bool parse_file(const std::string &file_name);
        bool parse_buffer(const std::string &buffer);
        bool parse_lines(const std::vector<std::string> &lines);

        Variable *find_variable(const std::string &name) {
            return variables_.find_variable(name);
        }
        const Variable *find_variable(const std::string &name) const {
            return variables_.find_variable(name);
        }

        bool has_constant(const std::string &c) const {
            return constants_.find_variable(c) != nullptr;
        }

        double get_image_height() const;
        double get_focal_length() const;
        double get_focal_length(int scenario) const;
        double get_stop_diameter(int scenario) const;
        double get_angle_of_view_in_degrees(int scenario) const;
        double get_f_number(int scenario) const;
        double get_half_angle_of_view_in_radians(int scenario) const;
        double get_half_angle_of_view_in_degrees(int scenario) const;

        const VarSet &get_descriptive_data() const { return descriptive_data_; }
        VarSet &get_descriptive_data() { return descriptive_data_; }
        const std::vector<LensSurface> &get_surfaces() const { return surfaces_; }
        const std::vector<std::unique_ptr<AsphericalData>> &get_aspherical_data() const {
            return aspherical_data_;
        }
        const VarSet &get_patent_info() const { return patent_info_; }
        VarSet &get_patent_info() { return patent_info_; }
        const VarSet &get_report_data() const { return report_data_; }
        VarSet &get_report_data() { return report_data_; }

    private:
        static std::vector<std::string> splitLine(const std::string &line);

        LensSurface *find_surface(int id);

        void parse_thickness(const std::string &value, LensSurface &surface_builder);
        void parse_diameter(const std::string &value, bool isApertureStop,
                            LensSurface &surface_builder);

        VarSet descriptive_data_;
        VarSet variables_;
        std::vector<LensSurface> surfaces_;
        /**
         * Held by unique_ptr because LensSurface keeps a bare pointer to an
         * element: a vector of values would move them on the next push_back and
         * leave that pointer dangling.
         */
        std::vector<std::unique_ptr<AsphericalData>> aspherical_data_;
        VarSet constants_;
        VarSet patent_info_;
        VarSet report_data_;
    };
};

} // namespace redukti::importers

#endif // REDUKTI_IMPORTERS_OPTICALBENCHDATAIMPORTER_H
