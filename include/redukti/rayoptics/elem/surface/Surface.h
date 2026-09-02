// C++ port of org.redukti.rayoptics.elem.surface.Surface
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_ELEM_SURFACE_SURFACE_H
#define REDUKTI_RAYOPTICS_ELEM_SURFACE_SURFACE_H

#include "redukti/rayoptics/elem/profiles/Spherical.h"
#include "redukti/rayoptics/elem/surface/Aperture.h"
#include "redukti/rayoptics/seq/Interface.h"

#include <memory>
#include <string>
#include <vector>

namespace redukti::rayoptics::elem::surface {

class Surface : public seq::Interface {
public:
    std::string label;
    /** Apertures are shared with the element model, so held by pointer. */
    std::vector<std::shared_ptr<Aperture>> clear_apertures;
    std::vector<std::shared_ptr<Aperture>> edge_apertures;

    Surface() : Surface("", seq::InteractMode::TRANSMIT) {}

    Surface(std::string label_, seq::InteractMode interact_mode_)
        : Surface(interact_mode_, 0.0, 1.0, nullptr, std::move(label_), nullptr) {}

    Surface(seq::InteractMode interact_mode_, double delta_n_, double max_ap,
            std::shared_ptr<DecenterData> decenter_, std::string label_,
            std::shared_ptr<profiles::SurfaceProfile> profile_)
        : seq::Interface(interact_mode_, delta_n_, max_ap, std::move(decenter_)),
          label(std::move(label_)) {
        if (profile_)
            this->profile = std::move(profile_);
        else
            this->profile = std::make_shared<profiles::Spherical>();
    }

    void update() override {
        seq::Interface::update();
        profile->update();
    }

    IntersectionResult intersect(const mathlib::Vector3 &p0, const mathlib::Vector3 &d,
                                 double eps, util::ZDir z_dir) const override {
        return profile->intersect(p0, d, eps, z_dir);
    }

    mathlib::Vector3 normal(const mathlib::Vector3 &p) const override {
        return profile->normal(p);
    }

    std::string toString() const override;

    double profile_cv() const override { return profile->cv; }

    double optical_power() const override { return delta_n * profile->cv; }

    void set_optical_power(double pwr) override {
        profile->cv = delta_n != 0.0 ? pwr / delta_n : 0.0;
    }

    void set_optical_power(double pwr, double n_before, double n_after) override {
        delta_n = n_after - n_before;
        set_optical_power(pwr);
    }

    void apply_scale_factor(double scale_factor) override {
        seq::Interface::apply_scale_factor(scale_factor);
        profile->apply_scale_factor(scale_factor);
        auto abs_scale_factor = std::abs(scale_factor);
        for (auto &e : edge_apertures)
            e->apply_scale_factor(abs_scale_factor);
        for (auto &ca : clear_apertures)
            ca->apply_scale_factor(abs_scale_factor);
    }

    std::vector<std::shared_ptr<Aperture>> get_ca_list() const {
        std::vector<std::shared_ptr<Aperture>> result;
        for (const auto &e : clear_apertures)
            if (!e->is_obscuration)
                result.push_back(e);
        return result;
    }

    double surface_od() const override {
        double od = 0.0;
        if (!edge_apertures.empty()) {
            for (const auto &e : edge_apertures) {
                double edg = e->max_dimension();
                if (edg > od)
                    od = edg;
            }
        } else if (!get_ca_list().empty()) {
            for (const auto &ca : get_ca_list()) {
                double ap = ca->max_dimension();
                if (ap > od)
                    od = ap;
            }
        } else {
            od = max_aperture;
        }
        return od;
    }

    bool point_inside(double x, double y, std::optional<double> fuzz) const override {
        double f = fuzz.has_value() ? *fuzz : 1e-5;
        bool is_inside = true;
        if (clear_apertures.size() > 0) {
            for (const auto &ca : clear_apertures) {
                is_inside = is_inside && ca->point_inside(x, y, f);
                if (!is_inside)
                    return is_inside;
            }
        } else
            return seq::Interface::point_inside(x, y, f);
        return is_inside;
    }

    mathlib::Vector2 edge_pt_target(const mathlib::Vector2 &rel_dir) const override {
        auto ca_list = get_ca_list();
        if (!ca_list.empty())
            return ca_list[0]->edge_pt_target(rel_dir);
        else
            return seq::Interface::edge_pt_target(rel_dir);
    }

    void set_max_aperture(double max_ap) override {
        seq::Interface::set_max_aperture(max_ap);
        auto ca_list = get_ca_list();
        for (auto &ap : ca_list) {
            ap->set_dimension(max_ap, max_ap);
        }
    }
};

} // namespace redukti::rayoptics::elem::surface

#endif // REDUKTI_RAYOPTICS_ELEM_SURFACE_SURFACE_H
