// C++ port of org.redukti.rayoptics.specs.FieldSpec
#include "redukti/rayoptics/specs/FieldSpec.h"

#include "redukti/Text.h"
#include "redukti/mathlib/M.h"
#include "redukti/mathlib/Matrix3.h"
#include "redukti/rayoptics/parax/ParaxTypes.h"
#include "redukti/rayoptics/raytr/Wideangle.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/rayoptics/util/Lists.h"

#include <cmath>
#include <cstdio>

namespace redukti::rayoptics::specs {

using mathlib::Matrix3;
using mathlib::Vector3;

FieldSpec::FieldSpec(OpticalSpecs *parent,
                     std::optional<util::Pair<ImageKey, ValueKey>> key_,
                     std::optional<double> value_,
                     std::optional<std::vector<double>> flds,
                     std::optional<bool> is_relative_,
                     std::optional<bool> is_wide_angle_, std::optional<bool> do_init)
    : key(SpecType::Field, ImageKey::Object, ValueKey::Angle) {
    util::Pair<ImageKey, ValueKey> kk =
        key_.has_value() ? *key_
                         : util::Pair<ImageKey, ValueKey>(ImageKey::Object,
                                                          ValueKey::Angle);
    double vv = value_.has_value() ? *value_ : 0.0;
    bool rel = is_relative_.has_value() ? *is_relative_ : false;
    bool wide = is_wide_angle_.has_value() ? *is_wide_angle_ : false;
    bool init = do_init.has_value() ? *do_init : true;

    optical_spec = parent;
    this->key = SpecKey(SpecType::Field, kk.first, kk.second);
    this->value = vv;
    this->is_relative = rel;
    this->is_wide_angle = wide;
    if (init) {
        std::vector<double> f = flds.has_value() ? *flds : std::vector<double>{0., 1.};
        set_from_list(f);
    } else {
        fields.clear();
    }
}

void FieldSpec::set_from_list(const std::vector<double> &flds) {
    fields.clear();
    fields.reserve(flds.size());
    for (std::size_t i = 0; i < flds.size(); i++) {
        auto f = std::make_unique<Field>(this);
        f->y = flds[i];
        fields.push_back(std::move(f));
    }
    value = max_field().first;
}

util::Triple<ImageKey, std::optional<ValueKey>, double>
FieldSpec::derive_parax_params() const {
    auto fov_oi_key = key.imageKey;
    auto fov_value_key = key.valueKey;
    auto fov_value = this->value != 0 ? this->value : 1.0;
    std::optional<ValueKey> field_key; // null in the Java unless a branch sets it
    double field_value = 0.0;
    if (ValueKey::Angle == fov_value_key) {
        auto slope_bar = parax::Etendue::ang2slp(fov_value);
        field_key = ValueKey::Slope;
        field_value = slope_bar;
    } else if (ValueKey::Height == fov_value_key) {
        auto height_bar = fov_value;
        field_key = ValueKey::Height;
        field_value = height_bar;
    } else if (ValueKey::RealHeight == fov_value_key) {
        auto height_bar = fov_value;
        field_key = ValueKey::Height;
        field_value = height_bar;
    }
    return util::Triple<ImageKey, std::optional<ValueKey>, double>(fov_oi_key, field_key,
                                                                   field_value);
}

bool FieldSpec::check_is_wide_angle(double angle_threshold) {
    is_wide_angle = false;
    if (key.imageKey == ImageKey::Image && key.valueKey == ValueKey::RealHeight &&
        optical_spec->conjugate_type(ImageKey::Object) == ConjugateType::INFINITE) {
        is_wide_angle = true;
    } else if (key.imageKey == ImageKey::Object && key.valueKey == ValueKey::Angle) {
        auto max_angle = max_field().second;
        is_wide_angle = (max_angle > angle_threshold);
    }
    return is_wide_angle;
}

void FieldSpec::update_model() {
    for (auto &f : fields) {
        f->update();
    }
    double field_norm;
    if (is_relative)
        field_norm = 1.0;
    else
        field_norm = (value == 0.0) ? 1.0 : 1.0 / value;
    std::vector<std::string> labels;
    char buf[64];
    for (auto &f : fields) {
        std::string fldx, fldy;
        if (f->x != 0.0) {
            std::snprintf(buf, sizeof(buf), "%5.2fx", field_norm * f->x);
            fldx = buf;
        }
        if (f->y != 0.0) {
            std::snprintf(buf, sizeof(buf), "%5.2fy", field_norm * f->y);
            fldy = buf;
        }
        labels.push_back(fldx + fldy);
    }
    labels[0] = "axis";
    if (labels.size() > 1)
        util::Lists::set(labels, -1, std::string("edge"));
    this->index_labels = labels;
}

Field *FieldSpec::with_index_label(const std::string &label) {
    for (std::size_t i = 0; i < index_labels.size(); i++) {
        if (index_labels[i] == label)
            return fields[i].get();
    }
    return nullptr;
}

void FieldSpec::apply_scale_factor(double scale_factor) {
    ValueKey value_key = key.valueKey;
    if (value_key == ValueKey::Height) {
        if (!is_relative) {
            for (auto &f : fields) {
                f->apply_scale_factor(scale_factor);
            }
        }
        value *= scale_factor;
    }
}

Coord FieldSpec::obj_coords(Field &fld) {
    ImageKey obj_img_key = key.imageKey;
    ValueKey value_key = key.valueKey;
    Vector3 fld_coord(fld.x, fld.y, 0.0);
    Vector3 rel_fld_coord(fld.x, fld.y, 0.0);
    if (is_relative)
        fld_coord = fld_coord.times(value);
    else if (value != 0.0)
        rel_fld_coord = rel_fld_coord.divide(value);
    auto opt_model = optical_spec->opt_model;
    auto &pr = optical_spec->parax_data->pr_ray;
    auto &fod = optical_spec->parax_data->fod;
    Vector3 obj_pt = Vector3::ZERO;
    Vector3 obj_dir = Vector3::ZERO;
    auto obj2enp_dist = fod.obj_dist + fod.enp_dist;
    Vector3 pt1(0.0, 0.0, obj2enp_dist);
    ConjugateType obj_conj = optical_spec->conjugate_type(ImageKey::Object);
    if (obj_conj == ConjugateType::INFINITE) {
        Vector3 fld_angle = Vector3::ZERO;
        if (obj_img_key == ImageKey::Image) {
            double max_field_ang;
            if (value_key == ValueKey::RealHeight) {
                double wvl = optical_spec->wvls->central_wvl();
                auto pkg = raytr::Wideangle::eval_real_image_ht(opt_model, fld, wvl);
                obj_pt = pkg.ray_data.pt;
                obj_dir = pkg.ray_data.dir;
                if (is_wide_angle) {
                    fld.z_enp = pkg.z_enp;
                    fld.aim_info = std::nullopt;
                } else {
                    auto del_z = fod.enp_dist - pkg.z_enp;
                    std::vector<double> aim_pt;
                    if (mathlib::M::isZero(obj_dir.z)) {
                        aim_pt = {0.0, 0.0};
                    } else {
                        aim_pt = {(obj_dir.x / obj_dir.z) * del_z,
                                  (obj_dir.y / obj_dir.z) * del_z};
                    }
                    fld.aim_info = aim_pt;
                    fld.z_enp = std::nullopt;
                }
                return Coord(obj_pt, obj_dir);
            } else {
                max_field_ang = std::atan(pr[0].slp);
                fld_angle = rel_fld_coord.times(max_field_ang);
            }
        } else {
            if (value_key == ValueKey::Angle) {
                fld_angle = fld_coord.deg2rad();
            } else {
                obj_pt = fld_coord;
                obj_dir = pt1.minus(obj_pt).normalize();
                return Coord(obj_pt, obj_dir);
            }
        }
        auto ang_x = fld_angle.x;
        auto ang_y = fld_angle.y;
        Vector3 dir_cos(std::sin(ang_x) * std::cos(ang_y), std::sin(ang_y),
                        std::cos(ang_x) * std::cos(ang_y));
        if (is_wide_angle) {
            auto rot_mat = Matrix3::rot_v1_into_v2(Vector3::vector3_001, dir_cos);
            obj_pt = rot_mat.multiply(pt1.negate()).plus(pt1);
        } else {
            obj_pt = Vector3(dir_cos.x / dir_cos.z, dir_cos.y / dir_cos.z, 0.0)
                         .times(obj2enp_dist);
        }
        obj_dir = dir_cos;
    } else if (obj_conj == ConjugateType::FINITE) {
        if (obj_img_key == ImageKey::Image) {
            if (value_key == ValueKey::RealHeight) {
                auto wvl = optical_spec->wvls->central_wvl();
                auto pkg = raytr::Wideangle::eval_real_image_ht(opt_model, fld, wvl);
                obj_pt = pkg.ray_data.pt;
                obj_dir = pkg.ray_data.dir;
                auto z_enp = pkg.z_enp;
                if (is_wide_angle) {
                    fld.z_enp = z_enp;
                } else { // compute offset at paraxial entrance pupil
                    auto del_z = fod.enp_dist - z_enp;
                    std::vector<double> aim_pt;
                    if (mathlib::M::isZero(obj_dir.z))
                        aim_pt = {0., 0.};
                    else {
                        aim_pt = {del_z * obj_dir.x / obj_dir.z,
                                  del_z * obj_dir.y / obj_dir.z};
                    }
                    fld.aim_info = aim_pt;
                }
                return Coord(obj_pt, obj_dir);
            } else {
                auto max_field_ht = pr[0].ht;
                obj_pt = rel_fld_coord.times(max_field_ht);
            }
        } else { // obj_img_key == 'object'
            if (value_key == ValueKey::Angle) {
                auto fld_angle = fld_coord.deg2rad();
                obj_dir = fld_angle.sin();
                auto z = std::sqrt(1.0 - obj_dir.x * obj_dir.x - obj_dir.y * obj_dir.y);
                obj_dir = Vector3(obj_dir.x, obj_dir.y, z);
                obj_pt = Vector3(obj_dir.x / obj_dir.z, obj_dir.y / obj_dir.z, 0.0)
                             .times(obj2enp_dist);
                return Coord(obj_pt, obj_dir);
            } else {
                obj_pt = fld_coord;
            }
        }
        obj_dir = pt1.minus(obj_pt).normalize();
    }
    return Coord(obj_pt, obj_dir);
}

util::Pair<double, int> FieldSpec::max_field() const {
    int max_fld = 0;
    double max_fld_sqrd = -1.0;
    for (std::size_t i = 0; i < fields.size(); i++) {
        const Field &f = *fields[i];
        double fld_sqrd = f.x * f.x + f.y * f.y;
        if (fld_sqrd > max_fld_sqrd) {
            max_fld_sqrd = fld_sqrd;
            max_fld = static_cast<int>(i);
        }
    }
    double max_fld_value = std::sqrt(max_fld_sqrd);
    if (is_relative)
        max_fld_value *= value;
    return util::Pair<double, int>(max_fld_value, max_fld);
}

void FieldSpec::clear_vignetting() {
    for (auto &f : fields) {
        f->clear_vignetting();
    }
}

std::string FieldSpec::toString() const {
    return "FieldSpec(key=" + key.toString() +
           ", max field=" + doubleToString(max_field().first) +
           ", is wide angle=" + (is_wide_angle ? "true" : "false") + ")";
}

void FieldSpec::list_str(std::string &sb) const {
    sb += std::string(name(key.type)) + ": " + name(key.imageKey) + " " +
          name(key.valueKey) + "; value = " + doubleToString(value) + "\n";
    bool has_x = false;
    bool has_y = false;
    std::string fmtstr = "";
    for (auto &fld : fields) {
        if (fld->x != 0. && fld->y != 0.) {
            has_y = true;
            has_x = true;
        } else if (fld->x == 0.0 && fld->y != 0.) {
            has_y = true;
            fmtstr = "y";
        } else if (fld->x != 0. && fld->y == 0.) {
            has_x = true;
            fmtstr = "x";
        }
    }
    if (has_x && has_y)
        fmtstr = "xy";
    for (auto &fld : fields)
        fld->list_str(sb, fmtstr);
    sb += std::string("is_relative=") + (is_relative ? "true" : "false") +
          ", is_wideangle=" + (is_wide_angle ? "true" : "false") + "\n";
}

} // namespace redukti::rayoptics::specs
