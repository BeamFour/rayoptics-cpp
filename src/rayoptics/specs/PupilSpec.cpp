// C++ port of org.redukti.rayoptics.specs.PupilSpec
#include "redukti/rayoptics/specs/PupilSpec.h"

#include "redukti/Text.h"
#include "redukti/rayoptics/parax/ParaxTypes.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"

namespace redukti::rayoptics::specs {

const std::vector<std::vector<double>> &PupilSpec::default_pupil_rays() {
    static const std::vector<std::vector<double>> v = {
        {0., 0.}, {1., 0.}, {-1., 0.}, {0., 1.}, {0., -1.}};
    return v;
}

const std::vector<std::string> &PupilSpec::default_ray_labels() {
    static const std::vector<std::string> v = {"00", "+X", "-X", "+Y", "-Y"};
    return v;
}

PupilSpec::PupilSpec(OpticalSpecs *parent,
                     std::optional<util::Pair<ImageKey, ValueKey>> k,
                     std::optional<double> value_)
    : key(SpecType::Aperture, ImageKey::Object, ValueKey::EPD) {
    util::Pair<ImageKey, ValueKey> kk =
        k.has_value() ? *k : util::Pair<ImageKey, ValueKey>(ImageKey::Object,
                                                            ValueKey::EPD);
    double vv = value_.has_value() ? *value_ : 1.0;
    this->optical_spec = parent;
    this->key = SpecKey(SpecType::Aperture, kk.first, kk.second);
    this->value = vv;
    this->pupil_rays = default_pupil_rays();
    this->ray_labels = default_ray_labels();
}

void PupilSpec::update_model() {
    if (pupil_rays.empty()) {
        pupil_rays = default_pupil_rays();
        ray_labels = default_ray_labels();
    }
}

void PupilSpec::apply_scale_factor(double scale_factor) {
    auto value_key = key.valueKey;
    if (value_key == ValueKey::EPD || value_key == ValueKey::PUPIL)
        value *= scale_factor;
}

util::Triple<ImageKey, std::optional<ValueKey>, double> PupilSpec::derive_parax_params() const {
    auto pupil_oi_key = key.imageKey;
    auto pupil_value_key = key.valueKey;
    auto pupil_value = this->value;
    std::optional<ValueKey> pupil_key;  // null in the Java unless a branch sets it
    if (ValueKey::NA == pupil_value_key) {
        auto obj_img_rindx = optical_spec->obj_img_rindex();
        auto n_obj = obj_img_rindx.first;
        auto n_img = obj_img_rindx.second;
        auto na = pupil_value;
        auto n = pupil_oi_key == ImageKey::Object ? n_obj : n_img;
        auto slope = parax::Etendue::na2slp(na, n);
        pupil_key = ValueKey::Slope;
        pupil_value = slope;
    } else if (ValueKey::Fnum == pupil_value_key) {
        auto fno = pupil_value;
        auto slope = -1.0 / (2.0 * fno);
        pupil_key = ValueKey::Slope;
        pupil_value = slope;
    } else if (ValueKey::EPD == pupil_value_key) {
        auto height = pupil_value / 2.0;
        pupil_key = ValueKey::Height;
        pupil_value = height;
    }
    return util::Triple<ImageKey, std::optional<ValueKey>, double>(
        pupil_oi_key, pupil_key, pupil_value);
}

std::string PupilSpec::toString() const {
    return "PupilSpec(key=" + key.toString() + ", value=" + doubleToString(value) + ")";
}

void PupilSpec::list_str(std::string &sb) const {
    sb += std::string(name(key.type)) + ": " + name(key.imageKey) + " " +
          name(key.valueKey) + "; value = " + doubleToString(value) + "\n";
}

} // namespace redukti::rayoptics::specs
