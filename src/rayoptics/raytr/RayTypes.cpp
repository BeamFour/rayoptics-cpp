// Out-of-line members of the raytr data types.
#include "redukti/rayoptics/raytr/RayTypes.h"

#include "redukti/Text.h"

namespace redukti::rayoptics::raytr {

std::string RaySeg::toString() const {
    return "RaySeg(p=" + p.toString() + ", d=" + d.toString() +
           ", dst=" + doubleToString(dst) + ", nrml=" + nrml.toString() + ")";
}

RayPkg::RayPkg(std::vector<RaySeg> ray_, double op_delta_, double wvl_,
               const specs::Field *fld_, std::optional<mathlib::Vector2> input_pupil_,
               std::optional<mathlib::Vector2> vig_pupil_)
    : ray(std::move(ray_)), op_delta(op_delta_), wvl(wvl_),
      input_pupil(input_pupil_), vig_pupil(vig_pupil_) {
    // Java snapshots the Field into a ReadOnlyField here rather than aliasing
    // it, so later mutation of the Field does not change a traced ray.
    fld = fld_ != nullptr ? std::make_shared<const specs::ReadOnlyField>(*fld_)
                          : nullptr;
}

std::shared_ptr<const RayPkg> RayPkg::with(
    const specs::Field *fld_, std::optional<mathlib::Vector2> input_pupil_,
    std::optional<mathlib::Vector2> vig_pupil_) const {
    return std::make_shared<const RayPkg>(this->ray, this->op_delta, this->wvl, fld_,
                                          input_pupil_, vig_pupil_);
}

std::string RayPkg::toString() const {
    std::string s = "RayPkg(ray=[";
    for (std::size_t i = 0; i < ray.size(); i++) {
        if (i > 0)
            s += ", ";
        s += ray[i].toString();
    }
    s += "], op_delta=" + doubleToString(op_delta) + ", wvl=" + doubleToString(wvl) + ")";
    return s;
}

} // namespace redukti::rayoptics::raytr
