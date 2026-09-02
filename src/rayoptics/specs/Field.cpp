// Out-of-line members of Field and ReadOnlyField. Defined here rather than in
// the header because their shared_ptr members point at raytr types that are
// only forward-declared there.
#include "redukti/rayoptics/specs/Field.h"

#include "redukti/Text.h"
#include "redukti/rayoptics/raytr/RayTypes.h"
#include "redukti/rayoptics/specs/FieldSpec.h"

#include <cstdio>

namespace redukti::rayoptics::specs {

Field::Field(FieldSpec *fov_) : fov(fov_) {}

Field::~Field() = default;

void Field::update() {
    aim_info = std::nullopt;
    z_enp = std::nullopt;
    chief_ray = nullptr;
    ref_sphere = nullptr;
}

bool Field::is_relative() const {
    if (fov != nullptr)
        return fov->is_relative;
    return false;
}

double Field::max_field() const {
    if (fov != nullptr)
        return fov->value;
    return 1.;
}

std::string Field::toString() const {
    return "Field(x=" + doubleToString(x) + ", y=" + doubleToString(y) + ")";
}

void Field::list_str(std::string &sb, const std::string &fmtstr) const {
    char buf[256];
    if (fmtstr == "x") {
        std::snprintf(buf, sizeof(buf),
                      "x =%7.3f (%5.2f) vlx=%6.3f vux=%6.3f vly=%6.3f vuy=%6.3f", xv(),
                      xf(), vlx, vux, vly, vuy);
        sb += buf;
    } else if (fmtstr == "y") {
        std::snprintf(buf, sizeof(buf),
                      "y =%7.3f (%5.2f) vlx=%6.3f vux=%6.3f vly=%6.3f vuy=%6.3f", yv(),
                      yf(), vlx, vux, vly, vuy);
        sb += buf;
    } else if (fmtstr.empty()) {
        std::snprintf(buf, sizeof(buf),
                      "x,y=%5.2f vlx=%6.3f vux=%6.3f vly=%6.3f vuy=%6.3f", yv(), vlx,
                      vux, vly, vuy);
        sb += buf;
    } else {
        std::snprintf(buf, sizeof(buf),
                      "xy=(%7.3f, %7.3f) (%5.2f, %5.2f) vlx=%6.3f vux=%6.3f vly=%6.3f "
                      "vuy=%6.3f",
                      xv(), yv(), xf(), yf(), vlx, vux, vly, vuy);
        sb += buf;
    }
    if (aim_info.has_value()) {
        sb += " aim_info: [" + doubleToString((*aim_info)[0]) + "," +
              doubleToString((*aim_info)[1]) + "]";
    }
    if (z_enp.has_value()) {
        sb += " z_enp: " + doubleToString(*z_enp);
    }
    sb += "\n";

}

ReadOnlyField::ReadOnlyField(const Field &fld)
    : x(fld.x), y(fld.y), vux(fld.vux), vuy(fld.vuy), vlx(fld.vlx), vly(fld.vly),
      wt(fld.wt), z_enp(fld.z_enp), chief_ray(fld.chief_ray),
      ref_sphere(fld.ref_sphere), fov(fld.fov) {
    if (fld.aim_info.has_value())
        aim_info = mathlib::Vector2((*fld.aim_info)[0], (*fld.aim_info)[1]);
    else
        aim_info = std::nullopt;
}

ReadOnlyField::~ReadOnlyField() = default;

} // namespace redukti::rayoptics::specs
