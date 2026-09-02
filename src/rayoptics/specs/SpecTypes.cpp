// C++ port of the small specs types.
#include "redukti/rayoptics/specs/SpecTypes.h"

#include "redukti/Text.h"
#include "redukti/rayoptics/specs/WvlSpec.h"

#include <algorithm>
#include <cctype>

namespace redukti::rayoptics::specs {

namespace {
bool equalsIgnoreCase(const std::string &a, const char *b) {
    std::string bs(b);
    if (a.size() != bs.size())
        return false;
    for (std::size_t i = 0; i < a.size(); i++) {
        if (std::tolower(static_cast<unsigned char>(a[i])) !=
            std::tolower(static_cast<unsigned char>(bs[i])))
            return false;
    }
    return true;
}
} // namespace

const char *name(SpecType v) {
    switch (v) {
    case SpecType::Field:
        return "Field";
    case SpecType::Aperture:
        return "Aperture";
    }
    return "";
}

const char *name(ImageKey v) {
    switch (v) {
    case ImageKey::Object:
        return "Object";
    case ImageKey::Image:
        return "Image";
    }
    return "";
}

const char *name(ValueKey v) {
    switch (v) {
    case ValueKey::RealHeight:
        return "RealHeight";
    case ValueKey::Height:
        return "Height";
    case ValueKey::Angle:
        return "Angle";
    case ValueKey::Fnum:
        return "Fnum";
    case ValueKey::NA:
        return "NA";
    case ValueKey::EPD:
        return "EPD";
    case ValueKey::Slope:
        return "Slope";
    case ValueKey::PUPIL:
        return "PUPIL";
    }
    return "";
}

const char *name(ConjugateType v) {
    switch (v) {
    case ConjugateType::FINITE:
        return "FINITE";
    case ConjugateType::INFINITE:
        return "INFINITE";
    }
    return "";
}

std::string SpecKey::toString() const {
    return std::string("SpecKey(type=") + name(type) + ", imageKey=" + name(imageKey) +
           ", valueKey=" + name(valueKey) + ")";
}

WvlWt::WvlWt(const std::string &wvl_name, double wt_)
    : wvl(WvlSpec::get_wavelength(wvl_name)), wt(wt_) {}

std::string Coord::toString() const {
    return "Coord(pt=" + pt.toString() + ", dir=" + dir.toString() + ")";
}

double SystemSpec::nm_to_sys_units(double nm) const {
    if (equalsIgnoreCase(dimensions, "m"))
        return 1e-9 * nm;
    else if (equalsIgnoreCase(dimensions, "cm"))
        return 1e-7 * nm;
    else if (equalsIgnoreCase(dimensions, "mm"))
        return 1e-6 * nm;
    else if (equalsIgnoreCase(dimensions, "in"))
        return 1e-6 * nm / 25.4;
    else if (equalsIgnoreCase(dimensions, "ft"))
        return 1e-6 * nm / 304.8;
    else
        return nm;
}

std::string FocusRange::toString() const {
    return "FocusRange{focus_shift=" + doubleToString(focus_shift) +
           ", defocus_range=" + doubleToString(defocus_range) + "}";
}

void FocusRange::list_str(std::string &sb) const {
    sb += "focus shift=" + doubleToString(focus_shift);
    if (defocus_range != 0.)
        sb += ", defocus range=" + doubleToString(defocus_range);
    sb += "\n";
}

} // namespace redukti::rayoptics::specs
