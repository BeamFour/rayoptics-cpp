// C++ port of the small specs types:
//   SpecType, ImageKey, ValueKey, ConjugateType, SpecKey, WvlWt, Coord,
//   SystemSpec, FocusRange.
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_SPECS_SPECTYPES_H
#define REDUKTI_RAYOPTICS_SPECS_SPECTYPES_H

#include "redukti/mathlib/Vector3.h"

#include <string>

namespace redukti::rayoptics::specs {

enum class SpecType {
    Field,
    Aperture,
};

enum class ImageKey {
    Object,
    Image,
};

enum class ValueKey {
    RealHeight,
    Height,
    Angle,
    Fnum,
    NA,
    EPD,
    Slope,
    PUPIL,
};

enum class ConjugateType {
    FINITE,
    INFINITE,
};

const char *name(SpecType v);
const char *name(ImageKey v);
const char *name(ValueKey v);
const char *name(ConjugateType v);

class SpecKey {
public:
    SpecType type;
    ImageKey imageKey;
    ValueKey valueKey;

    SpecKey(SpecType type_, ImageKey imageKey_, ValueKey valueKey_)
        : type(type_), imageKey(imageKey_), valueKey(valueKey_) {}

    std::string toString() const;
};

class WvlWt {
public:
    double wvl;
    double wt;

    WvlWt(double wvl_, double wt_) : wvl(wvl_), wt(wt_) {}
    /** Looks the wavelength up by spectral line name, e.g. "d". */
    WvlWt(const std::string &wvl_name, double wt_);
};

class Coord {
public:
    mathlib::Vector3 pt;
    mathlib::Vector3 dir;

    Coord(const mathlib::Vector3 &pt_, const mathlib::Vector3 &dir_)
        : pt(pt_), dir(dir_) {}

    std::string toString() const;
};

class SystemSpec {
public:
    std::string title;
    std::string initials;
    std::string dimensions;
    double temperature;
    double pressure;

    SystemSpec()
        : title(""), initials(""), dimensions("mm"), temperature(20.0), pressure(760.0) {}

    /** Convert nm to the system units. */
    double nm_to_sys_units(double nm) const;
};

/**
 * Focus range specification.
 *   focus_shift: focus shift (z displacement) from nominal image interface
 *   defocus_range: +/- half the total focal range, from the focus_shift position
 */
class FocusRange {
public:
    double focus_shift;
    double defocus_range;

    FocusRange(double focus_shift_, double defocus_range_)
        : focus_shift(focus_shift_), defocus_range(defocus_range_) {}

    FocusRange() : FocusRange(0.0, 0.0) {}

    std::string toString() const;

    void apply_scale_factor(double scale_factor) {
        focus_shift *= scale_factor;
        defocus_range *= scale_factor;
    }

    double get_focus(double fr) const { return focus_shift + fr * defocus_range; }
    double get_focus() const { return get_focus(0.0); }

    void list_str(std::string &sb) const;
};

} // namespace redukti::rayoptics::specs

#endif // REDUKTI_RAYOPTICS_SPECS_SPECTYPES_H
