// C++ port of org.redukti.mathlib.M
#include "redukti/mathlib/M.h"

namespace redukti::mathlib::M {

double cosd(double deg) { return std::cos(toRadians(deg)); }

double sind(double deg) { return std::sin(toRadians(deg)); }

double tand(double deg) { return std::tan(toRadians(deg)); }

} // namespace redukti::mathlib::M
