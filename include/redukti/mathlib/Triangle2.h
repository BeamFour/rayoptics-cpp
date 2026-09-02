// C++ port of org.redukti.mathlib.Triangle2
//
// The software is ported from Goptical, hence is licensed under the GPL.
// Copyright (c) 2021 Dibyendu Majumdar
// Copyright (C) 2010-2011 Free Software Foundation, Inc; Author: Alexandre Becoulet
#ifndef REDUKTI_MATHLIB_TRIANGLE2_H
#define REDUKTI_MATHLIB_TRIANGLE2_H

#include "redukti/mathlib/Vector2.h"

#include <array>

namespace redukti::mathlib {

class Triangle2 {
public:
    Triangle2(const Vector2 &a, const Vector2 &b, const Vector2 &c) : _v{a, b, c} {}

    Vector2 get_centroid() const { return _v[0].plus(_v[1]).plus(_v[2]).divide(3.); }

    std::array<Vector2, 3> as_array() const { return _v; }

private:
    static const int N = 3;
    std::array<Vector2, N> _v;
};

} // namespace redukti::mathlib

#endif // REDUKTI_MATHLIB_TRIANGLE2_H
