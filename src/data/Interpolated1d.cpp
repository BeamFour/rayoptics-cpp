// C++ port of org.redukti.data.Interpolated1d
#include "redukti/data/Interpolated1d.h"

#include "redukti/Exceptions.h"
#include "redukti/mathlib/M.h"
#include "redukti/mathlib/Vector2.h"
#include "redukti/util/ArrayIndex2D.h"

namespace redukti::data {

using mathlib::Vector2;
using mathlib::M::square;
using util::ArrayIndex2D;

void Interpolated1d::set_interpolation(Interpolation i) {
    this->_method = i;
    invalidate();
}

void Interpolated1d::resizePoly(int n) {
    for (int i = static_cast<int>(_poly.size()); i < n; i++)
        _poly.push_back(PolyS(0, 0, 0, 0));
}

double Interpolated1d::interpolate(double x, int d) {
    switch (_method) {
    case Interpolation::Nearest:
        if (_invalid) {
            _poly.clear();
            _invalid = false;
        }
        return interpolate_nearest(d, x);
    case Interpolation::Linear:
        if (_invalid) {
            _poly.clear();
            _invalid = false;
        }
        return interpolate_linear(d, x);
    case Interpolation::Quadratic:
        return _invalid ? update_quadratic(d, x) : interpolate_quadratic(d, x);
    case Interpolation::CubicSimple:
        return _invalid ? update_cubic_simple(d, x) : interpolate_cubic(d, x);
    case Interpolation::CubicDeriv:
        return _invalid ? update_cubic_deriv(d, x) : interpolate_cubic(d, x);
    case Interpolation::Cubic2Deriv:
        return _invalid ? update_cubic2_deriv(d, x) : interpolate_cubic(d, x);
    case Interpolation::CubicDerivInit:
        return _invalid ? update_cubic_deriv_init(d, x) : interpolate_cubic(d, x);
    case Interpolation::Cubic2DerivInit:
        return _invalid ? update_cubic2_deriv_init(d, x) : interpolate_cubic(d, x);
    case Interpolation::Cubic:
        return _invalid ? update_cubic(d, x) : interpolate_cubic(d, x);
    case Interpolation::Cubic2:
        return _invalid ? update_cubic2(d, x) : interpolate_cubic(d, x);
    default:
        throw IllegalStateException("invalid interpolation selected");
    }
}

void Interpolated1d::compute_cubic_2nd_deriv(cubic_2nd_deriv_init_e de, int n,
                                             std::vector<double> &dd, double d0,
                                             double dn) {
    std::vector<double> eq(static_cast<std::size_t>(n) * 3, 0.0);
    ArrayIndex2D idx(n, 3);
    auto un = static_cast<std::size_t>(n);
    switch (de) {
    case cubic_2nd_deriv_init_e::Cubic2ndDerivQuadratic:
        dd[0] = dd[un - 1] = 0.0;
        eq[static_cast<std::size_t>(idx.i(0, 1))] =
            eq[static_cast<std::size_t>(idx.i(n - 1, 1))] = 1.0;
        eq[static_cast<std::size_t>(idx.i(1, 0))] =
            eq[static_cast<std::size_t>(idx.i(n - 2, 2))] = -1.0;
        break;
    case cubic_2nd_deriv_init_e::Cubic2ndDerivFirst: {
        double x0 = _data_set->get_x_interval(0);
        double xn = _data_set->get_x_interval(n - 2);
        dd[0] = (_data_set->get_y_value(1) - _data_set->get_y_value(0)) / x0 - d0;
        dd[un - 1] =
            dn - (_data_set->get_y_value(n - 1) - _data_set->get_y_value(n - 2)) / xn;
        eq[static_cast<std::size_t>(idx.i(0, 1))] = x0 / 3.0;
        eq[static_cast<std::size_t>(idx.i(1, 0))] = x0 / 6.0;
        eq[static_cast<std::size_t>(idx.i(n - 2, 2))] = xn / 6.0;
        eq[static_cast<std::size_t>(idx.i(n - 1, 1))] = xn / 3.0;
        break;
    }
    case cubic_2nd_deriv_init_e::Cubic2ndDerivSecond:
        dd[0] = d0;
        dd[un - 1] = dn;
        eq[static_cast<std::size_t>(idx.i(1, 0))] =
            eq[static_cast<std::size_t>(idx.i(n - 2, 2))] = 0.0;
        eq[static_cast<std::size_t>(idx.i(0, 1))] =
            eq[static_cast<std::size_t>(idx.i(n - 1, 1))] = 1.0;
        break;
    }
    int i;
    for (i = 1; i < n - 1; i++) {
        eq[static_cast<std::size_t>(idx.i(i - 1, 2))] =
            _data_set->get_x_interval(i - 1) / 6.0;
        eq[static_cast<std::size_t>(idx.i(i, 1))] =
            _data_set->get_x_interval(i - 1, i + 1) / 3.0;
        eq[static_cast<std::size_t>(idx.i(i + 1, 0))] = _data_set->get_x_interval(i) / 6.0;
        dd[static_cast<std::size_t>(i)] =
            (_data_set->get_y_value(i + 1) - _data_set->get_y_value(i)) /
                _data_set->get_x_interval(i) -
            (_data_set->get_y_value(i) - _data_set->get_y_value(i - 1)) /
                _data_set->get_x_interval(i - 1);
    }
    for (i = 1; i < n; i++) {
        double f = eq[static_cast<std::size_t>(idx.i(i - 1, 2))] /
                   eq[static_cast<std::size_t>(idx.i(i - 1, 1))];
        eq[static_cast<std::size_t>(idx.i(i, 1))] -=
            f * eq[static_cast<std::size_t>(idx.i(i, 0))];
        dd[static_cast<std::size_t>(i)] -= f * dd[static_cast<std::size_t>(i - 1)];
    }
    double k = 0;
    for (i = n - 1; i >= 0; i--) {
        double ddi = (dd[static_cast<std::size_t>(i)] - k) /
                     eq[static_cast<std::size_t>(idx.i(i, 1))];
        dd[static_cast<std::size_t>(i)] = ddi;
        k = eq[static_cast<std::size_t>(idx.i(i, 0))] * ddi;
    }
}

void Interpolated1d::set_linear_poly(PolyS &p, double p1x, double p1y, double p2x,
                                     double p2y) {
    p.a = 0.0;
    p.b = 0.0;
    p.c = (p1y - p2y) / (p1x - p2x);
    p.d = (p2x * p1y - p1x * p2y) / (p2x - p1x);
}

void Interpolated1d::set_linear_poly(PolyS &p, double p1x, double p1y, double d1) {
    p.a = 0.0;
    p.b = 0.0;
    p.c = d1;
    p.d = p1y - d1 * p1x;
}

void Interpolated1d::set_quadratic_poly(PolyS &p, double p1x, double p1y, double p2x,
                                        double p2y, double p3x, double p3y) {
    double n = ((p2x - p1x) * (p3x - p1x) * (p3x - p2x));
    p.a = 0.0;
    p.b = (p3y * (p2x - p1x) + p2y * (p1x - p3x) + p1y * (p3x - p2x)) / n;
    double p1x2 = square(p1x);
    double p2x2 = square(p2x);
    double p3x2 = square(p3x);
    p.c = (p3y * (p1x2 - p2x2) + p2y * (p3x2 - p1x2) + p1y * (p2x2 - p3x2)) / n;
    p.d = (p3y * (p1x * p2x2 - p2x * p1x2) + p2y * (p3x * p1x2 - p1x * p3x2) +
           p1y * (p2x * p3x2 - p3x * p2x2)) /
          n;
}

void Interpolated1d::set_quadratic_poly(PolyS &p, double px, double py, double d,
                                        double dd) {
    p.a = 0;
    p.b = dd / 2.0;
    p.c = -px * dd + d;
    p.d = 0.5 * px * px * dd - px * d + py;
}

void Interpolated1d::set_cubic_poly(PolyS &p, double p1x, double p1y, double p2x,
                                    double p2y, double d1, double d2) {
    double x1 = p1x;
    double x2 = p2x;
    double y1 = p1y;
    double y2 = p2y;
    p.a = -(2. * y1 - 2. * y2 + (d2 + d1) * x2 - (d2 + d1) * x1) /
          (3. * x1 * x2 * x2 - x2 * x2 * x2 - 3. * x1 * x1 * x2 + x1 * x1 * x1);
    p.b = (x1 * ((d2 - d1) * x2 - 3. * y2) - 3. * x2 * y2 + (3. * x2 + 3. * x1) * y1 +
           (d2 + 2. * d1) * x2 * x2 - (2. * d2 + d1) * x1 * x1) /
          (3. * x1 * x2 * x2 - x2 * x2 * x2 - 3. * x1 * x1 * x2 + x1 * x1 * x1);
    p.c = -(x1 * ((2. * d2 + d1) * x2 * x2 - 6. * x2 * y2) + 6. * x1 * x2 * y1 +
            d1 * x2 * x2 * x2 - (d2 + 2. * d1) * x1 * x1 * x2 - d2 * x1 * x1 * x1) /
          (3. * x1 * x2 * x2 - x2 * x2 * x2 - 3. * x1 * x1 * x2 + x1 * x1 * x1);
    p.d = (x1 * x1 * ((d2 - d1) * x2 * x2 - 3. * x2 * y2) + x1 * x1 * x1 * (y2 - d2 * x2) +
           (3. * x1 * x2 * x2 - x2 * x2 * x2) * y1 + d1 * x1 * x2 * x2 * x2) /
          (3. * x1 * x2 * x2 - x2 * x2 * x2 - 3. * x1 * x1 * x2 + x1 * x1 * x1);
}

void Interpolated1d::set_cubic_poly2(PolyS &p, double p1x, double p1y, double p2x,
                                     double p2y, double dd1, double dd2) {
    p.a = (dd1 - dd2) / (6. * p1x - 6. * p2x);
    p.b = (dd2 * p1x - dd1 * p2x) / (2. * p1x - 2. * p2x);
    p.c = (6. * p1y - 6. * p2y + (dd2 + 2. * dd1) * p2x * p2x +
           (2. * dd1 - 2. * dd2) * p1x * p2x - (2. * dd2 + dd1) * p1x * p1x) /
          (6. * p1x - 6. * p2x);
    p.d = -(p1x * ((dd2 + 2. * dd1) * p2x * p2x - 6. * p2y) + 6. * p2x * p1y -
            (2. * dd2 + dd1) * p1x * p1x * p2x) /
          (6. * p1x - 6. * p2x);
}

double Interpolated1d::interpolate_nearest(int d, double x) {
    switch (d) {
    case 0:
        return _data_set->get_y_value(_data_set->get_nearest(x));
    default:
        return 0.0;
    }
}

double Interpolated1d::interpolate_linear(int d, double x) {
    int di = _data_set->get_interval(x);
    if (di == 0)
        di++;
    else if (di == _data_set->get_count())
        di--;
    switch (d) {
    case 0: {
        double mu =
            (x - _data_set->get_x_value(di - 1)) / (_data_set->get_x_interval(di - 1));
        return _data_set->get_y_value(di - 1) * (1.0 - mu) +
               _data_set->get_y_value(di) * mu;
    }
    case 1: {
        return (_data_set->get_y_value(di) - _data_set->get_y_value(di - 1)) /
               (_data_set->get_x_interval(di - 1));
    }
    default: {
        return 0.0;
    }
    }
}

double Interpolated1d::interpolate_quadratic(int d, double x) {
    const PolyS &p = _poly[static_cast<std::size_t>(_data_set->get_nearest(x))];
    switch (d) {
    case 0:
        return x * (p.b * x + p.c) + p.d;
    case 1:
        return 2.0 * p.b * x + p.c;
    case 2:
        return 2.0 * p.b;
    default:
        return 0.0;
    }
}

double Interpolated1d::update_quadratic(int d, double x) {
    if (_data_set->get_count() < 3)
        throw IllegalStateException("data set doesn't contains enough data");
    resizePoly(_data_set->get_count());
    set_linear_poly(_poly[0], _data_set->get_x_value(0), _data_set->get_y_value(0),
                    _data_set->get_x_value(1), _data_set->get_y_value(1));
    int i;
    for (i = 1; i < _data_set->get_count() - 1; i++) {
        double p1x = (_data_set->get_x_value(i - 1) + _data_set->get_x_value(i)) / 2.0;
        double p1y = (_data_set->get_y_value(i - 1) + _data_set->get_y_value(i)) / 2.0;
        double p3x = (_data_set->get_x_value(i) + _data_set->get_x_value(i + 1)) / 2.0;
        double p3y = (_data_set->get_y_value(i) + _data_set->get_y_value(i + 1)) / 2.0;
        set_quadratic_poly(_poly[static_cast<std::size_t>(i)], p1x, p1y,
                           _data_set->get_x_value(i), _data_set->get_y_value(i), p3x, p3y);
    }
    set_linear_poly(_poly[static_cast<std::size_t>(i)], _data_set->get_x_value(i - 1),
                    _data_set->get_y_value(i - 1), _data_set->get_x_value(i),
                    _data_set->get_y_value(i));
    _invalid = false;
    return interpolate_quadratic(d, x);
}

double Interpolated1d::interpolate_cubic(int d, double x) {
    const PolyS &p = _poly[static_cast<std::size_t>(_data_set->get_interval(x))];
    switch (d) {
    case 0:
        return ((p.a * x + p.b) * x + p.c) * x + p.d;
    case 1:
        return (3.0 * p.a * x + 2.0 * p.b) * x + p.c;
    case 2:
        return 6.0 * p.a * x + 2.0 * p.b;
    case 3:
        return 6.0 * p.a;
    default:
        return 0.0;
    }
}

double Interpolated1d::update_cubic_simple(int d, double x) {
    int n = _data_set->get_count();
    if (n < 4)
        throw IllegalStateException("data set doesn't contains enough data");
    resizePoly(n + 1);
    Vector2 vm1(_data_set->get_x_value(0), _data_set->get_y_value(0));
    Vector2 vm2 = vm1;
    Vector2 v(_data_set->get_x_value(1), _data_set->get_y_value(1));
    Vector2 vp1(_data_set->get_x_value(2), _data_set->get_y_value(2));
    double d1 = (v.y - vm1.y) / (v.x - vm1.x);
    double d2 = (vp1.y - vm1.y) / (vp1.x - vm1.x);
    set_linear_poly(_poly[0], vm1.x, vm1.y, d1);
    set_cubic_poly(_poly[1], vm1.x, vm1.y, v.x, v.y, d1, d2);
    for (int i = 2; i < n - 1; i++) {
        vm2 = vm1;
        vm1 = v;
        v = vp1;
        vp1 = Vector2(_data_set->get_x_value(i + 1), _data_set->get_y_value(i + 1));
        d1 = d2;
        d2 = (vp1.y - vm1.y) / (vp1.x - vm1.x);
        set_cubic_poly(_poly[static_cast<std::size_t>(i)], vm1.x, vm1.y, v.x, v.y, d1, d2);
    }
    (void)vm2; // assigned in the Java loop but never read, as here
    d1 = d2;
    d2 = (vp1.y - v.y) / (vp1.x - v.x);
    set_cubic_poly(_poly[static_cast<std::size_t>(n - 1)], v.x, v.y, vp1.x, vp1.y, d1, d2);
    set_linear_poly(_poly[static_cast<std::size_t>(n)], vp1.x, vp1.y, d2);
    _invalid = false;
    return interpolate_cubic(d, x);
}

double Interpolated1d::update_cubic(int d, double x) {
    int n = _data_set->get_count();
    if (n < 4)
        throw IllegalStateException("data set doesn't contains enough data");
    resizePoly(n + 1);
    double d0 =
        (_data_set->get_y_value(1) - _data_set->get_y_value(0)) / _data_set->get_x_interval(0);
    double dn = (_data_set->get_y_value(n - 1) - _data_set->get_y_value(n - 2)) /
                _data_set->get_x_interval(n - 2);
    std::vector<double> dd(static_cast<std::size_t>(n), 0.0);
    compute_cubic_2nd_deriv(cubic_2nd_deriv_init_e::Cubic2ndDerivFirst, n, dd, d0, dn);
    set_linear_poly(_poly[0], _data_set->get_x_value(0), _data_set->get_y_value(0), d0);
    for (int i = 1; i < n; i++)
        set_cubic_poly2(_poly[static_cast<std::size_t>(i)], _data_set->get_x_value(i - 1),
                        _data_set->get_y_value(i - 1), _data_set->get_x_value(i),
                        _data_set->get_y_value(i), dd[static_cast<std::size_t>(i - 1)],
                        dd[static_cast<std::size_t>(i)]);
    set_linear_poly(_poly[static_cast<std::size_t>(n)], _data_set->get_x_value(n - 1),
                    _data_set->get_y_value(n - 1), dn);
    _invalid = false;
    return interpolate_cubic(d, x);
}

double Interpolated1d::update_cubic2(int d, double x) {
    int n = _data_set->get_count();
    if (n < 4)
        throw IllegalStateException("data set doesn't contains enough data");
    resizePoly(n + 1);
    double d0 =
        (_data_set->get_y_value(1) - _data_set->get_y_value(0)) / _data_set->get_x_interval(0);
    double dn = (_data_set->get_y_value(n - 1) - _data_set->get_y_value(n - 2)) /
                _data_set->get_x_interval(n - 2);
    std::vector<double> dd(static_cast<std::size_t>(n), 0.0);
    compute_cubic_2nd_deriv(cubic_2nd_deriv_init_e::Cubic2ndDerivFirst, n, dd, d0, dn);
    set_quadratic_poly(_poly[0], _data_set->get_x_value(0), _data_set->get_y_value(0), d0,
                       dd[0]);
    for (int i = 1; i < n; i++)
        set_cubic_poly2(_poly[static_cast<std::size_t>(i)], _data_set->get_x_value(i - 1),
                        _data_set->get_y_value(i - 1), _data_set->get_x_value(i),
                        _data_set->get_y_value(i), dd[static_cast<std::size_t>(i - 1)],
                        dd[static_cast<std::size_t>(i)]);
    set_quadratic_poly(_poly[static_cast<std::size_t>(n)], _data_set->get_x_value(n - 1),
                       _data_set->get_y_value(n - 1), dn,
                       dd[static_cast<std::size_t>(n - 1)]);
    _invalid = false;
    return interpolate_cubic(d, x);
}

double Interpolated1d::update_cubic_deriv_init(int d, double x) {
    int n = _data_set->get_count();
    if (n < 4)
        throw IllegalStateException("data set doesn't contains enough data");
    resizePoly(n + 1);
    std::vector<double> dd(static_cast<std::size_t>(n), 0.0);
    double d0 = _data_set->get_d_value(0);
    double dn = _data_set->get_d_value(n - 1);
    compute_cubic_2nd_deriv(cubic_2nd_deriv_init_e::Cubic2ndDerivFirst,
                            _data_set->get_count(), dd, d0, dn);
    set_linear_poly(_poly[0], _data_set->get_x_value(0), _data_set->get_y_value(0), d0);
    for (int i = 1; i < n; i++)
        set_cubic_poly2(_poly[static_cast<std::size_t>(i)], _data_set->get_x_value(i - 1),
                        _data_set->get_y_value(i - 1), _data_set->get_x_value(i),
                        _data_set->get_y_value(i), dd[static_cast<std::size_t>(i - 1)],
                        dd[static_cast<std::size_t>(i)]);
    set_linear_poly(_poly[static_cast<std::size_t>(n)], _data_set->get_x_value(n - 1),
                    _data_set->get_y_value(n - 1), dn);
    _invalid = false;
    return interpolate_cubic(d, x);
}

double Interpolated1d::update_cubic2_deriv_init(int d, double x) {
    int n = _data_set->get_count();
    if (n < 4)
        throw IllegalStateException("data set doesn't contains enough data");
    resizePoly(n + 1);
    std::vector<double> dd(static_cast<std::size_t>(n), 0.0);
    double d0 = _data_set->get_d_value(0);
    double dn = _data_set->get_d_value(n - 1);
    compute_cubic_2nd_deriv(cubic_2nd_deriv_init_e::Cubic2ndDerivFirst,
                            _data_set->get_count(), dd, d0, dn);
    set_quadratic_poly(_poly[0], _data_set->get_x_value(0), _data_set->get_y_value(0), d0,
                       dd[0]);
    for (int i = 1; i < n; i++)
        set_cubic_poly2(_poly[static_cast<std::size_t>(i)], _data_set->get_x_value(i - 1),
                        _data_set->get_y_value(i - 1), _data_set->get_x_value(i),
                        _data_set->get_y_value(i), dd[static_cast<std::size_t>(i - 1)],
                        dd[static_cast<std::size_t>(i)]);
    set_quadratic_poly(_poly[static_cast<std::size_t>(n)], _data_set->get_x_value(n - 1),
                       _data_set->get_y_value(n - 1), dn,
                       dd[static_cast<std::size_t>(n - 1)]);
    _invalid = false;
    return interpolate_cubic(d, x);
}

double Interpolated1d::update_cubic2_deriv(int d, double x) {
    int n = _data_set->get_count();
    if (n < 4)
        throw IllegalStateException("data set doesn't contains enough data");
    resizePoly(n + 1);
    double dd0 =
        (_data_set->get_d_value(1) - _data_set->get_d_value(0)) / _data_set->get_x_interval(0);
    set_quadratic_poly(_poly[0], _data_set->get_x_value(0), _data_set->get_y_value(0),
                       _data_set->get_d_value(0), dd0);
    for (int i = 1; i < n; i++)
        set_cubic_poly(_poly[static_cast<std::size_t>(i)], _data_set->get_x_value(i - 1),
                       _data_set->get_y_value(i - 1), _data_set->get_x_value(i),
                       _data_set->get_y_value(i), _data_set->get_d_value(i - 1),
                       _data_set->get_d_value(i));
    double ddn = (_data_set->get_d_value(n - 1) - _data_set->get_d_value(n - 2)) /
                 _data_set->get_x_interval(n - 2);
    set_quadratic_poly(_poly[static_cast<std::size_t>(n)], _data_set->get_x_value(n - 1),
                       _data_set->get_y_value(n - 1), _data_set->get_d_value(n - 1), ddn);
    _invalid = false;
    return interpolate_cubic(d, x);
}

double Interpolated1d::update_cubic_deriv(int d, double x) {
    int n = _data_set->get_count();
    if (n < 4)
        throw IllegalStateException("data set doesn't contains enough data");
    resizePoly(n + 1);
    set_linear_poly(_poly[0], _data_set->get_x_value(0), _data_set->get_y_value(0),
                    _data_set->get_d_value(0));
    for (int i = 1; i < n; i++)
        set_cubic_poly(_poly[static_cast<std::size_t>(i)], _data_set->get_x_value(i - 1),
                       _data_set->get_y_value(i - 1), _data_set->get_x_value(i),
                       _data_set->get_y_value(i), _data_set->get_d_value(i - 1),
                       _data_set->get_d_value(i));
    set_linear_poly(_poly[static_cast<std::size_t>(n)], _data_set->get_x_value(n - 1),
                    _data_set->get_y_value(n - 1), _data_set->get_d_value(n - 1));
    _invalid = false;
    return interpolate_cubic(d, x);
}

} // namespace redukti::data
