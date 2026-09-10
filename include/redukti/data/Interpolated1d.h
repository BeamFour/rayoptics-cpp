// The software is ported from Goptical, hence is licensed under the GPL.
// Copyright (c) 2021 Dibyendu Majumdar
// Goptical: Copyright (C) 2010-2011 Free Software Foundation, Inc; Author: Alexandre Becoulet
// Licensed under the GNU General Public License, version 3 or later; see LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.data.Interpolation, InterpolatableDataSet and
// Interpolated1d.
#ifndef REDUKTI_DATA_INTERPOLATED1D_H
#define REDUKTI_DATA_INTERPOLATED1D_H

#include <vector>

namespace redukti::data {

/** Specifies data interpolation methods. Availability depends on data
 * container used. */
enum class Interpolation {
    /** 1d and 2d nearest interpolation */
    /** 1d linear and 2d bilinear interpolations */
    /** 1d quadratic interpolation */
    /** 1d cubic piecewise interpolation. It has continuous
     piecewise first derivative, non-continuous piecewise
     linear second derivative. Use segments slope as first
     derivative. Less accurate than other cubic interpolants
     but requires less computation on data set change. */
    /** 1d cubic piecewise interpolation. It has smooth first
     derivative and continuous piecewise linear second
     derivative. Derivatives for first and last entries are
     based on first and last segments slope. It uses linear
     extrapolation (continuous but non-smooth first derivative
     on both ends). */
    /** Same interpolation as Cubic, with quadratic extrapolation
     (continous and smooth first derivative on both ends). */
    /** Same as Cubic with first derivative prescribed for first
     and last entries only. */
    /** Same as Cubic2 with first derivative prescribed for first
     and last entries only. */
    /** 1d cubic piecewise interpolation. First derivatives must
     be provided for all entries. It uses linear extrapolation.*/
    /** 1d cubic piecewise interpolation. First derivatives must
     be provided for all entries. It uses quadratic extrapolation.*/
    /** 2d bicubic interpolation. Use smooth first derivative and
     continuous piecewise linear second derivative. Use 1d
     cubic curve to extract gradients (smooth first derivative
     and continuous piecewise linear second derivative). This
     is the best 2d interpolation when derivatives are
     non-prescribed. */
    /** 2d bicubic interpolation. Use numerical differencing to
     extract gradients. Less accurate than @ref Bicubic but
     requires less computation on data set change.*/
    /** 2d bicubic interpolation. x and y gradients must be
     provided. This is the best 2d interpolation when
     derivatives values are available. */
    Nearest,
    Linear,
    Quadratic,
    CubicSimple,
    Cubic,
    Cubic2,
    CubicDerivInit,
    Cubic2DerivInit,
    CubicDeriv,
    Cubic2Deriv,
    Bicubic,
    BicubicDiff,
    BicubicDeriv,
};

/** Java's `interface InterpolatableDataSet`. */
class InterpolatableDataSet {
public:
    virtual ~InterpolatableDataSet() = default;

    virtual double get_x_interval(int x) const = 0;
    virtual double get_x_interval(int x1, int x2) const = 0;
    virtual int get_interval(double x) const = 0;
    /** Get x data at index n in data set */
    virtual double get_x_value(int n) const = 0;
    /** Get y data stored at index n in data set */
    virtual double get_y_value(int n) const = 0;
    virtual double get_d_value(int n) const = 0;
    virtual int get_nearest(double x) const = 0;
    /** Get total number of data stored in data set */
    virtual int get_count() const = 0;
};

/**
 * Piecewise polynomial interpolation over an InterpolatableDataSet.
 *
 * The polynomials are rebuilt lazily: `invalidate()` marks them stale and the
 * next `interpolate()` recomputes, which is why each update_* method ends by
 * evaluating the point that triggered it.
 */
class Interpolated1d {
public:
    /** Java's nested `enum cubic_2nd_deriv_init_e`. */
    enum class cubic_2nd_deriv_init_e {
        Cubic2ndDerivQuadratic,
        Cubic2ndDerivFirst,
        Cubic2ndDerivSecond,
    };

    class PolyS {
    public:
        double a, b, c, d;

        PolyS(double a_, double b_, double c_, double d_) : a(a_), b(b_), c(c_), d(d_) {}
    };

    /** Borrowed; the data set owns this object, so it outlives it. */
    explicit Interpolated1d(InterpolatableDataSet *dataSet) : _data_set(dataSet) {}

    void invalidate() { _invalid = true; }

    /**
     * Re-point at the owning data set after that set has been copied or moved.
     * The constructor is handed `this` of the enclosing DiscreteSet, so a copy
     * would otherwise keep interpolating over the original.
     */
    void rebind(InterpolatableDataSet *dataSet) { _data_set = dataSet; }

    void set_interpolation(Interpolation i);

    double interpolate(double x) { return interpolate(x, 0); }

    double interpolate(double x, int d);

protected:
    InterpolatableDataSet *_data_set;
    std::vector<PolyS> _poly;
    bool _invalid = true;
    Interpolation _method = Interpolation::Linear;

private:
    void resizePoly(int n);

    void compute_cubic_2nd_deriv(cubic_2nd_deriv_init_e de, int n, std::vector<double> &dd,
                                 double d0, double dn);

    static void set_linear_poly(PolyS &p, double p1x, double p1y, double p2x, double p2y);
    static void set_linear_poly(PolyS &p, double p1x, double p1y, double d1);
    static void set_quadratic_poly(PolyS &p, double p1x, double p1y, double p2x,
                                   double p2y, double p3x, double p3y);
    static void set_quadratic_poly(PolyS &p, double px, double py, double d, double dd);
    static void set_cubic_poly(PolyS &p, double p1x, double p1y, double p2x, double p2y,
                               double d1, double d2);
    static void set_cubic_poly2(PolyS &p, double p1x, double p1y, double p2x, double p2y,
                                double dd1, double dd2);

    double interpolate_nearest(int d, double x);
    double interpolate_linear(int d, double x);
    double interpolate_quadratic(int d, double x);
    double interpolate_cubic(int d, double x);

    double update_quadratic(int d, double x);
    double update_cubic_simple(int d, double x);
    double update_cubic(int d, double x);
    double update_cubic2(int d, double x);
    double update_cubic_deriv_init(int d, double x);
    double update_cubic2_deriv_init(int d, double x);
    double update_cubic_deriv(int d, double x);
    double update_cubic2_deriv(int d, double x);
};

} // namespace redukti::data

#endif // REDUKTI_DATA_INTERPOLATED1D_H
