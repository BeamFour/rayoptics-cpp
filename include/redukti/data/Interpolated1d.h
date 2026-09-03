// C++ port of org.redukti.data.Interpolation, InterpolatableDataSet and
// Interpolated1d.
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_DATA_INTERPOLATED1D_H
#define REDUKTI_DATA_INTERPOLATED1D_H

#include <vector>

namespace redukti::data {

enum class Interpolation {
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
    virtual double get_x_value(int n) const = 0;
    virtual double get_y_value(int n) const = 0;
    virtual double get_d_value(int n) const = 0;
    virtual int get_nearest(double x) const = 0;
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
