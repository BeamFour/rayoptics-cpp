// C++ port of the org.redukti.data set hierarchy:
//   Range, Interpolation, InterpolatableDataSet, DataSet, Set1d,
//   DiscreteSetBase and DiscreteSet.
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_DATA_DATASET_H
#define REDUKTI_DATA_DATASET_H

#include "redukti/data/Interpolated1d.h"

#include <vector>

namespace redukti::data {

/** Java's `class Range`: a mutable (first, second) pair of doubles. */
class Range {
public:
    double first;
    double second;

    Range(double first_, double second_) : first(first_), second(second_) {}
};

/** Java's abstract `DataSet`, the n-dimensional base. */
class DataSet {
public:
    virtual ~DataSet() = default;

    virtual int get_dimensions() const = 0;
    virtual int get_count(int dim) const = 0;
    virtual double get_x_value(int n, int dim) const = 0;
    virtual double get_y_value(const std::vector<int> &x) const = 0;
    virtual double interpolate(const std::vector<double> &x) = 0;
    virtual double interpolate(const std::vector<double> &x, int deriv, int dim) = 0;
    virtual Range get_x_range(int dim) const = 0;

    Range get_y_range() const;

    int get_version() const { return _version; }

protected:
    int _version = 0;
    Interpolation _interpolation = Interpolation::Linear;
};

/** Java's abstract `Set1d`: a DataSet that is one-dimensional. */
class Set1d : public DataSet {
public:
    virtual int get_count() const = 0;
    virtual double get_x_value(int n) const = 0;
    virtual double get_y_value(int n) const = 0;
    virtual double interpolate(double x) = 0;
    virtual double interpolate(double x, int deriv) = 0;
    virtual Range get_x_range() const = 0;

    int get_dimensions() const override { return 1; }

    int get_count(int dimension) const override {
        (void)dimension; // Java asserts dimension == 0
        return get_count();
    }

    double get_x_value(int x, int dimension) const override {
        (void)dimension;
        return get_x_value(x);
    }

    double get_y_value(const std::vector<int> &x) const override {
        return get_y_value(x[0]);
    }

    Range get_x_range(int dimension) const override {
        (void)dimension;
        return get_x_range();
    }

    double interpolate(const std::vector<double> &x) override { return interpolate(x[0]); }

    double interpolate(const std::vector<double> &x, int deriv, int dimension) override {
        (void)dimension;
        return interpolate(x[0], deriv);
    }
};

/** Java's abstract `DiscreteSetBase`: sorted (x, y, derivative) samples. */
class DiscreteSetBase : public Set1d, public InterpolatableDataSet {
public:
    class EntryS {
    public:
        double x, y, d;

        EntryS(double x_, double y_, double d_) : x(x_), y(y_), d(d_) {}
    };

    void add_data(double x, double y, double d);
    void add_data(double x, double y) { add_data(x, y, 0.0); }

    void clear();

    // ---- InterpolatableDataSet ----
    double get_d_value(int n) const override { return _data[static_cast<std::size_t>(n)].d; }
    int get_count() const override { return static_cast<int>(_data.size()); }
    double get_x_value(int n) const override { return _data[static_cast<std::size_t>(n)].x; }
    double get_y_value(int n) const override { return _data[static_cast<std::size_t>(n)].y; }
    int get_interval(double x) const override;
    int get_nearest(double x) const override;
    double get_x_interval(int x) const override;
    double get_x_interval(int x1, int x2) const override;

    Range get_x_range() const override;

    // Set1d also declares get_x_value(int)/get_y_value(int)/get_count(); the
    // overrides above satisfy both bases, but the name lookup needs help.
    using Set1d::get_x_value;
    using Set1d::get_y_value;
    using Set1d::get_count;

protected:
    virtual void invalidate() = 0;

    std::vector<EntryS> _data;
};

/** Java's `DiscreteSet`: a DiscreteSetBase with an Interpolated1d attached. */
class DiscreteSet : public DiscreteSetBase {
public:
    DiscreteSet() : _interpolated_1d(this) {}

    /**
     * Copy and move have to re-point the interpolator: it holds a bare pointer
     * to the DiscreteSet that owns it, so the compiler-generated versions would
     * leave the new object interpolating over the old one. Plot borrows its
     * data sets, so they do get stored and passed around.
     */
    DiscreteSet(const DiscreteSet &o)
        : DiscreteSetBase(o), _interpolated_1d(o._interpolated_1d) {
        _interpolated_1d.rebind(this);
    }

    DiscreteSet(DiscreteSet &&o) noexcept
        : DiscreteSetBase(std::move(o)), _interpolated_1d(std::move(o._interpolated_1d)) {
        _interpolated_1d.rebind(this);
    }

    DiscreteSet &operator=(const DiscreteSet &o) {
        if (this != &o) {
            DiscreteSetBase::operator=(o);
            _interpolated_1d = o._interpolated_1d;
            _interpolated_1d.rebind(this);
        }
        return *this;
    }

    DiscreteSet &operator=(DiscreteSet &&o) noexcept {
        if (this != &o) {
            DiscreteSetBase::operator=(std::move(o));
            _interpolated_1d = std::move(o._interpolated_1d);
            _interpolated_1d.rebind(this);
        }
        return *this;
    }

    double interpolate(double x) override { return _interpolated_1d.interpolate(x); }

    double interpolate(double x, int deriv) override {
        return _interpolated_1d.interpolate(x, deriv);
    }

    void set_interpolation(Interpolation i) { _interpolated_1d.set_interpolation(i); }

    // Bring the vector-taking overloads from Set1d back into scope.
    using Set1d::interpolate;

protected:
    void invalidate() override { _interpolated_1d.invalidate(); }

    Interpolated1d _interpolated_1d;
};

} // namespace redukti::data

#endif // REDUKTI_DATA_DATASET_H
