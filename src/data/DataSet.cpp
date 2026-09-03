// C++ port of org.redukti.data.DataSet and DiscreteSetBase
#include "redukti/data/DataSet.h"

#include "redukti/Exceptions.h"

#include <limits>

namespace redukti::data {

Range DataSet::get_y_range() const {
    // Java seeds with (Double.MAX_VALUE, Double.MIN_VALUE). Double.MIN_VALUE is
    // the smallest *positive* value, not the most negative, so the upper seed is
    // 4.9e-324 -- denorm_min, not numeric_limits<double>::min(), which is the
    // smallest positive *normal*. The Java behaviour is reproduced as written.
    Range r(std::numeric_limits<double>::max(),
            std::numeric_limits<double>::denorm_min());
    int d = get_dimensions();
    std::vector<int> x(static_cast<std::size_t>(d));
    std::vector<int> c(static_cast<std::size_t>(d));
    for (int i = 0; i < d; i++) {
        if (get_count(i) == 0)
            throw IllegalStateException("data set contains no data");
        x[static_cast<std::size_t>(i)] = 0;
        c[static_cast<std::size_t>(i)] = get_count(i) - 1;
    }
    while (true) {
        double y = get_y_value(x);
        if (y < r.first)
            r.first = y;
        if (y > r.second)
            r.second = y;
        for (int i = 0;;) {
            if (x[static_cast<std::size_t>(i)] < c[static_cast<std::size_t>(i)]) {
                x[static_cast<std::size_t>(i)]++;
                break;
            } else {
                x[static_cast<std::size_t>(i++)] = 0;
                if (i == d)
                    return r;
            }
        }
    }
}

void DiscreteSetBase::add_data(double x, double y, double d) {
    EntryS e(x, y, d);
    _version++;
    int di = get_interval(x);
    if (di > 0 && (_data[static_cast<std::size_t>(di - 1)].x == x))
        _data[static_cast<std::size_t>(di - 1)] = e;
    else
        _data.insert(_data.begin() + di, e);
    invalidate();
}

void DiscreteSetBase::clear() {
    _data.clear();
    _version++;
    invalidate();
}

Range DiscreteSetBase::get_x_range() const {
    if (_data.empty())
        throw IllegalStateException("_data set contains no _data");
    return Range(_data[0].x, _data[_data.size() - 1].x);
}

int DiscreteSetBase::get_interval(double x) const {
    int min_idx = 0;
    int max_idx = static_cast<int>(_data.size()) + 1;
    while (max_idx - min_idx > 1) {
        int p = (max_idx + min_idx) / 2;
        if (x >= _data[static_cast<std::size_t>(p - 1)].x)
            min_idx = p;
        else
            max_idx = p;
    }
    return min_idx;
}

int DiscreteSetBase::get_nearest(double x) const {
    int min_idx = 0;
    int max_idx = static_cast<int>(_data.size());
    while (max_idx - min_idx > 1) {
        int p = (max_idx + min_idx) / 2;
        if (x + x >= _data[static_cast<std::size_t>(p - 1)].x +
                         _data[static_cast<std::size_t>(p)].x)
            min_idx = p;
        else
            max_idx = p;
    }
    return min_idx;
}

double DiscreteSetBase::get_x_interval(int x) const {
    return _data[static_cast<std::size_t>(x + 1)].x - _data[static_cast<std::size_t>(x)].x;
}

double DiscreteSetBase::get_x_interval(int x1, int x2) const {
    return _data[static_cast<std::size_t>(x2)].x - _data[static_cast<std::size_t>(x1)].x;
}

} // namespace redukti::data
