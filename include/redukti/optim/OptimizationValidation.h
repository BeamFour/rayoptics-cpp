// C++ port of org.redukti.optim.OptimizationValidation
#ifndef REDUKTI_OPTIM_OPTIMIZATIONVALIDATION_H
#define REDUKTI_OPTIM_OPTIMIZATIONVALIDATION_H

#include <algorithm>
#include <cmath>
#include <set>
#include <vector>

namespace redukti::optim {

/**
 * Domain rules shared by the Java API and text reader; callers supply contextual errors.
 *
 * The Java takes a Supplier of the exception; here `error` is any callable returning the
 * exception to throw, so each caller keeps its own message and exception type.
 */
namespace OptimizationValidation {

template <typename Error>
void fieldCount(const std::vector<double> &values, int count, Error error) {
    if (static_cast<int>(values.size()) != count)
        throw error();
}

/** Inclusive non-negative range; an infinite upper limit still rejects non-finite values. */
template <typename Error>
void range(const std::vector<double> &values, double upper, Error error) {
    for (double value : values)
        if (!std::isfinite(value) || value < 0.0 || value > upper)
            throw error();
}

template <typename Error>
void frequency(int value, const std::vector<int> &measured, Error error) {
    if (std::find(measured.begin(), measured.end(), value) == measured.end())
        throw error();
}

template <typename Error>
void positiveUnique(int value, std::set<int> &seen, Error error) {
    if (value <= 0 || !seen.insert(value).second)
        throw error();
}

} // namespace OptimizationValidation

} // namespace redukti::optim

#endif // REDUKTI_OPTIM_OPTIMIZATIONVALIDATION_H
