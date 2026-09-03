// C++ port of org.redukti.rayoptics.util.Lists
#ifndef REDUKTI_RAYOPTICS_UTIL_LISTS_H
#define REDUKTI_RAYOPTICS_UTIL_LISTS_H

#include "redukti/Exceptions.h"
#include "redukti/rayoptics/util/Tuples.h"

#include <algorithm>
#include <memory>
#include <optional>
#include <type_traits>
#include <vector>

namespace redukti::rayoptics::util {

/**
 * Python-style list slicing, carried over from the upstream ray-optics code.
 * Java passes boxed Integers so that null means "unbounded"; those become
 * std::optional<int> here.
 */
namespace Lists {

template <typename E>
std::vector<E> slice(const std::vector<E> &inputList, std::optional<int> start_,
                     std::optional<int> stop_, std::optional<int> step_) {
    std::vector<E> newList;
    // A slice can never contain more elements than its source. Reserving the
    // upper bound avoids repeatedly copying value types and incrementing
    // shared_ptr reference counts while the vector grows.
    newList.reserve(inputList.size());
    int step = step_.has_value() ? *step_ : 1;
    int length = static_cast<int>(inputList.size());
    int start, stop;
    if (!start_.has_value()) {
        start = (step < 0) ? length - 1 : 0;
    } else {
        start = (*start_ < 0) ? *start_ + length : *start_;
    }
    if (!stop_.has_value()) {
        stop = (step < 0) ? -1 : length;
    } else {
        stop = (*stop_ < 0) ? *stop_ + length : *stop_;
    }
    if (step < 0) {
        for (int i = std::min(start, length - 1); i >= std::max(stop, 0); i += step) {
            newList.push_back(inputList[static_cast<std::size_t>(i)]);
        }
    } else if (step > 0) {
        for (int i = std::max(0, start); i < std::min(stop, length); i += step) {
            newList.push_back(inputList[static_cast<std::size_t>(i)]);
        }
    } else {
        throw IllegalArgumentException();
    }
    return newList;
}

template <typename E> std::vector<E> from(const std::vector<E> &inputList, int start) {
    return slice(inputList, start, std::nullopt, std::nullopt);
}

template <typename E> std::vector<E> upto(const std::vector<E> &inputList, int stop) {
    return slice(inputList, std::nullopt, stop, std::nullopt);
}

template <typename E> std::vector<E> step(const std::vector<E> &inputList, int step_) {
    return slice(inputList, std::nullopt, std::nullopt, step_);
}

/** Negative indices count back from the end, as in Python. */
template <typename E> E &get(std::vector<E> &inputList, int i) {
    if (i < 0)
        i += static_cast<int>(inputList.size());
    return inputList[static_cast<std::size_t>(i)];
}

template <typename E> const E &get(const std::vector<E> &inputList, int i) {
    if (i < 0)
        i += static_cast<int>(inputList.size());
    return inputList[static_cast<std::size_t>(i)];
}

template <typename E> void set(std::vector<E> &inputList, int i, const E &e) {
    if (i < 0)
        i += static_cast<int>(inputList.size());
    inputList[static_cast<std::size_t>(i)] = e;
}

/**
 * Only types with a genuine "absent" state may be zipped. Java pads the short
 * list with null and callers test for it -- SequentialModel's path walk does
 * `if (g != null)` on the gap, because `ifcs` always has one more entry than
 * `gaps` (the image surface has no gap after it). Value-initialising a value
 * type would produce a default-constructed object that reads as present, which
 * compiles and is silently wrong, so those instantiations are rejected here
 * instead.
 */
template <typename T> struct is_nullable : std::false_type {};
template <typename T> struct is_nullable<T *> : std::true_type {};
template <typename T> struct is_nullable<std::shared_ptr<T>> : std::true_type {};
template <typename T> struct is_nullable<std::optional<T>> : std::true_type {};

/** Java pads the shorter list with null; here that is a null pointer or nullopt. */
template <typename T1, typename T2>
std::vector<Pair<T1, T2>> zip_longest(const std::vector<T1> &list1,
                                      const std::vector<T2> &list2) {
    static_assert(is_nullable<T1>::value && is_nullable<T2>::value,
                  "zip_longest pads the shorter list with null; instantiate it with "
                  "pointer, shared_ptr or optional elements so the padding is "
                  "distinguishable from a real value");
    std::vector<Pair<T1, T2>> list;
    std::size_t n = std::max(list1.size(), list2.size());
    for (std::size_t i = 0; i < n; i++) {
        T1 t1 = i < list1.size() ? list1[i] : T1{};
        T2 t2 = i < list2.size() ? list2[i] : T2{};
        list.push_back(Pair<T1, T2>(t1, t2));
    }
    return list;
}

} // namespace Lists
} // namespace redukti::rayoptics::util

#endif // REDUKTI_RAYOPTICS_UTIL_LISTS_H
