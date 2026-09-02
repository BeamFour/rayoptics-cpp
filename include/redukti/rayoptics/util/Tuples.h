// C++ port of org.redukti.rayoptics.util.{Pair,Triple,Quad,Quint}
//
// Four one-class Java files collapsed into one header: they are the same
// generic container at four arities and nothing else in the port depends on
// their file layout.
#ifndef REDUKTI_RAYOPTICS_UTIL_TUPLES_H
#define REDUKTI_RAYOPTICS_UTIL_TUPLES_H

namespace redukti::rayoptics::util {

template <typename T1, typename T2> struct Pair {
    T1 first;
    T2 second;

    Pair() = default;
    Pair(T1 first_, T2 second_) : first(first_), second(second_) {}

    bool operator==(const Pair &o) const {
        return first == o.first && second == o.second;
    }
};

template <typename T1, typename T2, typename T3> struct Triple {
    T1 first;
    T2 second;
    T3 third;

    Triple() = default;
    Triple(T1 first_, T2 second_, T3 third_)
        : first(first_), second(second_), third(third_) {}
};

template <typename T1, typename T2, typename T3, typename T4> struct Quad {
    T1 first;
    T2 second;
    T3 third;
    T4 fourth;

    Quad() = default;
    Quad(T1 first_, T2 second_, T3 third_, T4 fourth_)
        : first(first_), second(second_), third(third_), fourth(fourth_) {}
};

template <typename T1, typename T2, typename T3, typename T4, typename T5> struct Quint {
    T1 first;
    T2 second;
    T3 third;
    T4 fourth;
    T5 fifth;

    Quint() = default;
    Quint(T1 first_, T2 second_, T3 third_, T4 fourth_, T5 fifth_)
        : first(first_), second(second_), third(third_), fourth(fourth_), fifth(fifth_) {}
};

} // namespace redukti::rayoptics::util

#endif // REDUKTI_RAYOPTICS_UTIL_TUPLES_H
