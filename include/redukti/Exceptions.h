// C++ port of the Java rayoptics module.
// Original software https://github.com/mjhoptics/ray-optics (Michael J. Hayford)
// Java version by Dibyendu Majumdar
#ifndef REDUKTI_EXCEPTIONS_H
#define REDUKTI_EXCEPTIONS_H

#include <exception>
#include <string>
#include <utility>

namespace redukti {

/**
 * The Java exception hierarchy is reproduced here so that the ported catch
 * clauses keep their original meaning. In particular the ray tracer relies on
 * catching a common base (`RuntimeException`, and in the optimizer everything
 * derived from `Exception`) and mapping the failure onto a penalty value.
 */
class Exception : public std::exception {
public:
    Exception() = default;
    explicit Exception(std::string message) : message_(std::move(message)) {}
    const char *what() const noexcept override { return message_.c_str(); }
    const std::string &getMessage() const noexcept { return message_; }

private:
    std::string message_;
};

class RuntimeException : public Exception {
public:
    RuntimeException() = default;
    explicit RuntimeException(std::string message) : Exception(std::move(message)) {}
};

class IllegalArgumentException : public RuntimeException {
public:
    IllegalArgumentException() = default;
    explicit IllegalArgumentException(std::string message) : RuntimeException(std::move(message)) {}
};

class IllegalStateException : public RuntimeException {
public:
    IllegalStateException() = default;
    explicit IllegalStateException(std::string message) : RuntimeException(std::move(message)) {}
};

class NoSuchElementException : public RuntimeException {
public:
    NoSuchElementException() = default;
    explicit NoSuchElementException(std::string message) : RuntimeException(std::move(message)) {}
};

class IndexOutOfBoundsException : public RuntimeException {
public:
    IndexOutOfBoundsException() = default;
    explicit IndexOutOfBoundsException(std::string message) : RuntimeException(std::move(message)) {}
};

class NumberFormatException : public IllegalArgumentException {
public:
    NumberFormatException() = default;
    explicit NumberFormatException(std::string message) : IllegalArgumentException(std::move(message)) {}
};

class UnsupportedOperationException : public RuntimeException {
public:
    UnsupportedOperationException() = default;
    explicit UnsupportedOperationException(std::string message) : RuntimeException(std::move(message)) {}
};

} // namespace redukti

#endif // REDUKTI_EXCEPTIONS_H
