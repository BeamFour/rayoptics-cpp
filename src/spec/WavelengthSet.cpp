// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.spec.WavelengthSet
#include "redukti/spec/WavelengthSet.h"

#include "redukti/Exceptions.h"
#include "redukti/Text.h"
#include "redukti/rayoptics/specs/WvlSpec.h"
#include "redukti/spec/Prescription.h"

#include <algorithm>
#include <cmath>
#include <cstring>

namespace redukti::spec {

using rayoptics::specs::WvlSpec;

namespace {

/**
 * The spectral lines a set can be spelled with, longest wavelength first. The
 * wavelengths themselves come from WvlSpec::get_wavelength, which is
 * the ported table and the one authority for what a line's wavelength is; this list
 * only fixes the order a set is written in and the names a trial may use. Laser lines
 * `Nd` and `He-Ne` are resolvable there too but are not Fraunhofer lines,
 * so they are left out of the names a trial writes.
 *
 * Wider than what a catalog glass can answer, deliberately; see the class comment.
 */
const std::vector<std::string> &LINES() {
    static const std::vector<std::string> lines = {"t", "s",  "r", "C", "C'", "D", "d",
                                                   "e", "F", "F'", "g", "h",  "i"};
    return lines;
}

/** Within this of a line's wavelength is that line, for writing a set back out. */
constexpr double LINE_TOLERANCE = 1.0e-6;

/** Java's String.join(" ", ...). */
std::string joinWithSpaces(const std::vector<std::string> &values) {
    std::string result;
    for (const std::string &value : values) {
        if (!result.empty())
            result += " ";
        result += value;
    }
    return result;
}

/** Java's Arrays.toString(double[]). */
std::string arrayToString(const std::vector<double> &v) {
    std::string sb = "[";
    for (std::size_t i = 0; i < v.size(); i++) {
        if (i > 0)
            sb += ", ";
        sb += doubleToString(v[i]);
    }
    sb += "]";
    return sb;
}

/** Java's Arrays.equals(double[], double[]), which compares bit patterns. */
bool arraysEqual(const std::vector<double> &a, const std::vector<double> &b) {
    return a.size() == b.size() &&
           (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
}

} // namespace

WavelengthSet WavelengthSet::dLine() {
    return ofLetters({"d"});
}

WavelengthSet WavelengthSet::eLine() {
    return ofLetters({"e"});
}

WavelengthSet WavelengthSet::standardTriColor() {
    return ofLetters({"d", "F", "C"});
}

WavelengthSet WavelengthSet::standardFiveColor() {
    return ofLetters({"d", "C", "e", "F", "g"},
                     std::vector<double>{WT_d, WT_C, WT_e, WT_F, WT_g});
}

WavelengthSet WavelengthSet::k9Tester() {
    return ofLetters({"d", "C", "e", "F", "g"},
                     std::vector<double>{K9_WT_d, K9_WT_C, K9_WT_e, K9_WT_F, K9_WT_g});
}

WavelengthSet WavelengthSet::of(const std::vector<double> &wavelengths,
                                const std::vector<double> &weights) {
    if (wavelengths.empty())
        throw IllegalArgumentException("at least one wavelength is required");
    if (weights.size() != wavelengths.size())
        throw IllegalArgumentException(
            "expected one weight per wavelength: " +
            intToString(static_cast<int>(wavelengths.size())) + " wavelengths but " +
            intToString(static_cast<int>(weights.size())) + " weights");
    for (double wavelength : wavelengths)
        if (!std::isfinite(wavelength) || wavelength <= 0.0)
            throw IllegalArgumentException("wavelengths must be positive, in nanometres");
    for (double weight : weights)
        if (!std::isfinite(weight) || weight < 0.0)
            throw IllegalArgumentException(
                "wavelength weights must be finite and non-negative");
    return WavelengthSet(wavelengths, weights);
}

WavelengthSet WavelengthSet::of(const Prescription &prescription) {
    return of(prescription._wvls, prescription._wts);
}

WavelengthSet WavelengthSet::ofLetters(const std::vector<std::string> &letters,
                                       const std::optional<std::vector<double>> &weights) {
    if (letters.empty())
        throw IllegalArgumentException("at least one wavelength is required");
    std::vector<double> wavelengths(letters.size());
    for (std::size_t i = 0; i < wavelengths.size(); i++)
        wavelengths[i] = wavelengthOf(letters[i]);
    if (!weights.has_value())
        return of(wavelengths, std::vector<double>(wavelengths.size(), 1.0));
    return of(wavelengths, *weights);
}

double WavelengthSet::wavelengthOf(const std::string &letter) {
    const auto &lines = LINES();
    if (std::find(lines.begin(), lines.end(), letter) == lines.end())
        throw IllegalArgumentException("unknown wavelength '" + letter +
                                       "'; expected one of: " + joinWithSpaces(lines));
    return WvlSpec::get_wavelength(letter);
}

std::optional<std::string> WavelengthSet::letterOf(double wavelength) {
    for (const std::string &line : LINES())
        if (std::abs(WvlSpec::get_wavelength(line) - wavelength) <= LINE_TOLERANCE)
            return line;
    return std::nullopt;
}

bool WavelengthSet::hasEqualWeights() const {
    for (double weight : _weights)
        if (weight != 1.0)
            return false;
    return true;
}

bool WavelengthSet::hasStandardWeights() const {
    return *this == standardFiveColor();
}

std::string WavelengthSet::letters() const {
    std::vector<std::string> result;
    result.reserve(_wavelengths.size());
    for (double wavelength : _wavelengths) {
        std::optional<std::string> letter = letterOf(wavelength);
        if (!letter.has_value())
            throw IllegalStateException(
                "wavelength " + doubleToString(wavelength) +
                " nm is not a Fraunhofer line, so it cannot be written as a trial");
        result.push_back(*letter);
    }
    return joinWithSpaces(result);
}

bool WavelengthSet::operator==(const WavelengthSet &that) const {
    return arraysEqual(_wavelengths, that._wavelengths) &&
           arraysEqual(_weights, that._weights);
}

std::string WavelengthSet::toString() const {
    return "WavelengthSet{" + arrayToString(_wavelengths) + ", weights=" +
           arrayToString(_weights) + "}";
}

} // namespace redukti::spec
