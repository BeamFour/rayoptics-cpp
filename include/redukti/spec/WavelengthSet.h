// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.spec.WavelengthSet
#ifndef REDUKTI_SPEC_WAVELENGTHSET_H
#define REDUKTI_SPEC_WAVELENGTHSET_H

#include <optional>
#include <string>
#include <vector>

namespace redukti::spec {

class Prescription;

/**
 * The wavelengths a design is evaluated at, and the weight each one carries.
 *
 * The weights do two jobs, both of them the same decision: the analyses combine
 * per-wavelength results with them, and the optimizer multiplies them into the weight of
 * every per-wavelength residual. A wavelength weighted 0.15 therefore counts for 0.15 of a
 * wavelength weighted 1.0 in a polychromatic MTF and in the merit alike.
 *
 * Three sets cover the usual choices - dLine(), standardTriColor() and
 * standardFiveColor() - with k9Tester() for comparing against measured
 * curves. Nothing here is limited to them: of() takes any wavelengths and weights.
 *
 * <b>What the glasses can answer limits the choice.</b> A glass from a catalog holds
 * measured indices at d, C, F, e, C', F' and g only, and
 * rayoptics::seq::Glass::rindex(double) throws for any other wavelength.
 * A glass given as nd and vd computes its index from a dispersion model instead, so any
 * wavelength works there as an approximation. The other lines ofLetters() accepts
 * are for that case, and for when more wavelength data can be loaded.
 */
class WavelengthSet {
public:
    /**
     * The weights of standardFiveColor(), in its order. They are the relative
     * wavelength weights of the MTF design data in T. Steinich and V. Blahnik, "Optical
     * design of camera optics for mobile phones," Adv. Opt. Technol. <b>1</b>(1-2), 51-58
     * (2012), <a href="https://doi.org/10.1515/aot-2012-0002">doi:10.1515/aot-2012-0002</a>
     * - 151, 318, 312, 157 and 49 for C, d, e, F and g - divided by the d weight. That
     * paper also weights 404.66 nm at 13 (0.041); this set leaves it out.
     */
    static constexpr double WT_d = 1.0;
    static constexpr double WT_C = 0.475;
    static constexpr double WT_e = 0.98;
    static constexpr double WT_F = 0.49;
    static constexpr double WT_g = 0.15;

    /**
     * The weights of k9Tester(), in its order.
     *
     * Read off the normalized spectral weighting curve in the ZEISS MTF-Tester K9 data
     * sheet and divided by the value at d. The curve is what its halogen lamp at 3200 K,
     * its BG 38 filter and its Hamamatsu R446 photomultiplier come to between 380 and
     * 730 nm; it peaks near 555-560 nm, so e carries slightly more than d.
     *
     * <b>Treat these as approximate.</b> They were digitized from a plot in a single
     * product data sheet - one source, no tabulated values - so a weight is good to a
     * point or two at best. The data sheet is for the K9 service tester; whether the K8
     * that measures published MTF curves weights the same way is not established here.
     *
     * Source https://web.archive.org/web/20071118233344/http://www.zeiss.com/C12567A8003B8B6F/EmbedTitelIntern/ProspektK9_e_PDF/$File/K9_e_Version_021121.pdf
     */
    static constexpr double K9_WT_d = 1.0;
    static constexpr double K9_WT_C = 0.330;
    static constexpr double K9_WT_e = 1.046;
    static constexpr double K9_WT_F = 0.654;
    static constexpr double K9_WT_g = 0.235;

    /** The d line alone. */
    static WavelengthSet dLine();

    /** The e line alone. */
    static WavelengthSet eLine();

    /** d, F and C at equal weight. */
    static WavelengthSet standardTriColor();

    /** d, C, e, F and g at the standard weights; see hasStandardWeights(). */
    static WavelengthSet standardFiveColor();

    /**
     * The same five lines weighted as the ZEISS MTF-Tester K9 weights the spectrum, which
     * is broader than a photopic response - at g it passes 0.235 against photopic's 0.018.
     * Useful for comparing against measured MTF curves rather than for design. See
     * #K9_WT_d for where the numbers come from and how far to trust them.
     */
    static WavelengthSet k9Tester();

    /**
     * @param wavelengths in nanometres, at least one, each positive
     * @param weights     one per wavelength, each finite and non-negative
     */
    static WavelengthSet of(const std::vector<double> &wavelengths,
                            const std::vector<double> &weights);

    /** The set a prescription was built with. */
    static WavelengthSet of(const Prescription &prescription);

    /**
     * @param letters Fraunhofer letters, as letters() writes them
     * @param weights one per letter, or empty for equal weights of 1.0 (the Java's null)
     */
    static WavelengthSet ofLetters(const std::vector<std::string> &letters,
                                   const std::optional<std::vector<double>> &weights = std::nullopt);

    /** @throws IllegalArgumentException naming the lines that can be used */
    static double wavelengthOf(const std::string &letter);

    /** The letter for a wavelength, or empty if no spectral line is that wavelength. */
    static std::optional<std::string> letterOf(double wavelength);

    const std::vector<double> &wavelengths() const { return _wavelengths; }

    const std::vector<double> &weights() const { return _weights; }

    int size() const { return static_cast<int>(_wavelengths.size()); }

    /** True when every weight is 1.0, the case a trial writes without a weights clause. */
    bool hasEqualWeights() const;

    /** True when this is standardFiveColor(), the case a trial writes as `standard`. */
    bool hasStandardWeights() const;

    /**
     * The wavelengths as Fraunhofer letters, in order.
     *
     * @throws IllegalStateException if some wavelength is not a Fraunhofer line, which a
     *                               trial has no way to spell
     */
    std::string letters() const;

    bool operator==(const WavelengthSet &that) const;
    bool operator!=(const WavelengthSet &that) const { return !(*this == that); }

    std::string toString() const;

private:
    WavelengthSet(std::vector<double> wavelengths, std::vector<double> weights)
        : _wavelengths(std::move(wavelengths)), _weights(std::move(weights)) {}

    std::vector<double> _wavelengths;
    std::vector<double> _weights;
};

} // namespace redukti::spec

#endif
