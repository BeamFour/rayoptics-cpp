// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.analysis.DiffractionMTF.
#ifndef REDUKTI_ANALYSIS_DIFFRACTIONMTF_H
#define REDUKTI_ANALYSIS_DIFFRACTIONMTF_H

#include "redukti/mathlib/Vector2.h"

#include <array>
#include <optional>
#include <vector>

namespace redukti::rayoptics::optical {
class OpticalModel;
}

namespace redukti::rayoptics::specs {
class Field;
}

namespace redukti::rayoptics::raytr {
class RayGrid;
}

namespace redukti::analysis {

/**
 * Monochromatic diffraction MTF of one field, from the Fourier transform of the PSF.
 *
 * The wavefront is sampled on a RayGrid of n_pupil rays across the vignetted pupil,
 * its PSF computed by Analyses::calc_psf on a maxdim grid, and the MTF is the modulus
 * of the Fourier transform of the PSF, normalized to 1 at zero frequency. Sagittal is
 * the frequency along x, tangential along y.
 *
 * Frequency scale. The grid spans the pupil in n_pupil - 1 intervals in each
 * direction, so a shift of k intervals is k / (n_pupil - 1) of the pupil diameter and
 * corresponds to that fraction of the cutoff frequency. The cutoff is taken per
 * direction from the real rim rays: 2 n NA / lambda, where NA is half the difference
 * of the image-space direction cosines of the two rim rays in that direction. Off axis
 * and for anamorphic systems the two cutoffs differ. The grid is uniform in the
 * entrance pupil, so pupil aberration that distorts the exit pupil non-uniformly is
 * not accounted for.
 *
 * The FFT gives the MTF only at multiples of nu_cutoff / (n_pupil - 1).
 * #sagittal_otf and #tangential_otf evaluate the complex OTF at any frequency from the
 * line spread, with the phase relative to the reference image point.
 *
 * The approach, and the analytic reference #diffraction_limited, follow M. J. Hayford's
 * "Diffraction MTF" notebook for ray-optics issue #182, in
 * https://github.com/mjhoptics/ray-optics-notebooks. This is an independent implementation.
 */
class DiffractionMTF {
public:
    /** Borrowed. */
    rayoptics::specs::Field *fld;
    double wvl;
    double foc;
    int n_pupil;
    int maxdim;
    /** OTF phase reference in the image plane; empty for this wavelength's chief ray */
    std::optional<mathlib::Vector2> image_pt_2d;
    /** cycles per system length unit */
    double nu_cutoff_x;
    double nu_cutoff_y;
    /** maxdim x maxdim, first index along x, zero frequency at [maxdim/2][maxdim/2] */
    std::vector<std::vector<double>> mtf;

private:
    /** line spreads along x and y, unit sum, this wavelength's chief-ray image point at maxdim/2 */
    std::vector<double> lsf_x;
    std::vector<double> lsf_y;
    /** Chief-ray image point minus the requested OTF phase reference, in system length units. */
    mathlib::Vector2 reference_displacement{0.0, 0.0};

public:
    /**
     * How much of the PSF window the geometric spot needs; see #spot_to_window.
     * Above 1 the PSF wraps around in its window and the MTF is unreliable.
     */
    double spot_to_window_;

    /** Fraction of the grid rays the PSF window must hold for the MTF to be trusted. */
    static constexpr double WINDOW_RAY_FRACTION = 0.99;

    /**
     * @param n_pupil     rays across the pupil; even and at least 4
     * @param maxdim      FFT size; at least 2 * n_pupil so the MTF does not wrap
     * @param image_pt_2d image point the OTF phase is referenced to, x and y in the image
     *                    plane; empty for this wavelength's chief ray. Pass the same point,
     *                    e.g. #image_point of the reference wavelength, for every
     *                    wavelength of a field so that their OTF phases include lateral colour.
     *                    The PSF is computed about this wavelength's own chief-ray image point;
     *                    the reference displacement is applied analytically to the complex OTF.
     */
    DiffractionMTF(rayoptics::optical::OpticalModel *opm, rayoptics::specs::Field &fld, double wvl,
                   double foc, int n_pupil, int maxdim,
                   std::optional<mathlib::Vector2> image_pt_2d = std::nullopt);

    /** Cutoff frequency along x and y, 2 n NA / lambda, from the rim rays. */
    static std::array<double, 2> cutoff_frequencies(rayoptics::optical::OpticalModel *opm,
                                                     rayoptics::specs::Field &fld, double wvl);

    /**
     * Where the chief ray of wvl meets the image plane shifted by foc: the point RayGrid
     * references the wavefront to when given no image point.
     */
    static mathlib::Vector2 image_point(rayoptics::optical::OpticalModel *opm,
                                        rayoptics::specs::Field &fld, double wvl, double foc);

    /**
     * Extent of the geometric spot relative to the PSF window.
     *
     * The PSF covers (n_pupil - 1) / nu_cutoff in each direction, centred on the
     * chief-ray image point. A ray landing outside that window carries energy that the
     * periodic FFT folds back in from the opposite side. Each grid ray's distance from the
     * chief-ray image point is divided by the window half-width, taking the larger of x and y;
     * the result is the value that #WINDOW_RAY_FRACTION of the rays do not exceed.
     * At most 1 the window holds the spot.
     *
     * Java's field and static method share the name spot_to_window; the field is
     * spot_to_window_ here.
     *
     * @return 0 if no ray traced
     */
    static double spot_to_window(const rayoptics::raytr::RayGrid &grid, const mathlib::Vector2 &chief,
                                 double half_width_x, double half_width_y);

    /**
     * Smallest even n_pupil whose window would hold the spot measured at this one.
     * The window grows with n_pupil - 1; the spot does not change.
     */
    int suggested_n_pupil() const;

    /** MTF of an aberration-free circular pupil; 0 at and beyond the cutoff. */
    static double diffraction_limited(double nu, double nu_cutoff);

    /** Sagittal MTF at frequency nu in cycles per system length unit. */
    double sagittal(double nu) const;

    /** Tangential MTF at frequency nu in cycles per system length unit. */
    double tangential(double nu) const;

    /** MTF of the aberration-free pupil with the same cutoff, along x. */
    double sagittal_diffraction_limited(double nu) const;

    /** MTF of the aberration-free pupil with the same cutoff, along y. */
    double tangential_diffraction_limited(double nu) const;

    /**
     * Sagittal OTF at frequency nu in cycles per system length unit, as {re, im}.
     * Evaluated at nu exactly; see #otf.
     */
    std::array<double, 2> sagittal_otf(double nu) const;

    /**
     * Tangential OTF at frequency nu in cycles per system length unit, as {re, im}.
     * Evaluated at nu exactly; see #otf.
     */
    std::array<double, 2> tangential_otf(double nu) const;

    /**
     * OTF at frequency nu, as {re, im}: the Fourier transform of the line spread,
     * sum lsf[i] exp(-2 pi i nu x_i), with x_i measured from the reference
     * image point at lsf.size() / 2. The phase is relative to that point.
     *
     * The PSF of Analyses::calc_psf is mirrored with respect to the image plane:
     * index i lies at x_i = -(i - lsf.size() / 2) dx. With the sign taken
     * here, a PSF displaced by d in the image has its phase changed by -2 pi nu d.
     *
     * The FFT of the PSF has bins nu_cutoff / (n_pupil - 1) apart, so the PSF
     * samples are (n_pupil - 1) / (maxdim nu_cutoff) apart. At a bin frequency the
     * modulus equals the FFT value; between bins this is the band-limited interpolation of
     * the bins. 0 at and beyond the cutoff.
     *
     * @param lsf unit-sum line spread of maxdim samples
     */
    static std::array<double, 2> otf(const std::vector<double> &lsf, double nu, double nu_cutoff,
                                     int n_pupil);

    /** The PSF integrated across y (along x) or across x (along y), normalized to unit sum. */
    static std::vector<double> line_spread(const std::vector<std::vector<double>> &psf, bool x);

    /**
     * |FFT(psf)| with zero frequency at the centre, normalized to 1 there. The PSF is real,
     * so the modulus is symmetric and the sign convention of the transform does not matter.
     */
    static std::vector<std::vector<double>> modulus_of_transform(
        const std::vector<std::vector<double>> &psf);

private:
    /** Express an intrinsic complex OTF relative to a displaced image-plane origin. */
    static std::array<double, 2> translate_otf(std::array<double, 2> value, double nu,
                                               double displacement);

    /** Linear interpolation along the central row or column, which holds k = 0 at maxdim/2. */
    double along(double nu, double nu_cutoff, bool x) const;
};

} // namespace redukti::analysis

#endif // REDUKTI_ANALYSIS_DIFFRACTIONMTF_H
