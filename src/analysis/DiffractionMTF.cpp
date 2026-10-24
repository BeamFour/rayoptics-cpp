// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.analysis.DiffractionMTF.
#include "redukti/analysis/DiffractionMTF.h"

#include "redukti/Exceptions.h"
#include "redukti/mathlib/FFTUtil.h"
#include "redukti/mathlib/M.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/parax/ParaxTypes.h"
#include "redukti/rayoptics/raytr/Analyses.h"
#include "redukti/rayoptics/raytr/RayGrid.h"
#include "redukti/rayoptics/raytr/Trace.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/rayoptics/specs/PupilSpec.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace redukti::analysis {

using mathlib::FFTUtil;
using mathlib::Vector2;
using mathlib::Vector3;
using mathlib::M::PI;
using rayoptics::optical::OpticalModel;
using rayoptics::raytr::Analyses;
using rayoptics::raytr::RayGrid;
using rayoptics::raytr::RayPkg;
using rayoptics::raytr::Trace;
using rayoptics::raytr::TraceOptions;
using rayoptics::specs::Field;

namespace {

const Vector3 &image_direction(const RayPkg &pkg) { return pkg.ray[pkg.ray.size() - 2].d; }

/** Java's Double.compare ordering: NaN above everything, -0.0 below 0.0. */
bool java_double_less(double a, double b) {
    if (a < b) return true;
    if (a > b) return false;
    if (std::isnan(a) || std::isnan(b)) return !std::isnan(a) && std::isnan(b);
    return std::signbit(a) && !std::signbit(b);
}

} // namespace

DiffractionMTF::DiffractionMTF(OpticalModel *opm, Field &fld_, double wvl_, double foc_,
                               int n_pupil_, int maxdim_, std::optional<Vector2> image_pt_2d_) {
    if (n_pupil_ < 4 || n_pupil_ % 2 != 0)
        throw IllegalArgumentException("n_pupil must be even and at least 4");
    if (maxdim_ < 2 * n_pupil_)
        throw IllegalArgumentException("maxdim must be at least 2 * n_pupil");
    this->fld = &fld_;
    this->wvl = wvl_;
    this->foc = foc_;
    this->n_pupil = n_pupil_;
    this->maxdim = maxdim_;
    this->image_pt_2d = image_pt_2d_;

    auto cutoffs = cutoff_frequencies(opm, fld_, wvl_);
    nu_cutoff_x = cutoffs[0];
    nu_cutoff_y = cutoffs[1];

    reference_displacement = !image_pt_2d_.has_value()
                                 ? Vector2::vector2_0
                                 : image_point(opm, fld_, wvl_, foc_).minus(*image_pt_2d_);
    // Keep lateral colour out of the finite, periodic FFT window. It changes the
    // OTF's reference phase, not the physical image plane or this wavelength's blur.
    RayGrid grid(opm, fld_, wvl_, foc_, std::nullopt, std::nullopt, std::nullopt, std::nullopt,
                 n_pupil_, true, std::numeric_limits<double>::quiet_NaN(), 1.0);
    auto psf = Analyses::calc_psf(grid.grid[2], n_pupil_, maxdim_);
    spot_to_window_ = spot_to_window(grid, image_point(opm, fld_, wvl_, foc_),
                                     (n_pupil_ - 1) / (2 * nu_cutoff_x),
                                     (n_pupil_ - 1) / (2 * nu_cutoff_y));
    mtf = modulus_of_transform(psf);
    lsf_x = line_spread(psf, true);
    lsf_y = line_spread(psf, false);
}

std::array<double, 2> DiffractionMTF::cutoff_frequencies(OpticalModel *opm, Field &fld,
                                                         double wvl) {
    TraceOptions options;
    auto rim_rays = Trace::trace_boundary_rays_at_field(opm, fld, wvl, options);
    // default pupil rays: centre, +x, -x, +y, -y
    const auto &pupil_rays = opm->optical_spec->pupil->pupil_rays;
    int xp = -1, xm = -1, yp = -1, ym = -1;
    for (int i = 0; i < static_cast<int>(pupil_rays.size()); i++) {
        const auto &p = pupil_rays[static_cast<std::size_t>(i)];
        if (p[0] == 1.0 && p[1] == 0.0) xp = i;
        else if (p[0] == -1.0 && p[1] == 0.0) xm = i;
        else if (p[0] == 0.0 && p[1] == 1.0) yp = i;
        else if (p[0] == 0.0 && p[1] == -1.0) ym = i;
    }
    if (xp < 0 || xm < 0 || yp < 0 || ym < 0)
        throw IllegalArgumentException("pupil rays must include the four rim rays at +-1 in x and y");
    for (int i : {xp, xm, yp, ym})
        if (rim_rays[static_cast<std::size_t>(i)] == nullptr)
            throw IllegalArgumentException("a rim ray failed to trace");
    auto ray = [&](int i) -> const RayPkg & { return *rim_rays[static_cast<std::size_t>(i)]; };
    double na_x = 0.5 * std::abs(image_direction(ray(xp)).x - image_direction(ray(xm)).x);
    double na_y = 0.5 * std::abs(image_direction(ray(yp)).y - image_direction(ray(ym)).y);
    double n_img = std::abs(opm->optical_spec->parax_data->fod.n_img);
    double wl = opm->nm_to_sys_units(wvl);
    return {2 * n_img * na_x / wl, 2 * n_img * na_y / wl};
}

Vector2 DiffractionMTF::image_point(OpticalModel *opm, Field &fld, double wvl, double foc) {
    // Hold the package: on a field with no chief ray yet, it is the only owner.
    auto chief_ray_pkg = Trace::get_chief_ray_pkg(opm, fld, wvl, foc);
    const auto &last = chief_ray_pkg->chief_ray->ray.back();
    auto p = last.p.plus(last.d.times(foc / last.d.z));
    return Vector2(p.x, p.y);
}

double DiffractionMTF::spot_to_window(const RayGrid &grid, const Vector2 &chief,
                                      double half_width_x, double half_width_y) {
    std::vector<double> extents;
    for (const auto &row : grid.grid_pkg->grid)
        for (const auto &item : row) {
            if (item.ray_pkg == nullptr) continue;
            const auto &last = item.ray_pkg->ray.back();
            double x = last.p.x + last.d.x * (grid.foc / last.d.z);
            double y = last.p.y + last.d.y * (grid.foc / last.d.z);
            extents.push_back(std::max(std::abs(x - chief.x) / half_width_x,
                                       std::abs(y - chief.y) / half_width_y));
        }
    if (extents.empty()) return 0.0;
    std::sort(extents.begin(), extents.end(), java_double_less);
    int k = static_cast<int>(std::ceil(WINDOW_RAY_FRACTION * static_cast<double>(extents.size()))) - 1;
    return extents[static_cast<std::size_t>(std::max(k, 0))];
}

int DiffractionMTF::suggested_n_pupil() const {
    if (spot_to_window_ <= 1.0) return n_pupil;
    int n = static_cast<int>(std::ceil((n_pupil - 1) * spot_to_window_)) + 1;
    return n % 2 == 0 ? n : n + 1;
}

double DiffractionMTF::diffraction_limited(double nu, double nu_cutoff) {
    double cos_phi = std::abs(nu) / nu_cutoff;
    if (cos_phi >= 1.0) return 0.0;
    double phi = std::acos(cos_phi);
    return (2 / PI) * (phi - std::sin(phi) * cos_phi);
}

double DiffractionMTF::sagittal(double nu) const { return along(nu, nu_cutoff_x, true); }

double DiffractionMTF::tangential(double nu) const { return along(nu, nu_cutoff_y, false); }

double DiffractionMTF::sagittal_diffraction_limited(double nu) const {
    return diffraction_limited(nu, nu_cutoff_x);
}

double DiffractionMTF::tangential_diffraction_limited(double nu) const {
    return diffraction_limited(nu, nu_cutoff_y);
}

std::array<double, 2> DiffractionMTF::sagittal_otf(double nu) const {
    return translate_otf(otf(lsf_x, nu, nu_cutoff_x, n_pupil), nu, reference_displacement.x);
}

std::array<double, 2> DiffractionMTF::tangential_otf(double nu) const {
    return translate_otf(otf(lsf_y, nu, nu_cutoff_y, n_pupil), nu, reference_displacement.y);
}

std::array<double, 2> DiffractionMTF::translate_otf(std::array<double, 2> value, double nu,
                                                    double displacement) {
    if (displacement == 0.0 || (value[0] == 0.0 && value[1] == 0.0)) return value;
    double angle = -2 * PI * nu * displacement;
    double c = std::cos(angle), s = std::sin(angle);
    return {value[0] * c - value[1] * s, value[0] * s + value[1] * c};
}

std::array<double, 2> DiffractionMTF::otf(const std::vector<double> &lsf, double nu,
                                          double nu_cutoff, int n_pupil) {
    if (std::abs(nu) >= nu_cutoff) return {0.0, 0.0};
    int maxdim = static_cast<int>(lsf.size());
    int centre = maxdim / 2;
    double dx = (n_pupil - 1) / (maxdim * nu_cutoff);
    double re = 0, im = 0;
    for (int i = 0; i < maxdim; i++) {
        double angle = 2 * PI * nu * (i - centre) * dx;
        re += lsf[static_cast<std::size_t>(i)] * std::cos(angle);
        im += lsf[static_cast<std::size_t>(i)] * std::sin(angle);
    }
    return {re, im};
}

std::vector<double> DiffractionMTF::line_spread(const std::vector<std::vector<double>> &psf,
                                                bool x) {
    std::size_t n = psf.size();
    std::vector<double> lsf(n, 0.0);
    double sum = 0;
    for (std::size_t i = 0; i < n; i++)
        for (std::size_t j = 0; j < n; j++) {
            lsf[x ? i : j] += psf[i][j];
            sum += psf[i][j];
        }
    for (std::size_t i = 0; i < n; i++)
        lsf[i] /= sum;
    return lsf;
}

double DiffractionMTF::along(double nu, double nu_cutoff, bool x) const {
    double k = std::abs(nu) * (n_pupil - 1) / nu_cutoff;
    if (k >= n_pupil - 1) return 0.0;
    std::size_t centre = static_cast<std::size_t>(maxdim / 2);
    auto lo = static_cast<std::size_t>(std::floor(k));
    double t = k - static_cast<double>(lo);
    double a = x ? mtf[centre + lo][centre] : mtf[centre][centre + lo];
    double b = x ? mtf[centre + lo + 1][centre] : mtf[centre][centre + lo + 1];
    return a + t * (b - a);
}

std::vector<std::vector<double>> DiffractionMTF::modulus_of_transform(
    const std::vector<std::vector<double>> &psf) {
    std::size_t n = psf.size();
    std::size_t shift = n / 2;
    // the PSF as interleaved complex rows
    FFTUtil::ComplexArray x(n, std::vector<double>(2 * n, 0.0));
    for (std::size_t i = 0; i < n; i++)
        for (std::size_t j = 0; j < n; j++)
            x[i][2 * j] = psf[i][j];
    auto F = FFTUtil::fftshift(FFTUtil::fft2(FFTUtil::fftshift(x)));
    // take the modulus, normalize to the zero frequency
    std::vector<std::vector<double>> mtf(n, std::vector<double>(n));
    for (std::size_t i = 0; i < n; i++)
        for (std::size_t j = 0; j < n; j++)
            mtf[i][j] = std::hypot(F[i][2 * j], F[i][2 * j + 1]);
    double dc = mtf[shift][shift];
    if (!(dc > 0) || !std::isfinite(dc))
        throw IllegalArgumentException("PSF must have finite, positive total intensity");
    for (std::size_t i = 0; i < n; i++)
        for (std::size_t j = 0; j < n; j++)
            mtf[i][j] /= dc;
    return mtf;
}

} // namespace redukti::analysis
