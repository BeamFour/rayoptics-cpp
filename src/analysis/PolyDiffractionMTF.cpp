// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.analysis.PolyDiffractionMTF.
#include "redukti/analysis/PolyDiffractionMTF.h"

#include "redukti/Exceptions.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/rayoptics/specs/WvlSpec.h"

#include <cmath>

namespace redukti::analysis {

PolyDiffractionMTF::PolyDiffractionMTF(rayoptics::optical::OpticalModel *opm,
                                       rayoptics::specs::Field &fld_,
                                       const std::vector<double> &wvls,
                                       const std::vector<double> &weights_, double foc,
                                       int n_pupil, int maxdim)
    : fld(&fld_), image_pt_2d(0.0, 0.0) {
    if (wvls.size() != weights_.size())
        throw IllegalArgumentException("one weight is needed per wavelength");
    image_pt_2d = DiffractionMTF::image_point(opm, fld_, opm->optical_spec->wvls->central_wvl(), foc);
    for (std::size_t i = 0; i < wvls.size(); i++) {
        if (weights_[i] == 0.0) continue;
        mono.emplace_back(opm, fld_, wvls[i], foc, n_pupil, maxdim, image_pt_2d);
        weights.push_back(weights_[i]);
    }
    if (mono.empty())
        throw IllegalArgumentException("at least one wavelength must have a non-zero weight");
}

double PolyDiffractionMTF::sagittal(double nu) const { return along(nu, true); }

double PolyDiffractionMTF::tangential(double nu) const { return along(nu, false); }

double PolyDiffractionMTF::along(double nu, bool x) const {
    double re = 0, im = 0, sum = 0;
    for (std::size_t i = 0; i < mono.size(); i++) {
        auto otf = x ? mono[i].sagittal_otf(nu) : mono[i].tangential_otf(nu);
        re += weights[i] * otf[0];
        im += weights[i] * otf[1];
        sum += weights[i];
    }
    return std::hypot(re, im) / sum;
}

} // namespace redukti::analysis
