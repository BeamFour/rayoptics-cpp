// C++ port of org.redukti.optim.Analysis
#ifndef REDUKTI_OPTIM_ANALYSIS_H
#define REDUKTI_OPTIM_ANALYSIS_H

#include "redukti/rayoptics/analysis/ContrastAnalysis.h"
#include "redukti/rayoptics/analysis/MTF.h"
#include "redukti/rayoptics/analysis/RayAberrations.h"
#include "redukti/rayoptics/analysis/SpotAnalysis.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/spec/Prescription.h"

#include <memory>
#include <optional>
#include <vector>

namespace redukti::optim {

class Analysis {
public:
    /**
     * Borrowed, never owned: the optimizer mutates this prescription in place
     * between computes, and every Var writes back to it, so the Analysis has to
     * see the same object rather than a copy.
     */
    spec::Prescription *_prescription;
    std::vector<double> _fields;
    std::vector<double> _pfo;
    /**
     * Null in the Java when spots were not computed, and GoalSpotDeviation
     * tests for exactly that, so the emptiness has to be distinguishable from
     * a computed-but-empty result.
     */
    std::optional<std::vector<rayoptics::analysis::SpotAnalysisResult::SpotResultsForField>>
        _spots;
    /** Null in the Java when ray aberrations were not computed. */
    std::optional<rayoptics::analysis::RayAberrationResult> _ray_aberrations;
    /** Null in the Java when MTF was not computed. */
    std::optional<std::vector<rayoptics::analysis::MTFResultByFreq>> _mtfs;
    /** Never null in the Java; empty when no contrast frequencies are configured. */
    std::vector<rayoptics::analysis::ContrastAnalysisResult> _contrasts;
    std::vector<int> _freqs;
    std::vector<int> _contrast_freqs;
    const int _scenario;
    int _spot_pattern = rayoptics::analysis::SpotOptions::PATTERN_HEXAPOLAR;
    int _num_rays = 64;
    int _num_rings = 14;
    int _num_spokes = 20;
    /** Normalized central obscuration radius for the ordinary spot quadrature. */
    double _inner_pupil_radius = 0.0;
    bool _append_failed_spot_rays = false;
    /** Whether ordinary spot rays are rejected by physical surface apertures. */
    bool _check_spot_apertures = true;
    int _contrast_num_rings = 3;
    int _contrast_num_spokes = 6;
    /** See ContrastOptions::calibrate_frequency(bool); off by default. */
    bool _contrast_calibrate_frequency = false;
    /** See ContrastOptions::aim_exit_pupil(bool); off by default. */
    bool _contrast_aim_exit_pupil = false;
    /** See ContrastOptions::center_residuals(bool); off by default. */
    bool _contrast_center_residuals = false;
    bool _compute_spots = true;
    bool _compute_ray_aberrations = true;
    bool _compute_mtf = true;
    /** How each rebuilt model establishes its vignetting; see vignetting(VigType). */
    spec::VigType _vig_type = spec::VigType::SetPupil;
    /** See freezing_vignetting(bool); off by default. */
    bool _freeze_vignetting = false;

    static constexpr int NUM_TRANSVERSE_RAYS = 10;

    /**
     * Systems and spot analysis setup for each field.
     *
     * Owned here: the Java lets the garbage collector keep the previous model
     * alive for whatever still points at it, whereas here each compute()
     * replaces it and the old one is destroyed. Nothing outlives a compute.
     */
    std::unique_ptr<rayoptics::optical::OpticalModel> _opt_model;

    /**
     * When optimizing a prescription that has multiple scenarios configured
     * use this constructor and set the scenario. At present optimization must
     * be performed for each scenario independently.
     */
    Analysis(spec::Prescription *prescription, std::vector<double> fields,
             std::vector<int> freqs, int scenario);
    virtual ~Analysis() = default;
    Analysis(spec::Prescription *prescription, std::vector<double> fields,
             std::vector<int> freqs)
        : Analysis(prescription, std::move(fields), std::move(freqs), 0) {}

    Analysis &using_gauss_quadrature_pattern(int num_rings, int num_spokes) {
        return using_gauss_quadrature_pattern(num_rings, num_spokes, 0.0);
    }
    Analysis &using_gauss_quadrature_pattern(int num_rings, int num_spokes,
                                             double innerPupilRadius);
    Analysis &retaining_failed_spot_rays(bool value) {
        _append_failed_spot_rays = value;
        return *this;
    }
    Analysis &checking_spot_apertures(bool value) {
        _check_spot_apertures = value;
        return *this;
    }
    Analysis &using_hexapolar_pattern(int num_rays) {
        _spot_pattern = rayoptics::analysis::SpotOptions::PATTERN_HEXAPOLAR;
        _num_rays = num_rays;
        return *this;
    }
    Analysis &using_contrast_analysis(const std::vector<int> &frequencies, int num_rings,
                                      int num_spokes) {
        _contrast_freqs = frequencies;
        _contrast_num_rings = num_rings;
        _contrast_num_spokes = num_spokes;
        return *this;
    }
    Analysis &calibrating_contrast_frequency(bool value) {
        _contrast_calibrate_frequency = value;
        return *this;
    }
    Analysis &aiming_contrast_at_exit_pupil(bool value) {
        _contrast_aim_exit_pupil = value;
        return *this;
    }
    Analysis &centering_contrast_residuals(bool value) {
        _contrast_center_residuals = value;
        return *this;
    }
    Analysis &required_analyses(bool spots, bool rayAberrations, bool mtf) {
        _compute_spots = spots || mtf;
        _compute_ray_aberrations = rayAberrations;
        _compute_mtf = mtf;
        return *this;
    }

    /**
     * How every rebuilt optical model establishes its vignetting factors.
     *
     * VigType::SetPupil is the default and what all existing regression values
     * were generated under: it resizes the pupil so the axial marginal ray meets
     * the stop edge, then measures all four factors with real rays.
     * VigType::SetVig measures the factors the same way without the resize, and
     * agrees closely - within 0.005 of pupil half-width and 3-4 MTF decimals on
     * both test lenses.
     *
     * VigType::Paraxial is cheaper but sets only the y factors: a paraxial ray is
     * meridional and says nothing about the sagittal pupil, so x comes out
     * unvignetted at every field. That makes the pupil an ellipse even on axis,
     * where sagittal and tangential MTF must be equal - measured 0.148 apart at
     * 40 cyc/mm on the Leica 75/2 and 0.010 on the Otus. It also means sagittal
     * is optimized over a superset of the real aperture and tangential over a
     * subset.
     */
    Analysis &vignetting(spec::VigType vigType) {
        _vig_type = vigType;
        discard_frozen_vignetting();
        return *this;
    }

    /**
     * Measure the vignetting factors once, then hold them fixed for the rest of
     * the run.
     *
     * Apertures are never optimization variables, but vignetting is not
     * therefore constant: it is where rays land on those fixed apertures, and 28
     * of 29 variables on the Leica 75/2 move a factor within a single Jacobian
     * step. The drift is smooth, so it does not corrupt the finite difference,
     * but it does mean the solver differentiates the design and the pupil
     * together - and a more heavily vignetted lens has less aberration, so
     * shrinking the pupil is a way to improve an MTF-like merit that costs
     * nothing in the merit and real light in the lens. Freezing removes that
     * route and gives every iteration the same pupil to be compared on.
     *
     * The cost is staleness: the factors describe the design as it was at
     * capture, so the further the solve travels the more the assumed pupil
     * diverges from the real one. Call discard_frozen_vignetting() between
     * solver restarts to re-measure.
     *
     * With VigType::SetPupil the captured pupil value is held too, since factors
     * measured at one working f/# do not describe another. That pins fod.fno, so
     * a GoalParax on ParaxHelper::Fno becomes inert in that combination.
     */
    Analysis &freezing_vignetting(bool value) {
        _freeze_vignetting = value;
        if (!value)
            discard_frozen_vignetting();
        return *this;
    }

    /** Drop captured factors so the next compute() measures them afresh. */
    Analysis &discard_frozen_vignetting() {
        _frozen_vignetting.reset();
        _frozen_pupil_value.reset();
        return *this;
    }

    /** The frozen factors as [field][vux, vlx, vuy, vly], or null if not frozen yet. */
    std::optional<std::vector<std::vector<double>>> frozen_vignetting() const {
        return _frozen_vignetting;
    }

    /** Virtual because the Java is: LMDerMeritFunctionTest overrides it. */
    virtual void compute();

private:
    /** Captured on the first compute when frozen: [field][vux, vlx, vuy, vly]. */
    std::optional<std::vector<std::vector<double>>> _frozen_vignetting;
    std::optional<double> _frozen_pupil_value;

    std::unique_ptr<rayoptics::optical::OpticalModel> build_model(spec::VigType vigType);
    std::unique_ptr<rayoptics::optical::OpticalModel> build_vignetted_model();
    void capture_vignetting();
};

} // namespace redukti::optim

#endif // REDUKTI_OPTIM_ANALYSIS_H
