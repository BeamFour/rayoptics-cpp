// C++ port of the org.redukti.optim.Goal subclasses: GoalParax, GoalSpotRMS,
// GoalSpotMaxRadius, GoalSpotDeviation, GoalRayAberration, GoalGeoMTF,
// GoalMTFProxy, GoalContrast and GoalContrastBalance.
#ifndef REDUKTI_OPTIM_GOALS_H
#define REDUKTI_OPTIM_GOALS_H

#include "redukti/optim/Goal.h"

#include <vector>

namespace redukti::optim {

/**
 * Paraxial goals are helpful in anchoring the system so that
 * the optimizer does not make massive changes to focal length etc.
 */
class GoalParax : public Goal {
public:
    const int _parax_id;

    GoalParax(Analysis *analysis, int paraxId, double target, double weight)
        : Goal(analysis, target, weight), _parax_id(paraxId) {}

    double value() override;
    std::string toString() override;
};

class GoalSpotRMS : public Goal {
public:
    const int _field;

    GoalSpotRMS(Analysis *analysis, int field, double target, double weight)
        : Goal(analysis, target, weight), _field(field) {}

    double value() override;
    std::string toString() override;
};

class GoalSpotMaxRadius : public Goal {
public:
    const int _field;

    GoalSpotMaxRadius(Analysis *analysis, int field, double target, double weight)
        : Goal(analysis, target, weight), _field(field) {}

    double value() override;
    std::string toString() override;
};

/** One signed, Gaussian-weighted image-plane ray deviation for RMS spot optimization. */
class GoalSpotDeviation : public Goal {
public:
    /** Zero based, as stored; the constructor takes a one based field. */
    const int _field;
    const int _wavelength_index;
    const int _sample_index;
    const int _orientation;

    GoalSpotDeviation(Analysis *analysis, int field, int wavelength_index,
                      int sample_index, int orientation, double weight);

    double value() override;
    std::string toString() override;
};

/**
 * Ray aberration for a field / orientation / pos / wvl.
 */
class GoalRayAberration : public Goal {
public:
    const int _field;
    const int _orientation;
    const int _pos;
    const double _wvl;

    GoalRayAberration(Analysis *analysis, int field, int orientation, int pos, double wvl,
                      double target, double weight);

    double value() override;
    std::string toString() override;
};

/**
 * Targets sagittal or tangential MTF at given field and freq.
 */
class GoalGeoMTF : public Goal {
public:
    const int _freq;
    const int _orientation;
    const int _field;

    /**
     * @param field Fields start at 1
     * @param orientation Orientation::SAGITTAL or Orientation::TANGENTIAL
     * @param freq  MTF frequency
     */
    GoalGeoMTF(Analysis *analysis, int field, int orientation, int freq, double target,
               double weight);

    double value() override;
    std::string toString() override;
};

/**
 * MTF proxy for a field / orientation / pos / wvl.
 * Based on Kidger paper.
 */
class GoalMTFProxy : public Goal {
public:
    const int _field;
    const int _orientation;
    const int _pos;
    const double _wvl;
    const double _freq;

    GoalMTFProxy(Analysis *analysis, int field, int orientation, int pos, double wvl,
                 double freq, double target, double weight);

    double value() override;
    std::string toString() override;
};

/** A single weighted pupil phase-difference residual for contrast optimization. */
class GoalContrast : public Goal {
public:
    const int _contrast_index;
    const int _frequency;
    const int _field;
    const int _wavelength_index;
    const int _sample_index;
    const int _orientation;

    GoalContrast(Analysis *analysis, int contrast_index, int frequency, int field,
                 int wavelength_index, int sample_index, int orientation, double weight);

    double value() override;
    std::string toString() override;
};

/**
 * Holds sagittal and tangential contrast in balance at one field and frequency.
 *
 * The contrast merit minimizes sum(sagittal^2) + sum(tangential^2), which at a
 * fixed total barely discriminates how astigmatism is split between the two
 * meridians. A designer discriminates sharply: on the Leica 75/2 a solve
 * produced 0.156 sagittal against 0.725 tangential at field 0.8, 50 cyc/mm,
 * while the sum of the two MTFs stayed within 15% of its value everywhere else
 * across the field. The lens was not worse in that zone, it was lopsided, and
 * nothing in the merit had an opinion about that.
 *
 * This supplies the opinion. The value is the difference between what the two
 * orientations contribute to the merit,
 *
 *   sum_wavelengths w * ( w_sagittal   * sum_samples r_sagittal^2
 *                       - w_tangential * sum_samples r_tangential^2 )
 *
 * against a target of zero. Defining it on the residuals rather than on the raw
 * wavefront differences means it automatically follows whatever those residuals
 * already account for - residual centering, quadrature weights, wavelength and
 * orientation weights, and the frequency calibration.
 *
 * It is smooth and quadratic in the wavefront differences, with no modulus and
 * no square root, so unlike the phasor MTF goal that was tried and reverted it
 * has no kink to fall into.
 *
 * The orientation weights are the ones the ordinary contrast goals use, so the
 * ratio w_sagittal / w_tangential IS the instruction for how the two meridians
 * should differ: balance is reached when the weighted contributions match, not
 * when the two residual energies do. Weights of 0.5 and 0.1 ask for tangential
 * to carry five times the energy of sagittal, and the goal will deliver that.
 *
 * NOT MEANINGFUL ON AXIS. At field zero the meridians are identical by
 * rotational symmetry, so there is nothing to balance and the value reduces to
 * (w_sagittal - w_tangential) * S, where S is the axial residual energy. Equal
 * weights make that exactly zero and the goal inert; unequal weights turn it
 * into a second axial contrast goal whose strength is a number that usually fell
 * out of a field taper rather than a decision. On the Leica 75/2 with weights 8
 * and 4 it came to 6.3% of the merit, none of it balance. It is also shaped
 * unlike the goals it shadows - S is already a sum of squares, so this residual
 * is quadratic where the per-sample ones are linear, and it fades quadratically
 * as the design improves. Prefer raising the field-zero contrast weights if
 * axial emphasis is what is wanted.
 *
 * Two further things to be clear about. This is a design preference, not a
 * correction: it tells the optimizer something it cannot infer, rather than
 * fixing an error. And it is one residual against the thousands in a contrast
 * merit, so its weight has to be set deliberately - see
 * OptimizationBuilder::contrastBalanceGoals.
 */
class GoalContrastBalance : public Goal {
public:
    const int _contrast_index;
    const int _frequency;
    /** Zero based, as stored; the constructor takes a one based field. */
    const int _field;

    /**
     * @param field              one based, as for every other field-addressed goal
     * @param wavelengthWeights  one per wavelength, pooling the per-wavelength blocks
     */
    GoalContrastBalance(Analysis *analysis, int contrast_index, int frequency, int field,
                        const std::vector<double> &wavelengthWeights, double weight)
        : GoalContrastBalance(analysis, contrast_index, frequency, field,
                              wavelengthWeights, 1.0, 1.0, weight) {}

    /**
     * @param sagittalWeight   field/orientation weight used by the sagittal contrast goals
     * @param tangentialWeight field/orientation weight used by the tangential contrast goals
     */
    GoalContrastBalance(Analysis *analysis, int contrast_index, int frequency, int field,
                        const std::vector<double> &wavelengthWeights,
                        double sagittalWeight, double tangentialWeight, double weight);

    double value() override;
    std::string toString() override;

private:
    std::vector<double> _wavelength_weights;
    double _sagittal_weight;
    double _tangential_weight;
};

} // namespace redukti::optim

#endif // REDUKTI_OPTIM_GOALS_H
