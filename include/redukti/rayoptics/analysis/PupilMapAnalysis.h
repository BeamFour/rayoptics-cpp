// C++ port of org.redukti.rayoptics.analysis.PupilMapAnalysis
#ifndef REDUKTI_RAYOPTICS_ANALYSIS_PUPILMAPANALYSIS_H
#define REDUKTI_RAYOPTICS_ANALYSIS_PUPILMAPANALYSIS_H

#include <optional>
#include <string>
#include <vector>

namespace redukti::rayoptics {
namespace optical {
class OpticalModel;
}
namespace raytr {
class TraceOptions;
}
namespace specs {
class Field;
}
} // namespace redukti::rayoptics

namespace redukti::rayoptics::analysis {

/**
 * Which part of the pupil a field can actually use, measured rather than modelled.
 *
 * A dense grid of relative pupil coordinates is traced with the physical surface
 * apertures checked and the field's vignetting factors NOT applied, so every sample is a
 * raw pupil coordinate and the passed ones map out the bundle the lens really transmits.
 * The factors are then a hypothesis about that region which the map can be compared
 * against - see piecewise_quality and ellipse_quality.
 *
 * The grid deliberately reaches beyond the nominal pupil, because the factors do:
 * `scale = 1 - factor`, so a negative factor asks for samples outside the unit
 * circle, and on a wide angle lens the real bundle is indeed wider than the nominal pupil
 * off axis. A map that stopped at the unit circle could not show that.
 */
class PupilMapAnalysis {
public:
    /** One sampled pupil coordinate and what became of the ray traced through it. */
    struct Sample {
        double x;
        double y;
        bool passed;
        int blocked_by;
    };

    /**
     * How well a mapping of the four vignetting factors describes the measured bundle.
     *
     * sampled: fraction of the mapped region that the lens passes; the rest is rays the
     *          mapping spends on light that is blocked.
     * covered: fraction of the passed bundle that falls inside the mapped region; the
     *          rest is light the mapping never samples.
     */
    struct MappingQuality {
        double sampled;
        double covered;

        /** As the Java record prints itself, for messages. */
        std::string toString() const;
    };

    class PupilMapForField {
    public:
        int fi;
        /** Borrowed from the model the map was measured on. */
        specs::Field *fld;
        /** Half width of the sampled square, in nominal pupil radii. */
        double reach;
        /** Samples per axis; samples is column major, index = i * num_samples + j. */
        int num_samples;
        std::vector<Sample> samples;
        /** Largest passed |x| and |y|, in nominal pupil radii. */
        double passed_x;
        double passed_y;
        /** Fraction of the nominal pupil the lens passes. */
        double nominal_passed;
        MappingQuality piecewise_quality;
        MappingQuality ellipse_quality;

        PupilMapForField(int fi, specs::Field *fld, double reach, int num_samples,
                         std::vector<Sample> samples);

        const Sample &sample(int i, int j) const {
            return samples[static_cast<std::size_t>(i * num_samples + j)];
        }

        /** Pupil coordinate of grid index i or j along an axis. */
        double coordinate(int index) const {
            return -reach + 2 * reach * index / (num_samples - 1.0);
        }
    };

    class PupilMapResult {
    public:
        std::vector<PupilMapForField> maps;

        std::string toString() const;
    };

    /** Default samples per axis: enough to place a boundary to about 1% of the pupil radius. */
    static constexpr int DEFAULT_NUM_SAMPLES = 121;

    static PupilMapResult eval(optical::OpticalModel *opm);

    /**
     * Maps the pupil of each requested field.
     * The square grows by 50% while any sampled boundary ray passes. The sample
     * count stays fixed, so expanding the square increases the grid spacing.
     *
     * @throws IllegalStateException if passing boundary rays remain after eight
     *                               expansions; no truncated map is returned
     *
     * @param num_samples samples per axis across the sampled square
     * @param fields      indices into the field spec, or empty for every field
     */
    static PupilMapResult eval(optical::OpticalModel *opm, int num_samples,
                               const std::optional<std::vector<int>> &fields);

    /** Initial reach enclosing the candidate mappings; tracing may require a larger square. */
    static double reach_for(const specs::Field &fld);

    /** The factor as Field::apply_vignetting uses it; a negative factor expands. */
    static double scale(double factor) { return factor == 0.0 ? 1.0 : 1.0 - factor; }

    /** Whether a pupil coordinate lies in the region the factors sample as ray-optics maps them. */
    static bool inside_piecewise(double x, double y, const specs::Field &fld);

    /**
     * Whether a pupil coordinate lies in the single translated ellipse whose four extremes
     * are the measured factors: `x' = (a - b)/2 + x (a + b)/2` with `a = 1 - vux` and
     * `b = 1 - vlx`, and likewise in y.
     *
     * The same family as the piecewise map and through the same four points, but without
     * the kink on the axes - the experiment recorded in Documentation/OPTIMIZER.md. Kept
     * here so a map can be compared against both candidates whether or not the tracer
     * offers it.
     */
    static bool inside_ellipse(double x, double y, const specs::Field &fld);

    static double ellipse_scale(double lower, double upper) {
        return 1.0 - 0.5 * (lower + upper);
    }

    static double ellipse_offset(double lower, double upper) { return 0.5 * (lower - upper); }

private:
    static constexpr int MAX_REACH_EXPANSIONS = 8;

    static std::vector<Sample> sample_grid(optical::OpticalModel *opm, specs::Field &fld,
                                           double wvl, const raytr::TraceOptions &options,
                                           double reach, int num_samples);
};

} // namespace redukti::rayoptics::analysis

#endif // REDUKTI_RAYOPTICS_ANALYSIS_PUPILMAPANALYSIS_H
