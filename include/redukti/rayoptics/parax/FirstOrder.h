// C++ port of org.redukti.rayoptics.parax.FirstOrder
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_PARAX_FIRSTORDER_H
#define REDUKTI_RAYOPTICS_PARAX_FIRSTORDER_H

#include "redukti/mathlib/Matrix2.h"
#include "redukti/rayoptics/parax/ParaxTypes.h"
#include "redukti/rayoptics/seq/SurfaceData.h"
#include "redukti/rayoptics/util/Tuples.h"

#include <memory>
#include <optional>
#include <vector>

namespace redukti::rayoptics::optical {
class OpticalModel;
}

namespace redukti::rayoptics::parax {

class FirstOrder {
public:
    /** Compute the first order properties of the optical model. */
    static std::shared_ptr<ParaxData> compute_first_order(
        optical::OpticalModel *opt_model, std::optional<int> stop, double wvl);

    /** The paraxial ray transfer matrix and its inverse at surface `kth`. */
    static util::Pair<mathlib::Matrix2, mathlib::Matrix2> get_parax_matrix(
        const std::vector<ParaxComponent> &p_ray,
        const std::vector<ParaxComponent> &q_ray, int kth, double n_k);

    static PrincipalPointsInfo compute_principle_points(
        const std::vector<seq::PathSeg> &path, double oal, std::optional<double> n_0,
        std::optional<double> n_k, std::optional<int> os_idx,
        std::optional<int> is_idx);

    /** Trace the axial and chief paraxial rays through the path. */
    static util::Pair<std::vector<ParaxComponent>, std::vector<ParaxComponent>>
    paraxial_trace(const std::vector<seq::PathSeg> &path, int start,
                   const ParaxComponent &start_yu, const ParaxComponent &start_yu_bar);
};

} // namespace redukti::rayoptics::parax

#endif // REDUKTI_RAYOPTICS_PARAX_FIRSTORDER_H
