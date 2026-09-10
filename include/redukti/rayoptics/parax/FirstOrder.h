// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
// Java version by Dibyendu Majumdar
// See LICENSE-ray-optics.txt
//
// C++ port of org.redukti.rayoptics.parax.FirstOrder
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
    /**
     * Returns paraxial axial and chief rays, plus first order data.
     *
     * @param opt_model
     * @param stop
     * @param wvl
     * @return
     */
    /** Compute the first order properties of the optical model. */
    static std::shared_ptr<ParaxData> compute_first_order(
        optical::OpticalModel *opt_model, std::optional<int> stop, double wvl);

    /**
     * Calculate transfer matrix and inverse from 1st to kth surface.
     */
    /** The paraxial ray transfer matrix and its inverse at surface `kth`. */
    static util::Pair<mathlib::Matrix2, mathlib::Matrix2> get_parax_matrix(
        const std::vector<ParaxComponent> &p_ray,
        const std::vector<ParaxComponent> &q_ray, int kth, double n_k);

    /**
     * Returns paraxial p and q rays, plus partial first order data.
     *
     *     Args:
     *         path: an iterator containing interfaces and gaps to be traced.
     *               for each iteration, the sequence or generator should return a
     *               list containing: **Intfc, Gap, Trfm, Index, Z_Dir**
     *         oal: overall geometric length of the gaps in `path`
     *         n_0: refractive index preceding the first interface
     *         n_k: refractive index following last interface
     *
     *     Returns:
     *         (p_ray, q_ray, (efl, fl_obj, fl_img, pp1, ppk, pp_sep, ffl, bfl))
     *
     *         - p_ray: [ht, slp, aoi], [1, 0, -]
     *         - q_ray: [ht, slp, aoi], [0, 1, -]
     *         - power: optical power of system
     *         - efl: effective focal length
     *         - fl_obj: object space focal length, f
     *         - fl_img: image space focal length, f'
     *         - pp1: distance from the 1st interface to the front principle plane
     *         - ppk: distance from the last interface to the rear principle plane
     *         - pp_sep: distance from the front principle plane to the rear
     *                   principle plane
     *         - ffl: front focal length, distance from the 1st interface to the
     *                front focal point
     *         - bfl: back focal length, distance from the last interface to the back
     *                focal point
     */
    static PrincipalPointsInfo compute_principle_points(
        const std::vector<seq::PathSeg> &path, double oal, std::optional<double> n_0,
        std::optional<double> n_k, std::optional<int> os_idx,
        std::optional<int> is_idx);

    /**
     * Perform a paraxial raytrace of 2 linearly independent rays
     */
    /** Trace the axial and chief paraxial rays through the path. */
    static util::Pair<std::vector<ParaxComponent>, std::vector<ParaxComponent>>
    paraxial_trace(const std::vector<seq::PathSeg> &path, int start,
                   const ParaxComponent &start_yu, const ParaxComponent &start_yu_bar);
};

} // namespace redukti::rayoptics::parax

#endif // REDUKTI_RAYOPTICS_PARAX_FIRSTORDER_H
