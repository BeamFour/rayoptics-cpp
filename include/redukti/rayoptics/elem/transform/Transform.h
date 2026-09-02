// C++ port of org.redukti.rayoptics.elem.transform.Transform
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_ELEM_TRANSFORM_TRANSFORM_H
#define REDUKTI_RAYOPTICS_ELEM_TRANSFORM_TRANSFORM_H

#include "redukti/rayoptics/math/Tfm3d.h"
#include "redukti/rayoptics/raytr/RayTypes.h"
#include "redukti/rayoptics/seq/Gap.h"
#include "redukti/rayoptics/seq/Interface.h"
#include "redukti/rayoptics/util/Tuples.h"

#include <memory>
#include <optional>
#include <vector>

namespace redukti::rayoptics::seq {
class SequentialModel;
}

namespace redukti::rayoptics::elem::transform {

using IfcGapPair =
    util::Pair<std::shared_ptr<seq::Interface>, std::shared_ptr<seq::Gap>>;

class Transform {
public:
    /**
     * Return global surface coordinates (rot, t) wrt the `glo` surface.
     *
     * `origin` is the transform from the desired global origin to the global
     * surface `glo`; null means the identity.
     */
    static std::vector<math::Tfm3d> compute_global_coords(
        seq::SequentialModel *seq_model, std::optional<int> glo,
        std::optional<math::Tfm3d> origin);

    /** Return forward surface coordinates (r.T, t) for each interface. */
    static std::vector<math::Tfm3d> compute_local_transforms(
        seq::SequentialModel *seq_model, const std::vector<IfcGapPair> *seq, int step);

    static math::Tfm3d forward_transform(const seq::Interface &s1, double zdist,
                                         const seq::Interface &s2);

    static math::Tfm3d reverse_transform(const seq::Interface &s2, double zdist,
                                         const seq::Interface &s1);

    /** Transform a ray segment from after a surface back to before it. */
    static raytr::RayData transform_after_surface(const seq::Interface &ifc,
                                                  const raytr::RayData &ray_seg);
};

} // namespace redukti::rayoptics::elem::transform

#endif // REDUKTI_RAYOPTICS_ELEM_TRANSFORM_TRANSFORM_H
