// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
// Java version by Dibyendu Majumdar
// See LICENSE-ray-optics.txt
//
// C++ port of org.redukti.rayoptics.seq.Gap
#ifndef REDUKTI_RAYOPTICS_SEQ_GAP_H
#define REDUKTI_RAYOPTICS_SEQ_GAP_H

#include "redukti/rayoptics/seq/Medium.h"

#include <memory>
#include <string>

namespace redukti::rayoptics::seq {

/**
 * Gap container class.
 *
 *     The gap class represents the space between 2 surfaces. It contains the
 *     media definition for the space and a (z) displacement between the
 *     adjacent surfaces.
 *
 *     The most common use case is an optical system with surfaces centered on a
 *     common axis. The Gap structure implements this case in the simplest manner.
 *     More complicated transformations between surfaces are implemented using
 *     transformations associated with the surfaces themselves.
 *
 *     Attributes:
 *         thi: the length (along z) of the gap
 *         medium: a :class:`~optical.medium.Medium` or a catalog glass instance
 */
class Gap {
public:
    double thi;
    std::shared_ptr<Medium> medium;

    Gap(double thi_, std::shared_ptr<Medium> medium_)
        : thi(thi_), medium(std::move(medium_)) {}

    Gap() : Gap(0.0, Air::INSTANCE()) {}

    void apply_scale_factor(double scale_factor) { thi *= scale_factor; }

    std::string toString() const;
};

} // namespace redukti::rayoptics::seq

#endif // REDUKTI_RAYOPTICS_SEQ_GAP_H
