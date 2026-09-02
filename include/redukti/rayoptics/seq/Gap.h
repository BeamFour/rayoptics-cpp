// C++ port of org.redukti.rayoptics.seq.Gap
#ifndef REDUKTI_RAYOPTICS_SEQ_GAP_H
#define REDUKTI_RAYOPTICS_SEQ_GAP_H

#include "redukti/rayoptics/seq/Medium.h"

#include <memory>
#include <string>

namespace redukti::rayoptics::seq {

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
