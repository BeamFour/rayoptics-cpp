// C++ port of org.redukti.rayoptics.parax.ParaxModel
#ifndef REDUKTI_RAYOPTICS_PARAX_PARAXMODEL_H
#define REDUKTI_RAYOPTICS_PARAX_PARAXMODEL_H

#include "redukti/rayoptics/parax/ParaxTypes.h"
#include "redukti/rayoptics/seq/SurfaceData.h"
#include "redukti/rayoptics/util/Tuples.h"

#include <optional>
#include <vector>

namespace redukti::rayoptics::optical {
class OpticalModel;
}
namespace redukti::rayoptics::seq {
class SequentialModel;
}

namespace redukti::rayoptics::parax {

class ParaxModel {
public:
    /** Back-reference to the owning OpticalModel; borrowed. */
    optical::OpticalModel *opt_model = nullptr;
    std::optional<double> opt_inv;
    /**
     * A sibling of this object, not a child: the OpticalModel owns both, and
     * the Java caches the reference in the constructor. Borrowed, never owned.
     */
    seq::SequentialModel *seq_model = nullptr;

    class ParaxItem {
    public:
        double power;
        double reduced_thickness;
        double index;
        seq::InteractMode refract_mode;

        ParaxItem(double power_, double reduced_thickness_, double index_,
                  seq::InteractMode refract_mode_)
            : power(power_), reduced_thickness(reduced_thickness_), index(index_),
              refract_mode(refract_mode_) {}
    };

    std::vector<ParaxItem> sys;
    std::vector<ParaxComponent> ax;
    std::vector<ParaxComponent> pr;

    ParaxModel(optical::OpticalModel *opt_model_, std::optional<double> opt_inv_);

    static std::vector<ParaxItem> seq_path_to_paraxial_lens(
        const std::vector<seq::PathSeg> &path);

    void update_model();

    /**
     * Returns ((min lower vignetting, its surface index),
     *          (min upper vignetting, its surface index)). The indices are
     * null in the Java when no surface vignettes.
     */
    util::Pair<util::Pair<double, std::optional<int>>,
               util::Pair<double, std::optional<int>>>
    paraxial_vignetting(std::optional<double> rel_fov);

private:
    void build_lens();
};

} // namespace redukti::rayoptics::parax

#endif // REDUKTI_RAYOPTICS_PARAX_PARAXMODEL_H
