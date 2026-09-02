// C++ port of org.redukti.mathlib.Transform3
#include "redukti/mathlib/Transform3.h"

namespace redukti::mathlib {

std::string Transform3::toString() const {
    return "{translation=" + this->translation.toString() + ",rmat=" +
           this->rotation_matrix.toString() +
           ",use_rmat=" + (this->use_rotation_matrix ? "true" : "false") + "}";
}

} // namespace redukti::mathlib
