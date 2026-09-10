// The software is ported from Goptical, hence is licensed under the GPL.
// Copyright (c) 2021 Dibyendu Majumdar
// Goptical: Copyright (C) 2010-2011 Free Software Foundation, Inc; Author: Alexandre Becoulet
// Licensed under the GNU General Public License, version 3 or later; see LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.mathlib.Transform3
#include "redukti/mathlib/Transform3.h"

namespace redukti::mathlib {

std::string Transform3::toString() const {
    return "{translation=" + this->translation.toString() + ",rmat=" +
           this->rotation_matrix.toString() +
           ",use_rmat=" + (this->use_rotation_matrix ? "true" : "false") + "}";
}

} // namespace redukti::mathlib
