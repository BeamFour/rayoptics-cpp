// Shared prescriptions for the ported integration tests.
//
// The Java keeps US003549241Example05 as a standalone builder class because two
// test methods share it; the same split is kept here.
#ifndef REDUKTI_TESTS_INTEGRATIONMODELS_H
#define REDUKTI_TESTS_INTEGRATIONMODELS_H

#include "redukti/rayoptics/optical/OpticalModel.h"

#include <memory>

namespace integration {

/** C++ port of the Java US003549241Example05 builder. */
std::unique_ptr<redukti::rayoptics::optical::OpticalModel> build_US003549241Example05();

} // namespace integration

#endif // REDUKTI_TESTS_INTEGRATIONMODELS_H
