// C++ port of org.redukti.optim.SetupAssertions.
//
// The Java walks each object's fields by reflection and prints every plain value it
// finds. There is no reflection here, so each type spells out the same values by hand:
// the numbers, flags and enums a variable or goal holds, without the references to the
// prescription and the analysis, which differ between two setups by identity only.
#ifndef REDUKTI_TESTS_SETUPASSERTIONS_H
#define REDUKTI_TESTS_SETUPASSERTIONS_H

#include "redukti/optim/Analysis.h"
#include "redukti/optim/Goal.h"
#include "redukti/optim/OptimizationBuilder.h"
#include "redukti/optim/Var.h"

#include <string>
#include <vector>

namespace redukti::test {

std::string describe(const optim::Var &variable);
std::string describe(optim::Goal &goal);
std::string describe(const optim::Analysis &analysis);

/** The same variables and goals, in the same order, and the same analysis settings. */
void assertSameSetup(const char *what, optim::OptimizationBuilder::OptimizationSetup &expected,
                     optim::OptimizationBuilder::OptimizationSetup &actual);

} // namespace redukti::test

#endif // REDUKTI_TESTS_SETUPASSERTIONS_H
