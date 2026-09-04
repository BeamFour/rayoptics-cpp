// The Nikkor Z 58mm f/0.95 S model that MtfTest and ContrastAnalysisTest both
// build, from org.redukti.rayoptics.analysis.MtfTest.buildTestModel.
//
// Shared because the Java shares it -- ContrastAnalysisTest calls
// MtfTest.buildTestModel() directly -- and because the golden MTF and contrast
// values in both depend on this exact prescription, down to the stop
// semi-diameter.
#ifndef REDUKTI_TESTS_NIKKORZ58MODEL_H
#define REDUKTI_TESTS_NIKKORZ58MODEL_H

#include "redukti/rayoptics/optical/OpticalModel.h"

#include <memory>

namespace redukti::test {

std::unique_ptr<rayoptics::optical::OpticalModel> buildNikkorZ58TestModel();

} // namespace redukti::test

#endif // REDUKTI_TESTS_NIKKORZ58MODEL_H
