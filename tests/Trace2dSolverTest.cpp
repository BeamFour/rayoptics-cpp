// C++ port of org.redukti.rayoptics.raytr.Trace2dSolverTest.
//
// The 2D aiming solver, on a single plane or weakly curved surface where the
// answer is known in closed form.
#include "TestHarness.h"

#include "redukti/mathlib/Vector3.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/raytr/RayTypes.h"
#include "redukti/rayoptics/raytr/Trace.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/seq/SurfaceData.h"

#include <memory>
#include <vector>

namespace {

using redukti::mathlib::Vector3;
using redukti::rayoptics::optical::OpticalModel;
using redukti::rayoptics::raytr::RayResultWithStartCoord;
using redukti::rayoptics::raytr::Trace;
using redukti::rayoptics::seq::SequentialModel;
using redukti::rayoptics::seq::SurfaceData;

/**
 * The Java returns the SequentialModel and lets the collector keep its
 * OpticalModel alive; the model owns the seq model here, so the caller has to
 * hold it.
 */
std::unique_ptr<OpticalModel> surfaceModel(double curvature) {
    auto opticalModel = std::make_unique<OpticalModel>();
    auto seqModel = opticalModel->seq_model.get();
    seqModel->gaps[0]->thi = 100.0;
    SurfaceData sd(curvature, 10.0);
    sd.max_aperture(100.0);
    seqModel->add_surface(sd);
    seqModel->do_apertures = false;
    seqModel->update_model();
    return opticalModel;
}

std::unique_ptr<OpticalModel> planeSurfaceModel() { return surfaceModel(0.0); }

void checkStartCoords(const RayResultWithStartCoord &result,
                      const std::vector<double> &target, double tol) {
    CHECK(result.start_coords.has_value());
    if (!result.start_coords.has_value())
        return;
    CHECK_EQ(static_cast<int>(result.start_coords->size()),
             static_cast<int>(target.size()));
    for (std::size_t i = 0; i < target.size() && i < result.start_coords->size(); i++)
        CHECK_CLOSE((*result.start_coords)[i], target[i], tol);
}

} // namespace

TEST(trace2d_aims_an_off_axis_ray_in_both_coordinates) {
    auto opm = planeSurfaceModel();
    std::vector<double> target{2.5, -3.75};

    auto result = Trace::get_2d_solution(opm->seq_model.get(), 1, Vector3::ZERO, 100.0,
                                         587.5618, target, true);

    checkStartCoords(result, target, 1.0e-9);
    CHECK(result.rr.pkg != nullptr);
    if (result.rr.pkg == nullptr)
        return;
    CHECK_CLOSE(result.rr.pkg->ray[1].p.x, target[0], 1.0e-9);
    CHECK_CLOSE(result.rr.pkg->ray[1].p.y, target[1], 1.0e-9);
}

TEST(trace2d_raw_solver_aims_an_off_axis_ray_in_both_coordinates) {
    auto opm = planeSurfaceModel();
    std::vector<double> target{-1.25, 4.5};

    auto result = Trace::get_2d_solution_raw(opm->seq_model->path(), 1, Vector3::ZERO,
                                             100.0, 587.5618, target, true);

    checkStartCoords(result, target, 1.0e-9);
    CHECK(result.rr.pkg != nullptr);
    if (result.rr.pkg == nullptr)
        return;
    CHECK_CLOSE(result.rr.pkg->ray[1].p.x, target[0], 1.0e-9);
    CHECK_CLOSE(result.rr.pkg->ray[1].p.y, target[1], 1.0e-9);
}

TEST(trace2d_agrees_with_the_1d_solver_for_a_symmetric_ray) {
    auto opm = surfaceModel(0.01);
    double yTarget = 3.25;

    auto oneDimensional = Trace::get_1d_solution(opm->seq_model.get(), 1, Vector3::ZERO,
                                                 100.0, 587.5618, yTarget, true);
    auto twoDimensional =
        Trace::get_2d_solution(opm->seq_model.get(), 1, Vector3::ZERO, 100.0, 587.5618,
                               std::vector<double>{0.0, yTarget}, true);

    CHECK(oneDimensional.start_coords.has_value());
    CHECK(twoDimensional.start_coords.has_value());
    if (!oneDimensional.start_coords.has_value() ||
        !twoDimensional.start_coords.has_value())
        return;
    CHECK_CLOSE((*twoDimensional.start_coords)[0], 0.0, 1.0e-12);
    CHECK_CLOSE((*twoDimensional.start_coords)[1], (*oneDimensional.start_coords)[1],
                1.0e-8);
    CHECK_CLOSE(twoDimensional.rr.pkg->ray[1].p.y, oneDimensional.rr.pkg->ray[1].p.y,
                1.0e-8);
}

TEST(trace2d_accepts_an_exact_initial_solution) {
    auto opm = planeSurfaceModel();

    auto result = Trace::get_2d_solution(opm->seq_model.get(), 1, Vector3::ZERO, 100.0,
                                         587.5618, std::vector<double>{0.0, 0.0}, true);

    checkStartCoords(result, {0.0, 0.0}, 0.0);
    CHECK(result.rr.pkg != nullptr);
    if (result.rr.pkg == nullptr)
        return;
    CHECK_EQ(result.rr.pkg->ray[1].p.x, 0.0);
    CHECK_EQ(result.rr.pkg->ray[1].p.y, 0.0);
}
