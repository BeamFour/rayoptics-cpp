// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.tools.AsphericOptimizationTest.
//
// A radial asphere stores its omitted odd powers as zeros, so a coefficient index is the
// power it multiplies and not the column it came from. These check that the optimizer
// varies the term the trial names, and that a trial or pipeline written back out keeps
// every other coefficient where it was.
#include "TestHarness.h"

#include "redukti/importers/OpticalBenchDataImporter.h"
#include "redukti/optim/OptimizationBuilder.h"
#include "redukti/optim/OptimizationTrial.h"
#include "redukti/optim/Var.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/spec/Prescription.h"
#include "redukti/tools/LensTool2.h"
#include "redukti/util/Args.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

namespace {

namespace fs = std::filesystem;

using redukti::importers::OpticalBenchDataImporter;
using redukti::optim::OptimizationBuilder;
using redukti::optim::OptimizationTrial;
using redukti::optim::Var;
using redukti::optim::VarAsphCoeff;
using redukti::optim::VarAsphK;
using redukti::spec::Prescription;
using redukti::tools::LensTool2;
using redukti::util::Args;

const std::string LENS = "[descriptive data]\n"
                         "title\tMixed asphere regression\n"
                         "[constants]\n"
                         "AsphericalOddCount\t1\n"
                         "[variable distances]\n"
                         "Focal Length\t50.85\n"
                         "F-Number\t5\n"
                         "Angle of View\t10\n"
                         "[lens data]\n"
                         "1\t50\t5\t1.5\t12\t50\n"
                         "2\t-50\t48\t\t12\n"
                         "[aspherical data]\n"
                         "1\t50\t0.1\t1e-6\t2e-7\t3e-9\t4e-11\t5e-13\t6e-15\n";

Prescription prescription(const std::string &text) {
    OpticalBenchDataImporter::LensSpecifications specs;
    specs.parse_buffer(text);
    return Prescription::build_prescription(specs, false);
}

/** The coefficient indices the setup varies, in order. */
std::vector<int> indices(const std::vector<std::shared_ptr<Var>> &variables) {
    std::vector<int> out;
    for (const auto &v : variables) {
        if (const auto *coeff = dynamic_cast<const VarAsphCoeff *>(v.get()))
            out.push_back(coeff->_index);
    }
    return out;
}

std::string readFile(const std::string &path) {
    std::ifstream in(path, std::ios::binary);
    CHECK(in.good());
    return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
}

bool contains(const std::string &text, const std::string &part) {
    return text.find(part) != std::string::npos;
}

} // namespace

TEST(aspheric_existing_means_nonzero_normalized_terms_at_build_time) {
    auto p = prescription(LENS);
    auto builder = OptimizationBuilder::builder(&p)
                       .fields(std::vector<double>{0})
                       .mtfFrequencies(std::vector<int>{20})
                       .varyExistingAspherics()
                       .rayAberrationGoals();
    auto variables = builder.build().variables();
    CHECK_EQ(variables.size(), static_cast<std::size_t>(7));
    CHECK(dynamic_cast<const VarAsphK *>(variables[0].get()) != nullptr);
    CHECK(indices(variables) == (std::vector<int>{2, 3, 5, 7, 9, 11}));

    // A supplied zero and an omitted power are both excluded, as is K=0.
    (*p.get_surfaces()[0]._coeffs)[5] = 0;
    p.get_surfaces()[0]._k = 0;
    auto rebuilt = builder.build().variables();
    CHECK(indices(rebuilt) == (std::vector<int>{2, 3, 7, 9, 11}));
    CHECK_EQ(rebuilt.size(), static_cast<std::size_t>(5));
    CHECK(indices(variables) == (std::vector<int>{2, 3, 5, 7, 9, 11}));
}

TEST(aspheric_coefficient_variable_changes_a6_and_analysis_rebuilds_the_sag) {
    auto p = prescription(LENS);
    auto setup = OptimizationBuilder::builder(&p)
                     .fields(std::vector<double>{0})
                     .mtfFrequencies(std::vector<int>{20})
                     .varyAsphericCoefficient(0, 5)
                     .rayAberrationGoals()
                     .build();
    auto *variable = dynamic_cast<VarAsphCoeff *>(setup.variables()[0].get());
    CHECK(variable != nullptr);
    if (variable == nullptr)
        return;
    CHECK_EQ(variable->get_scaling_factor(), 1e9);
    CHECK_CLOSE(variable->read_from_prescription(), 3, 1e-14);
    setup.analysis()->compute();
    double before = setup.analysis()->_opt_model->seq_model->ifcs[1]->profile->sag(0, 2);
    std::vector<double> expected = *p.get_surfaces()[0]._coeffs;
    variable->set_scaled_value(4);
    variable->write_to_prescription();
    expected[5] = 4e-9;
    const auto &actual = *p.get_surfaces()[0]._coeffs;
    CHECK_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size() && i < expected.size(); i++)
        CHECK_CLOSE(actual[i], expected[i], 1e-24);
    setup.analysis()->compute();
    double after = setup.analysis()->_opt_model->seq_model->ifcs[1]->profile->sag(0, 2);
    CHECK_CLOSE(after - before, 1e-9 * std::pow(2, 6), 1e-16);
}

TEST(aspheric_explicit_terms_override_existing_and_introduce_an_omitted_odd_power) {
    std::string text = LENS +
                       "\n[trial 1]\nfields 0\nfrequencies 20\nvary aspherics existing\n"
                       "vary aspherics 0 4\ngoal spot-rms 1\n";
    auto builder = OptimizationTrial::read(text, 1, false);
    auto setup = builder.builder.build();
    CHECK_EQ(setup.variables().size(), static_cast<std::size_t>(1));
    auto *variable = dynamic_cast<VarAsphCoeff *>(setup.variables()[0].get());
    CHECK(variable != nullptr);
    if (variable == nullptr)
        return;
    CHECK_EQ(variable->_index, 4);              // A5, not the fifth input column
    CHECK_EQ(variable->get_scaling_factor(), 1e4); // round(log10(6^5))
    variable->set_scaled_value(0.01);
    variable->write_to_prescription();
    std::string saved;
    builder.prescription->to_opt_bench_str(saved);
    saved += "\n" + builder.builder.toTrial(1);
    auto restored = OptimizationTrial::read(saved, 1, false);
    CHECK_EQ(restored.prescription->_aspherical_odd_count, 2);
    const auto &before = *builder.prescription->get_surfaces()[0]._coeffs;
    const auto &after = *restored.prescription->get_surfaces()[0]._coeffs;
    CHECK_EQ(after.size(), before.size());
    for (std::size_t i = 0; i < after.size() && i < before.size(); i++)
        CHECK_EQ(after[i], before[i]);
    CHECK(indices(restored.builder.build().variables()) == (std::vector<int>{4}));
}

TEST(aspheric_lenstool_trials_and_pipelines_preserve_coefficient_powers) {
    std::string text = LENS + "\n"
                              "[trial 1]\n"
                              "fields 0\n"
                              "frequencies 20\n"
                              "vary aspherics 0 5\n"
                              "goal spot-rms 1\n"
                              "goal spot sampling hexapolar 2\n"
                              "solver max-evaluations 4\n"
                              "[pipeline 2]\n"
                              "trials 1 1\n";
    auto initialPrescription = prescription(LENS);
    std::vector<double> initial = *initialPrescription.get_surfaces()[0]._coeffs;
    fs::path dir = fs::temp_directory_path() / "rayoptics-aspheric-optimization";
    std::error_code ec;
    fs::remove_all(dir, ec);
    fs::create_directories(dir, ec);
    for (int number : {1, 2}) {
        fs::path input = dir / ("mixed" + std::to_string(number) + ".txt");
        std::ofstream out(input, std::ios::binary | std::ios::trunc);
        out << text;
        out.close();
        Args args = Args::parseArguments({"--specfile", input.string(), "--optimize",
                                          std::to_string(number)});
        std::string output = LensTool2::runOptimizationTrial(text, args);
        CHECK_STR_EQ(readFile(*args.specfile), output);
        auto restored = prescription(output);
        CHECK_EQ(restored._aspherical_odd_count, 1);
        const auto &actual = *restored.get_surfaces()[0]._coeffs;
        CHECK_EQ(actual.size(), initial.size());
        CHECK(std::isfinite(actual[5]));
        CHECK(actual[5] != initial[5]);
        for (std::size_t i = 0; i < initial.size() && i < actual.size(); i++) {
            if (i != 5)
                CHECK_EQ(actual[i], initial[i]);
        }
        CHECK(indices(OptimizationTrial::read(output, 1, false).builder.build().variables()) ==
              (std::vector<int>{5}));
        CHECK(contains(LensTool2::startREADME(restored), " P12 |"));
    }
    fs::remove_all(dir, ec);
}
