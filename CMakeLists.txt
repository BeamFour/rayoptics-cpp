cmake_minimum_required(VERSION 3.20)
project(rayoptics_cpp LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (MSVC)
    add_compile_options(/W3 /permissive- /EHsc)
else()
    add_compile_options(-Wall -Wextra)
endif()

# ---------------------------------------------------------------------------
# Ryu (Apache-2.0 / Boost-1.0) supplies shortest round-tripping decimal digits.
# Vendored rather than using std::to_chars: floating-point to_chars landed very
# late in libc++, so Apple Clang on macOS cannot be relied on to have it.
# ---------------------------------------------------------------------------
add_library(ryu STATIC third_party/ryu/d2s.c)
target_include_directories(ryu PUBLIC third_party)
if (NOT MSVC)
    target_compile_options(ryu PRIVATE -w)
endif()

# ---------------------------------------------------------------------------
# The library. Small value types (Vector2/3, pairs, ...) are header-only so the
# arithmetic inlines; anything larger gets a .h/.cpp split.
# ---------------------------------------------------------------------------
add_library(rayoptics
    src/Text.cpp
    src/mathlib/M.cpp
    src/mathlib/Vector2.cpp
    src/mathlib/Vector3.cpp
    src/mathlib/Matrix3.cpp
    src/mathlib/Quaternion.cpp
    src/mathlib/Transform3.cpp
    src/mathlib/RootResult.cpp
    src/mathlib/BrentSolver.cpp
    src/mathlib/SecantSolver.cpp
    src/mathlib/Derivatives.cpp
    src/mathlib/LMLSolver.cpp
    src/mathlib/MinPack.cpp
    src/mathlib/fftpack/ComplexDoubleFFT_Mixed.cpp
    src/mathlib/fftpack/ComplexDoubleFFT.cpp
    src/rayoptics/elem/profiles/Profiles.cpp
    src/rayoptics/seq/Seq.cpp
    src/rayoptics/seq/Glass.cpp
    src/rayoptics/seq/GlassCatalog.cpp
    src/rayoptics/seq/SequentialModel.cpp
    src/rayoptics/specs/Field.cpp
    src/rayoptics/specs/SpecTypes.cpp
    src/rayoptics/specs/WvlSpec.cpp
    src/rayoptics/specs/PupilSpec.cpp
    src/rayoptics/specs/FieldSpec.cpp
    src/rayoptics/specs/OpticalSpecs.cpp
    src/rayoptics/parax/ParaxTypes.cpp
    src/rayoptics/parax/FirstOrder.cpp
    src/rayoptics/parax/ParaxModel.cpp
    src/rayoptics/parax/ThirdOrder.cpp
    src/rayoptics/optical/OpticalModel.cpp
    src/rayoptics/elem/transform/Transform.cpp
    src/rayoptics/raytr/RayTypes.cpp
    src/rayoptics/raytr/RayTrace.cpp
    src/rayoptics/raytr/WaveAbr.cpp
    src/rayoptics/raytr/ExitPupilAiming.cpp
    src/rayoptics/raytr/Trace.cpp
    src/rayoptics/raytr/Wideangle.cpp
    src/rayoptics/raytr/VigCalc.cpp
    src/render/Renderer.cpp
    src/render/RendererSvg.cpp
    src/render/Plot.cpp
    src/data/Interpolated1d.cpp
    src/data/DataSet.cpp
    src/rayoptics/analysis/SpotIntercepts.cpp
    src/rayoptics/analysis/Histogram.cpp
    src/rayoptics/analysis/MTF.cpp
    src/rayoptics/analysis/SpotAnalysis.cpp
    src/rayoptics/analysis/RayAberrations.cpp
    src/rayoptics/analysis/ContrastAnalysis.cpp
)
target_include_directories(rayoptics PUBLIC include)
target_link_libraries(rayoptics PRIVATE ryu)

# ---------------------------------------------------------------------------
# Tests. Zero dependencies: tests/TestHarness.h self-registers each TEST().
# ---------------------------------------------------------------------------
enable_testing()

add_executable(mathlib_tests
    tests/TestMain.cpp
    tests/VectorTest.cpp
    tests/TextTest.cpp
    tests/MatrixTest.cpp
    tests/SolverTest.cpp
    tests/MinPackTest.cpp
    tests/MinPackHybrdBatteryTest.cpp
    tests/MinPackLmderBatteryTest.cpp
    tests/FftTest.cpp
    tests/TrigProbeTest.cpp
    tests/SeqTest.cpp
    tests/RaytrTest.cpp
    tests/AnalysisTest.cpp
    tests/DataTest.cpp
    tests/DecimalFormatTest.cpp
    tests/RenderTest.cpp
    tests/integration/IntegrationModels.cpp
    tests/integration/Leica50mmSummicronRTest.cpp
    tests/integration/Nikkor14mmTest.cpp
    tests/integration/Nikkor58mmTest.cpp
    tests/integration/US003549241Example05Test.cpp
    tests/upstream/JP2018189733Example01PUpstreamTest.cpp
    tests/upstream/JP2026105585Example01UpstreamTest.cpp
    tests/upstream/NikkorZ58mmf095Ex1UpstreamTest.cpp
    tests/upstream/US004123144Example08PUpstreamTest.cpp
    tests/upstream/US004717245Example02PUpstreamTest.cpp
)
target_link_libraries(mathlib_tests PRIVATE rayoptics)
add_test(NAME mathlib COMMAND mathlib_tests)
