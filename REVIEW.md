# C++ Port Review

Review date: 2026-09-04

This review compares the standalone C++ repository with the current Java
implementation in `Beam43/rayoptics`. It covers the ported library, command-line
tool, build, and tests. It is separate from `PERFORMANCE.md`: this document tracks
correctness, safety, parity, and coverage issues.

## Verification baseline

- The normal C++ Release suite passes in approximately 30 seconds.
- The C++ test executable registers 252 tests, but three optimizer tests return
  immediately unless `RAYOPTICS_RUN_SLOW_TESTS` is enabled.
- The Java Maven suite ran 197 tests with zero failures and one skipped test in
  approximately 3 minutes 25 seconds.
- The skipped Java test is `TestOptimLMder`.
- Java reproduced the stored direct Otus values exactly:
  `initialRms=0.07481204946808326` and
  `finalRms=0.011873429221888566`.
- Enabling the gated C++ Nikkor optimizer test consumed more than fifteen
  CPU-minutes without reaching its first optimizer result, so that diagnostic
  run was stopped.

## Open findings

### Fixed: optimizer recovery missed `std::bad_optional_access`

Locations:

- `src/optim/LMDer.cpp`, in `computeResiduals()`
- `src/optim/LMDer.cpp`, in `evaluate()`
- `tests/OptimSolveTest.cpp`, which records the observed failure

The Java optimizer catches `Exception` while evaluating a trial point. That
includes runtime failures such as `NullPointerException`, and the trial is
converted into an invalid residual vector rather than aborting the solve.

The C++ port catches only `redukti::Exception`. Some intentionally nullable Java
values are represented by `std::optional`, and `.value()` throws
`std::bad_optional_access`. That exception is outside the redukti hierarchy, so
it escapes the optimizer. The disabled Nikkor test records exactly this failure.

Implemented fix:

1. `std::bad_optional_access` is caught alongside `redukti::Exception` at both
   optimizer evaluation boundaries.
2. The exception policy is centralized in `tryEvaluation()` so the paths cannot
   drift.
3. Focused tests verify residual and nudge recovery and ensure an unrelated
   `std::logic_error` remains visible.
4. Longer term, replace `.value()` at expected ray-failure sites with an explicit
   project exception such as `InvalidRayResultException`. This gives the
   optimizer one deliberate exception family to recover from.

Do not catch all `std::exception` at this boundary. That would also hide
`std::bad_alloc` and unrelated programming errors, potentially turning a broken
program into a misleading optimizer result.

### P1: `LMLSolver::gaussj()` can overwrite its stack

Location: `src/mathlib/LMLSolver.cpp`, around the fixed `ik[100]` and `jk[100]`
arrays.

For `N > 100`, C++ writes beyond the arrays. Java also allocates 100 elements,
but Java raises a bounds exception rather than corrupting memory.

Recommended fix: replace both arrays with `std::vector<int>` sized to `N`.
Alternatively reject `N > 100` explicitly, although dynamic arrays remove an
unnecessary limitation.

### Fixed: analysis results contained fragile borrowed and self-referential pointers

Locations:

- `include/redukti/rayoptics/analysis/SpotIntercepts.h`
- `include/redukti/rayoptics/analysis/SpotAnalysis.h`
- `include/redukti/rayoptics/analysis/MTF.h`
- `include/redukti/rayoptics/analysis/ContrastAnalysis.h`
- `include/redukti/rayoptics/raytr/RayTypes.h`

Spot, contrast, and ray-fan results retain raw `Field*` values owned by the
optical model. The model must outlive the results. Java field references keep the
field objects alive, so the C++ lifetime requirement is stronger and is not
consistently visible in the API.

There are also self-reference problems:

- `SpotIntercepts` points into its owner's `trace_results`. Copying the owner
  copies the vector but leaves the pointer referring to the original object.
- `MonochromaticGeometricMTF::mtf` points to its sibling `h2d`. Its implicit
  copy/move operations do not rebind that pointer.
- `PolyChromaticGeometricMTF` has the same sibling-pointer concern after a move
  once `mtf` has been computed.

Implemented fix:

- Removed `SpotIntercepts::trace_data`; centroid operations use owned arrays.
- Removed `MTF::h2d`; the histogram is consumed during construction only.
  Monochromatic copies and moves and polychromatic moves no longer contain
  sibling pointers that need rebinding.
- Spot, contrast and ray-fan results own const `FieldSnapshot` metadata rather
  than borrowing the model's `Field`. Coordinates, weights, vignetting and the
  report label describe the field at result creation and survive model changes
  or destruction. The snapshot has no `FieldSpec` or ray-cache dependencies.
- Added regression coverage for source destruction, spot-result copying and
  MTF owner movement. Consumers needing mutable fields must use the model;
  result metadata intentionally exposes only the snapshot API.
- Java now uses the same `FieldSnapshot` fields (`x`, `y`, `vux`, `vuy`,
  `vlx`, `vly`, `wt`) and captured report label, and likewise removes the
  construction-only grid and histogram references. Both `GeoMTFPlot` APIs accept
  either a snapshot or a live field and retain owned snapshot metadata.

### P2: configuration access uses unchecked indexing

Locations:

- `src/spec/SurfaceType.cpp`
- `src/spec/RayOpticsModelBuilder.cpp`

Negative or oversized scenario/configuration numbers are converted to
`std::size_t` and used with unchecked `operator[]`. Java raises a bounds
exception, whereas C++ invokes undefined behavior.

Recommended fix: introduce one checked configuration accessor and use it for
f-number, angle of view, thickness, and diameter arrays. Test negative, exactly
past-the-end, and mismatched per-surface configuration counts.

### P2: malformed CLI scenario values silently select scenario zero

Location: `src/util/Args.cpp`.

Java uses `Integer.parseInt()`. C++ uses `std::atoi()`, for which
`--scenario rubbish` returns zero and silently selects the base configuration.

Recommended fix: use strict full-string parsing with range checking, consistent
with the existing integer parsers in the project, and add CLI parsing tests for
garbage, trailing characters, overflow, and a missing value.

### P2: default C++ results overstate optimizer test coverage

Locations:

- `tests/OptimSolveTest.cpp`
- `tests/OtusOptimTest.cpp`
- `CMakeLists.txt`

Three registered tests simply return when slow tests are disabled, so CTest
reports them as passed even though their bodies did not execute. Both Otus tests
run normally in Java; only `TestOptimLMder` is skipped there.

Recommended fix:

- Split slow optimizer cases into a separate test executable or individually
  registered CTest cases with a `slow` label.
- Make the skipped state visible rather than reporting a pass.
- Run the two Otus parity tests in a scheduled or explicit CI job.
- Keep the Nikkor test separate until its exception failure and runtime are
  understood.

### P3: direct glass-catalog access can observe an empty catalog

Locations:

- `include/redukti/rayoptics/seq/Glass.h`
- `src/rayoptics/seq/Glass.cpp`

Java initializes the catalog before any static field access. C++ lookup methods
load it lazily, but the public `Glass::glasses()` accessor itself does not call
`ensureCatalogLoaded()`. Its existing test calls the initializer first and
therefore masks this difference.

Recommended fix: separate private catalog storage from the public accessor. The
public accessor should ensure initialization; the initialization routine should
write directly to the private storage to avoid recursion. Add an isolated test
that makes `glasses()` the first catalog call.

## Numerical differences from Java

The `sin`/`cos` explanation is supported, but it is not a complete explanation
for every numerical difference.

- On the reviewed Windows/MSVC system, the focused probe found one of 24 sampled
  `sin` results different from JDK 25 and no differing sampled `cos` results.
  The difference was within one unit in the last place.
- Hexapolar and Gaussian pupil generation call `sin` and `cos`.
- Pre-optimization ray and spot results generally agree at approximately
  `1e-12` to `1e-14`.
- The Otus contrast merit starts approximately `3e-13` relatively away from
  Java.
- Histogram bin assignment and adaptive FFT sizing are discontinuous. A
  last-bit coordinate difference can move a complete ray weight to another bin
  or alter the FFT size, increasing the derived MTF difference to approximately
  `1e-5`.
- A long Levenberg-Marquardt solve can amplify a tiny initial difference and
  follow a different path through a shallow valley, producing percent-level
  differences in the converged prescription.
- The MINPACK battery separately records path divergence involving `pow` and
  `exp`, so transcendental differences are not limited to `sin` and `cos`.

The most reliable optimizer parity checks are therefore the pre-solve merit
vector, the finite-difference Jacobian, the first few accepted iterations, and
physical post-solve metrics. An exact final prescription is not generally a
portable target across math libraries.

## Port completeness

The central prescription, sequential model, ray tracing, analysis,
plotting/rendering, optimizer, OpticalBench importer, and Zemax export paths are
present. The repository is not yet a complete Java API port.

Material areas without C++ counterparts include:

- AGF importing and its dispersion formula classes;
- Beam42 export;
- Java, Python, and rayoptics model/test writers;
- real-valued FFT variants;
- several tool/example applications and lens-model entry points.

Some apparent Java filename omissions are not missing functionality: several
small Java records/classes and test classes have deliberately been consolidated
into broader C++ headers and test files.

## Previously fixed during review

- Builder-created merit functions and solvers now retain shared ownership of
  their `Analysis`.
- Public optimizer variables and constraints validate surface, scenario, and
  aspheric coefficient indices.
- Wavelength documentation now records the deliberate case-sensitive distinction
  between sodium `D` and helium `d`, while other names retain the
  case-insensitive fallback.
