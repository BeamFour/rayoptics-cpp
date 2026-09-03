# Performance: C++ port vs the Java original

First measurement of the ported tool against the Java it was ported from,
taken 2026-09-03. **The C++ is currently about 1.5x slower than the Java.**

This file records what was measured and the leading explanation. Nothing here
has been acted on — no optimisation work has been done, and the numbers below
are the baseline to beat if it ever is.

## Method

Both versions run `LensTool2` over the same prescription and write the same 29
report artifacts. Wall clock, three consecutive runs each, nothing else running.

| | |
|---|---|
| Lens | `Examples/jfotoptix/canon-rf70-200mm-f2.8LZ/US20250155694_Example01P.txt` |
| | a two-configuration zoom, ~35 surfaces |
| Machine | AMD Ryzen 7 PRO 4750U, 16 threads, Windows 11 |
| C++ | MSVC, Visual Studio 17 2022, Release (`/O2 /Ob2 /DNDEBUG`), no LTCG |
| Java | OpenJDK Temurin 25.0.1+8 LTS, default JVM options |

Both were run against a scratch copy of the prescription, not the repository
copy. **`LensTool2` writes its `.zmx` next to the spec file regardless of
`--outdir`**, so running it directly against a folder under `Examples/`
overwrites committed files. `LensTool2Test` copies the spec to a scratch
directory for this reason; do the same by hand.

## Results

| | run 1 | run 2 | run 3 |
|---|---|---|---|
| C++ Release | 59.4 s | 60.1 s | 60.2 s |
| Java (JDK 25) | 38.3 s | 40.8 s | 40.9 s |

The Java figures include JVM startup and JIT warm-up. Both produce byte-identical
output apart from the six artifacts affected by the documented hexapolar
sin/cos divergence, and the README generation date.

## What was measured inside the C++ run

Temporary instrumentation, since reverted:

| | |
|---|---|
| `SequentialModel::path()` calls | **3,026,267** — one per ray trace |
| `sizeof(PathSeg)` | **168 bytes** |
| `sizeof(Tfm3d)` | **104 bytes** |
| `TraceException` constructions | **117,857** |

## Assessment

**Exceptions are not the problem.** The port uses exceptions as control flow, so
they were the first suspect. At a pessimistic 50 microseconds per throw on
Windows x64, 118k throws is about 6 seconds — at most a tenth of the 20 second
gap, and probably far less. Ruled out.

**Copying the ray path is the leading candidate.** Every ray trace rebuilds the
whole path:

- ~37 interfaces, so each `path()` call builds a ~37-element `vector<PathSeg>`,
  about 6 KB
- two `shared_ptr` atomic increments per element, 74 per call
- over 3.03M calls: roughly **19 GB of construction and copying, and ~450M
  atomic refcount operations**

All of it rebuilding an identical, immutable path for each ray.

The root cause is a direct consequence of the faithful port. In Java,
`PathSeg.Tfrm` is a *reference* — 8 bytes copied. Here it became an inline
`std::optional<Tfm3d>` — 104 bytes copied. Java moves ~40 bytes per segment with
no atomics; this moves 168 with two. Multiplied by ~112M segment constructions,
that is the shape of the gap.

Java also gets whole-program inlining from the JIT at runtime, which this build
does not have.

## Candidate fixes, in rough order of expected value

**None of these have been measured.** They are sized by inspection only.

1. **Cache the path instead of rebuilding it per ray.** `RayTrace::trace` calls
   `seq_model->path(...)` on every trace. The path is immutable between model
   updates, so it can be computed once per (wavelength, model version). Likely
   the single biggest win and semantically transparent, though it departs from
   the Java call structure.
2. **Hold `const Tfm3d *` in `PathSeg`** rather than a 104-byte copy, pointing
   into the `SequentialModel`'s `lcl_tfrms`. Restores Java's reference semantics
   for the field that dominates the struct.
3. **Use raw `Interface *` / `Gap *` in `PathSeg`.** The `SequentialModel` owns
   them and outlives every trace, so the atomic refcounting buys nothing.
4. **Enable `/GL /LTCG`** in the Release configuration for cross-translation-unit
   inlining.

Options 2 and 3 change a type the whole raytr layer reads, so they need the
full test suite as cover — which is the point of having it.

## Caveat on the comparison

This measures one lens on one machine. The workload is dominated by ray tracing
through the hexapolar spot analysis (`SpotOptions` defaults to 64 rings, and the
tool runs two spot analyses per configuration, twice over). A different lens or
a different analysis mix could shift the balance.
