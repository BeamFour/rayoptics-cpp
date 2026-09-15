# Performance: C++ port vs the Java original

First measurement of the ported tool against the Java it was ported from,
taken 2026-09-03. **The initial C++ baseline was about 1.5x slower than the Java.
After caching the path it ran at about 0.73x the Java.**

**Re-measured on 2026-09-14, after the path cache and list presizing were
backported to the Java, the C++ runs at about 0.95x the Java: the two are
now comparable.** See [Re-measurement, 2026-09-14](#re-measurement-2026-09-14).

This file records what was measured, the leading explanation, and the
optimisation experiments in the order they were made. The original numbers
below remain the baseline.

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
copy. Without `--outdir`, `LensTool2` writes every artifact next to the spec
file, so running it directly against a folder under `Examples/` overwrites
committed files. `LensTool2Test` copies the spec to a scratch directory for
this reason; do the same by hand. (At the time of the first measurements the
`.zmx` was written next to the spec even when `--outdir` was given; it now
follows `--outdir` like everything else.)

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

## First optimisation experiment: reserve vector capacity

`SequentialModel::path()` creates five sliced vectors, a selected refractive-index
vector, and the final `PathSeg` vector. None reserved capacity, so vector growth
repeatedly copied transforms and `PathSeg` values and repeatedly adjusted
`shared_ptr` reference counts.

Capacity is now reserved for these vectors. On the same machine and prescription,
one before/after pair measured:

| | wall clock |
|---|---:|
| Before | 64.17 s |
| After, profiling disabled | 45.08 s |
| After, path profiling enabled | 45.58 s |

This is a **29.7% reduction** for the unprofiled pair. It is one pair rather than
a full benchmark series, so the exact percentage should not be over-interpreted,
but the size of the change confirms that path materialisation was a major cost.

With the opt-in path timer enabled, the post-change run reported:

| | |
|---|---:|
| `SequentialModel::path()` calls | 3,026,267 |
| Cumulative time in `path()` | 17.149 s |
| Average time per call | 5.667 us |

Set `RAYOPTICS_PROFILE_PATH` to a non-empty value other than `0` when running an
executable to enable this timing. The result is written to standard error at
normal process exit. Timing is disabled by default.

## Assessment

**Exceptions are unlikely to explain the entire gap, but remain unmeasured.**
The port uses exceptions as control flow, so they were the first suspect. At an
assumed 50 microseconds per throw on Windows x64, 118k throws would be about 6
seconds. That estimate is not a substitute for measuring exception-heavy and
exception-free workloads separately.

**Copying the ray path is a confirmed major cost.** Every ray trace rebuilds the
whole path:

- ~37 interfaces, so each `path()` call builds a ~37-element `vector<PathSeg>`,
  about 6 KB
- two `shared_ptr` atomic increments per element, 74 per call
- over 3.03M calls: at least roughly **19 GB of final `PathSeg` construction and
  ~450M atomic refcount operations**

The original implementation did still more work during vector growth; the first
optimisation experiment above removed that portion. Most calls nevertheless
continue to rebuild an identical, immutable path for each ray.

The root cause is a direct consequence of the faithful port. In Java,
`PathSeg.Tfrm` is a *reference* — 8 bytes copied. Here it became an inline
`std::optional<Tfm3d>` — 104 bytes copied. Java moves ~40 bytes per segment with
no atomics; this moves 168 with two. Multiplied by ~112M segment constructions,
that is the shape of the gap.

Java also gets whole-program inlining from the JIT at runtime, which this build
does not have.

## Second optimisation experiment: bounded cache of computed paths

A separate experiment -- borrowing the interfaces, gaps and transforms as raw
pointers instead of holding `shared_ptr`s and copying the transform -- was
measured at about 18% and then reverted. It is preserved in commit 262d1f69,
reverted by a07c0201. The work below replaces it and does not depend on it.

`RayTrace::trace` calls `path()` once per ray, and every call rebuilt an
identical vector. The path depends only on the wavelength and the surface range
asked for, and cannot change while the model is unchanged, so it can be built
once and handed back.

`SequentialModel` now holds a ring of at most `PATH_CACHE_CAPACITY` (16) entries
keyed on the four arguments as given. A miss builds the path and takes the next
slot, overwriting the oldest once the ring is full. `path()` returns a
`const std::vector<PathSeg> &` into that ring rather than a fresh vector, which
is what removes the per-ray copy; `RayTrace::trace` binds it by reference.

### Invalidation

A `PathSeg` copies the transform, refractive index and z direction out of the
model, so a cached entry goes stale when `lcl_tfrms`, `rndx`, `z_dir` or
`wvlns` is rewritten, or when the interface and gap lists are spliced. It does
*not* go stale when an `Interface` or `Gap` object is edited in place, because
the entry holds `shared_ptr`s to those same objects -- that is what keeping the
`shared_ptr`s buys.

Three methods write those arrays, and nothing outside the class touches them, so
`invalidate_path_cache()` is called at the end of each: `initialize_arrays()`,
`insert()` (and so `add_surface`), and `update_model()`. Invalidating last
rather than first means anything reached during the update itself cannot leave a
stale entry behind.

The ring is reserved to full capacity before the first insert. That is for
correctness, not speed: growing the vector must never reallocate, or a reference
handed out earlier would dangle.

### Results

The decisive measurement is the internal profile, which is immune to the machine
drift discussed below:

| | calls | cache hits | hit rate | cumulative | average |
|---|---:|---:|---:|---:|---:|
| Reserve only | 3,026,267 | -- | -- | 17.149 s | 5.667 us |
| With the cache | 3,032,053 | 3,032,005 | **99.998%** | **0.099 s** | **0.033 us** |

48 misses in three million calls. The ring never thrashes -- the misses are
essentially one per distinct key per model rebuild -- so 16 slots is ample, and
the linear scan over them costs nothing measurable.

Wall clock, before and after measured back to back in one sitting:

| | run 1 | run 2 | run 3 |
|---|---:|---:|---:|
| Before (reserve only) | 58.50 s | 50.42 s | 50.52 s |
| After (path cache) | 27.96 s | 27.77 s | 27.42 s |

Run 1 of the "before" set is a cold-cache outlier. Against the settled ~50.5 s
that is a **~45% reduction**.

### A caveat on the wall-clock numbers

Absolute times drifted upward substantially over this session: Java, unchanged
throughout, ran 39-40 s early on and 48-54 s an hour later. Any before/after
pair spanning that drift is meaningless. Only the back-to-back pair above, and
the ratio below, should be read as measurements.

Interleaved C++/Java pairs, run alternately so drift affects both equally:

| | C++ | Java | ratio |
|---|---:|---:|---:|
| pair 1 | 33.90 s | 48.61 s | 0.70 |
| pair 2 | 38.92 s | 53.80 s | 0.72 |
| pair 3 | 38.06 s | 48.17 s | 0.79 |

The C++ runs at roughly **0.73x the Java**, consistently, whatever the absolute
level. Against the original 1.5x-slower baseline, that is the whole gap closed
and then some.

The full test suite passes unchanged in Release and Debug -- 168 tests -- and
that includes the end-to-end `LensTool2` run comparing all 29 artifacts against
the committed Examples.

### Cache implementation caveats

The cache is designed for the current single-threaded use of a
`SequentialModel`. It is not thread-safe: concurrent cache misses, invalidation,
or tracing while the model changes would race and could invalidate a returned
reference. Keep a model confined to one thread unless the cache is later given
explicit synchronization and reference-lifetime rules.

The arrays from which a path is copied are currently public. No code in this
repository writes them outside `SequentialModel`, but the type system does not
enforce that rule. Any future direct mutation of `lcl_tfrms`, `rndx`, `z_dir`,
`wvlns`, `ifcs`, or `gaps` must also invalidate the cache.

The hit and miss counters currently perform a relaxed atomic increment on every
lookup even when path profiling is disabled. This makes the reported hit rate
available, but leaves a small profiling cost in normal execution. Gating those
counters behind `RAYOPTICS_PROFILE_PATH` is a separate experiment so it is not
mixed into other measurements.

## Third optimisation experiment: reserve traced-ray capacity

`RayTrace::trace_raw()` builds a fresh `std::vector<RaySeg>` for every ray. A
completed trace emits at most one segment per path entry, while blocked,
filtered, missed-surface and total-internal-reflection paths emit fewer. The
implementation now reserves `path.size()` entries before tracing, which is a
safe upper bound and prevents repeated vector growth in the normal case.

This experiment deliberately changes only the fundamental ray vector. Other
known opportunities are being kept separate so their effects remain
attributable:

- the three exception handlers copy the traced segment vector into `RayPkg`,
  even though the local vector is not used after the exception is rethrown;
- fan results can reserve `num_rays`;
- rectangular-grid results can reserve `num_rays * num_rays`;
- ring-grid and contrast results can reserve their generated point count;
- the hexapolar and ordinary ring point generators can reserve their exact or
  upper-bound output sizes.

One pre-change run and two post-change runs, taken in the same session, measured:

| | run 1 | run 2 |
|---|---:|---:|
| Before | 35.06 s | -- |
| After reserving `RaySeg` capacity | 23.31 s | 23.41 s |

Against the mean post-change time of 23.36 s, this is a **33.4% reduction**.
Only one pre-change run was captured, so a longer interleaved series would be
needed for a precise percentage. The two post-change runs agree closely and the
effect is much larger than the short-term run-to-run variation observed here.

## Re-measurement, 2026-09-14

Since the experiments above, two of the C++ changes were backported to the
Java (Beam43 commits 26e5ae31, list presizing including the `RaySeg` list, and
efaf5c84, the path cache), so the 0.73x ratio no longer described the pair.

### Method

As above -- same lens, same machine, a scratch copy of the prescription -- with
these differences:

- Browsers and IDEs were closed, and nothing else was run during a series.
- Runs are strictly one process at a time, C++ and Java alternating, so that
  drift affects both equally.
- Two workloads: the plain run, which writes the same 29 artifacts as before,
  and a run with `--output-ray-aberration-plots --output-wavelength-mtfs
  --do-wideangle-layout`, which writes 227.

| | |
|---|---|
| C++ | rayoptics-cpp dd8e10b3, MSVC 19.40.33811, Release (`/O2 /Ob2 /DNDEBUG`), no LTCG |
| Java | Beam43 9d11c8bd classes, OpenJDK Temurin 25.0.1+8 LTS, default JVM options |

### Results

Only ratios are recorded. The absolute times in this sitting are suspect: every
run wrote its artifacts under `%TEMP%`, and Windows Defender real-time scanning
was active, which may add a cost unrelated to the code. Because the two versions
alternate and write the same files, that cost applies to both and the ratio
still compares them.

C++ time as a fraction of the Java's:

| | plain run, 29 artifacts | with the output flags, 227 artifacts |
|---|---:|---:|
| pair 1 | 0.97 | 0.99 |
| pair 2 | 0.96 | 0.93 |
| pair 3 | 0.94 | 1.02 |
| ratio of means | **0.95** | **0.98** |

The extra plots add under 10% to either version's time; the work is in the spot
and MTF analyses both workloads share.

The path cache still behaves as measured above:

| | calls | cache misses | hit rate |
|---|---:|---:|---:|
| Plain run | 3,117,809 | 48 | 99.998% |
| With the output flags | 3,123,595 | 48 | 99.998% |

With the output flags, 185 of the 227 artifacts are byte-identical between the
two versions (ignoring line endings). The 42 that differ are all MTF SVGs and
the spot reports -- the documented hexapolar sin/cos divergence, now also
reaching 36 of the 110 per-field MTF plots.

The full C++ test suite, now 300 tests, passes.

### The C++ did not get slower

The absolute C++ times in this sitting were much longer than those recorded
after the third experiment. To tell a regression from a change in the
environment, the 2026-09-04 tree (84b1b439, which already has the path cache
and the `RaySeg` reservation) was built with the same compiler and flags and
timed against the current build, alternating:

| | current time / 84b1b439 time |
|---|---:|
| pair 1 | 0.99 |
| pair 2 | 0.97 |
| pair 3 | 1.02 |
| ratio of means | **0.99** |

The two builds are indistinguishable. The old build makes 3,026,267 `path()`
calls, as recorded above; the current one makes about 3% more, from the
separate 21-ring spot-diagram analysis. The longer absolute times therefore
come from the environment, not the code; Defender scanning the output files is
the leading suspect but was not isolated. That confirms the caveat on
wall-clock numbers: **compare only runs made back to back in one sitting**, and
read ratios rather than absolute times across sessions.

The change in ratio, from 0.73x to 0.95x, therefore comes from the Java side.
The backported path cache and presizing are the likely cause, but the Java
before and after them was not timed in this sitting, so their individual
effects are not measured here.

## Rejected experiment: pocketfft in place of fftpack

On 2026-09-11 pocketfft (header-only C++, BSD-3) was trialled as a drop-in for
`ComplexDoubleFFT` in `BaseMTF::compute_fft`.

| | |
|---|---|
| FFT alone, at MTF sizes (1024-4096) | about 2x faster |
| Round-trip accuracy | about 20x better |
| Agreement with fftpack | to about 7e-16; 2.8% of doubles bit-identical |
| Tests and the 111 `LensTool2` outputs across 5 Examples | all passed; outputs byte-identical |
| End-to-end `LensTool2` time | no measurable change |

The FFT is a negligible share of a run. The trial was reverted: it would add a
dependency, diverge from the Java's fftpack, and take MTF out of the bit-exact
comparison against the JVM that `FftTest` provides.

## Further candidate fixes, in rough order of expected value

**None of these further changes have been measured.** They are sized by
inspection only.

1. **Enable `/GL /LTCG`** in the Release configuration for cross-translation-unit
   inlining, which the JIT gets for free.
2. **Revisit exceptions.** Still unmeasured, and now a larger share of a much
   smaller total: 117,857 throws against a run that no longer spends meaningful
   time building paths.
3. **Reconsider borrowing in `PathSeg`** (the reverted experiment). With
   `path()` down to 0.099 s cumulative, it would now be optimising something
   that has stopped mattering. Its value was in what it revealed about where the
   time went, not in what is left to win.

Items 1 and 2 of the original list -- caching the path, and holding
`const Tfm3d *` in `PathSeg` -- are done and reverted-as-superseded
respectively; see the experiment above.

The smaller follow-ups listed under the third experiment -- the exception
handlers' copy of the traced segments, and reserving capacity for fan, grid,
ring and contrast results and the point generators -- are also still open as of
2026-09-14. With the C++ and Java now comparable, none is needed to meet the
project's goal; they remain candidates only.

## Caveat on the comparison

This measures one lens on one machine. The workload is dominated by ray tracing
through the hexapolar spot analyses. For each configuration the tool now runs
one 21-ring analysis for the spot diagrams and three 64-ring analyses: the spot
report, the MTF and the wavelength MTFs. A different lens, a different
`--spot-pattern`, or a different analysis mix could shift the balance.
