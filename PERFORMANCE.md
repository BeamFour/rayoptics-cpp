# Performance: C++ port vs the Java original

First measurement of the ported tool against the Java it was ported from,
taken 2026-09-03. **The initial C++ baseline was about 1.5x slower than the Java.
After caching the path it runs at about 0.73x the Java.**

This file records what was measured, the leading explanation, and the first
optimisation experiment. The original numbers below remain the baseline.

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
measured at about 18% and then reverted. It is preserved in commit 99345790,
reverted by 71229d20. The work below replaces it and does not depend on it.

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

## Caveat on the comparison

This measures one lens on one machine. The workload is dominated by ray tracing
through the hexapolar spot analysis (`SpotOptions` defaults to 64 rings, and the
tool runs two spot analyses per configuration, twice over). A different lens or
a different analysis mix could shift the balance.
