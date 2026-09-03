# Performance: C++ port vs the Java original

First measurement of the ported tool against the Java it was ported from,
taken 2026-09-03. **The initial C++ baseline was about 1.5x slower than the Java.
After two optimisations it is now slightly faster than the Java.**

This file records what was measured, the leading explanation, and the
optimisation experiments in the order they were done. The original numbers below
remain the baseline.

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

Note that the 64.17 s "before" here is higher than the 59.4-60.2 s in the
baseline table above, so the two sets are not directly comparable — machine
state differed between sittings. Re-measured later in one sitting, the reserved
build ran 44.27 / 44.21 / 43.94 s against the ~60 s baseline, about 26%. Compare
within a block of runs, not across them.

With the opt-in path timer enabled, the post-change run reported:

| | |
|---|---:|
| `SequentialModel::path()` calls | 3,026,267 |
| Cumulative time in `path()` | 17.149 s |
| Average time per call | 5.667 us |

Set `RAYOPTICS_PROFILE_PATH` to a non-empty value other than `0` when running an
executable to enable this timing. The result is written to standard error at
normal process exit. Timing is disabled by default.

## Second optimisation experiment: borrow instead of copying

`PathSeg` held `shared_ptr<Interface>`, `shared_ptr<Gap>` and an inline
`optional<Tfm3d>`. None of those objects change once the model is built, and the
`SequentialModel` that owns them outlives every trace run against it, so the
reference counting bought nothing and the transform copy was pure overhead.

`PathSeg` now borrows all three as raw pointers, and the slices that feed
`zip_longest` collect pointers rather than copying elements (`Lists::slice_ptrs`
and `Lists::slice_addrs`). That removes, per segment, two atomic refcount
operations and a 104-byte transform copy -- twice over, since the old code
copied once into the slice and again into the `PathSeg`.

| | |
|---|---:|
| `sizeof(PathSeg)` before | 168 bytes |
| `sizeof(PathSeg)` after | **48 bytes** |

Measured on the same machine and prescription, three runs each:

| | run 1 | run 2 | run 3 |
|---|---:|---:|---:|
| Before (reserve only) | 44.27 s | 44.21 s | 43.94 s |
| After (borrowed pointers) | 41.14 s | 35.84 s | 36.55 s |

Run 1 after the change is an outlier, most likely a cold file cache right after
the rebuild; the settled figure is about 36.2 s, a further **~18% reduction**.

The path profiler confirms where it went:

| | calls | cumulative | average |
|---|---:|---:|---:|
| Before | 3,026,267 | 17.149 s | 5.667 us |
| After | 3,026,267 | **8.963 s** | **2.962 us** |

Time in `path()` roughly halved, and the 8.2 s saved there accounts for
essentially all of the 7.9 s saved overall — the two agree, which is the check
that the attribution is right.

### The lifetime rule this introduces

Borrowing is only safe because the model is immutable while a path is in use. A
`PathSeg` is now valid only until the model is rebuilt, so it must not be held
across `update_model`, `add_surface`, or anything else touching `ifcs`, `gaps`
or `lcl_tfrms`. `reverse_path()` computes its transforms on the fly, so they are
parked in a member (`reverse_path_tfrms_`) rather than a local that would dangle
the moment it returned; its single caller consumes the path immediately.

The full test suite passes unchanged, including the end-to-end `LensTool2` run
that compares all 29 artifacts byte-for-byte against the committed Examples.

## Where it stands now

Three runs each, same machine, same sitting:

| | run 1 | run 2 | run 3 |
|---|---:|---:|---:|
| C++ Release (both optimisations) | 41.14 s | 35.84 s | 36.55 s |
| Java (JDK 25) | 39.31 s | 38.08 s | 40.34 s |

From ~60 s to ~36 s is a **40% reduction against the original baseline**, and
puts the C++ slightly ahead of the Java rather than 1.5x behind it.

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

The original implementation did still more work during vector growth, which the
first experiment removed, and copied every interface, gap and transform, which
the second removed. What remains is that every ray still rebuilds an identical,
immutable path -- now as 3.03M allocations of a 48-byte-per-segment vector
rather than 168.

The root cause was a direct consequence of the faithful port. In Java,
`PathSeg.Tfrm` is a *reference* — 8 bytes copied. Here it became an inline
`std::optional<Tfm3d>` — 104 bytes copied. Java moved ~40 bytes per segment with
no atomics; this moved 168 with two. Multiplied by ~112M segment constructions,
that was the shape of the gap. The second experiment restored the reference
semantics and closed it.

Java also gets whole-program inlining from the JIT at runtime, which this build
does not have.

## Further candidate fixes, in rough order of expected value

**None of these further changes have been measured.** They are sized by
inspection only.

1. **Cache the path instead of rebuilding it per ray.** `RayTrace::trace` calls
   `seq_model->path(...)` on every trace, 3.03M times, and still spends ~9 s
   there. The path is immutable between model updates, so it could be computed
   once per (wavelength, model version). Now the largest remaining item, but it
   requires reliable invalidation whenever model geometry or refractive-index
   data changes — a correctness risk the two experiments above did not carry.
2. **Enable `/GL /LTCG`** in the Release configuration for cross-translation-unit
   inlining, which the JIT gets for free.
3. **Revisit exceptions.** Still unmeasured, and now a larger share of a smaller
   total.

Items 2 and 3 of the original list — borrowing the transform, and dropping the
shared_ptr refcounting — were done; see the second experiment above.

## Caveat on the comparison

This measures one lens on one machine. The workload is dominated by ray tracing
through the hexapolar spot analysis (`SpotOptions` defaults to 64 rings, and the
tool runs two spot analyses per configuration, twice over). A different lens or
a different analysis mix could shift the balance.
