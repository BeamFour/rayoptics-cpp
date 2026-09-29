# Tools

Scripts used to maintain the port. None is part of the build: the CMake build lists its
sources explicitly and ignores this directory. Python 3 and, for the benchmark, Windows
PowerShell.

Several of them reach into a Beam43 checkout, which is this port's ground truth. They
default to `../Beam43`, beside this repository, and take an argument when it is elsewhere.

| Tool | Does |
| --- | --- |
| `regen_expected.py` | Rewrites one `EXPECTED_*` array in a `tests/*Expected.h` from a Java dump |
| `regen_raw.py` | The same for a block written as one raw string, `R"JV(...)JV"` |
| `sync_upstream_tests.py` | Carries regenerated Java upstream-test expectations into `tests/upstream/*.cpp` |
| `licence_audit.py` | Compares licence headers here with the Java headers they port. Read-only |
| `benchmark.ps1` | Times this build against the Java, as [PERFORMANCE.md](../PERFORMANCE.md) describes |

## Regenerating expected values

The port is verified by asserting values dumped from the Java, so when Java output
changes deliberately, the C++ expectation is regenerated rather than edited. The dump
programs live on the Java side, in `rayoptics/src/test/java/org/redukti/cppport/`, one
per C++ test that carries a large expectation:

| Java program | Feeds | Applied with |
| --- | --- | --- |
| `DumpSpec` | `tests/SpecExpected.h` | `regen_expected.py` |
| `DumpLayout` | `tests/LayoutPlotterExpected.h`, one file per array | `regen_expected.py` |
| `DumpAnalysis` | `tests/AnalysisExpected.h`, one file per block | `regen_raw.py` |

Which of the two appliers to use depends on how the header stores the block: an array of
one string per line, or a single raw string. The table says which.

Run one from the Beam43 checkout, having compiled it (`mvn -o test-compile`):

```
java -cp "rayoptics/target/classes;rayoptics/target/test-classes" \
    org.redukti.cppport.DumpSpec spec-dump.txt
```

Then apply it here, with `--dry` first to see what would change:

```
python tools/regen_expected.py tests/SpecExpected.h EXPECTED_SPEC_LINES spec-dump.txt --dry
python tools/regen_expected.py tests/SpecExpected.h EXPECTED_SPEC_LINES spec-dump.txt
```

The dump programs read their prescriptions from this checkout's `Examples`, since those
are the copies the tests assert against.

A regeneration on unchanged Java reports `lines differing=0`, since these expectations are
the Java's own output. Anything else means either the Java output really moved, or the
options in the dump program have drifted from the C++ test's -- check which before writing.

The `sin`/`cos` divergence between the JVM and MSVC that `PERFORMANCE.md` notes appears
when the *C++* is run against these values, not here.

## Upstream tests

`tests/upstream/*.cpp` mirror the Java tests of the same name, which assert against
upstream ray-optics. When the Java expectations are regenerated, pass the revision that
holds the values the C++ currently matches -- the commit *before* that regeneration:

```
python tools/sync_upstream_tests.py --old-rev abc1234^ --dry
```

Each assertion is matched by its label, and a value is replaced only when the C++ still
agrees with the old Java. Anything else is reported and nothing is written, so a test
whose meaning moved gets looked at by hand.

Otherwise it behaves like the Java regeneration: every tolerance becomes the Java's.
A C++-only widening -- where the C++ misses the upstream value by more than the Java does,
typically the JVM/MSVC `sin`/`cos` difference amplified -- is reset to the Java tolerance,
listed in the output, and its "Manually widened" comment removed. Build and run the
tests afterwards and re-widen, with a comment, any assertion that then fails, as is done
on the Java side. The Java's own widened tolerances and their comments are carried across.

## Licence headers

Attributions here must match the Java files they port; Beam43's `AGENTS.md` section 4 has
the rules and the categories.

```
python tools/licence_audit.py            # only the disagreements
python tools/licence_audit.py --all      # every file
```

It reports a file with no Java counterpart as well. That is usually a C++-only file, which
is fine, so read the output rather than acting on the count.

## Benchmark

```
powershell -File tools\benchmark.ps1                   # plain run, 3 pairs
powershell -File tools\benchmark.ps1 -Workload flags -Profile
```

It alternates C++ and Java runs one process at a time and reports the ratio, because
absolute times on a laptop are not comparable between sittings. Close editors and
browsers first. Its work directory is `cmake-build-perf/`, deliberately outside `%TEMP%`,
where Defender's scanning of newly written files inflates the times.
