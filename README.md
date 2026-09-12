# RayOptics-Cpp

This is a C++ port of [Beam42](https://github.com/BeamFour/Beam42) and Michael Hayford's [RayOptics](https://github.com/mjhoptics/ray-optics).
The functionality of this port is almost identical to Beam42's RayOptics component. For full details please refer
to that project. Differences are noted below.

## Executables

* LensTool2 as documented in [LensTool2](https://github.com/BeamFour/Beam42/blob/main/Documentation/LENSTOOL2.md). See exceptions below.
* An optimization run is driven from `[trial n]` and `[pipeline n]` sections of the prescription, as documented in
  [Optimizer](https://github.com/BeamFour/Beam42/blob/main/Documentation/OPTIMIZER.md); `--optimize n` runs one.

## Building

The build uses CMake and a C++20 compiler, and needs no third-party dependencies beyond the bundled Ryu.

```bash
cmake -S . -B build
cmake --build build --config Release
ctest --test-dir build -C Release
```

Single-configuration generators (Makefiles, Ninja) default to a Release build. Multi-configuration generators
such as Visual Studio need the configuration given explicitly, as above.

The tests read the prescriptions under `Examples/`. The optimizer solve tests are slow and off by default;
set `RAYOPTICS_RUN_SLOW_TESTS=1` to run them.

## Differences from Beam42

* The LensTool2 utility does not support obtaining prescriptions directly from [PhotonsToPhotos Optical Bench](https://www.photonstophotos.net/GeneralTopics/Lenses/OpticalBench/OpticalBenchHub.htm).
* Michael Lampton's BeamFour product is not included.
* Some utilities used during development and verification are not included.

## Porting Strategy

* Beam42 remains the primary development project. Development occurs there first and is then ported here.
* The C++ port was created using Claude and Codex and aims to be a faithful replica of the Java version, except for details such as memory management.
* The performance of both projects is comparable. See [Performance](PERFORMANCE.md) doc for details. Extreme performance is not a goal for either project.

## License

The project includes code derived from several open-source projects. See the individual license notices in the source code
and the LICENSE files:

* [LICENSE-GPL-3.0.txt](LICENSE-GPL-3.0.txt) - the overall license, and the license of the code derived from Goptical.
* [LICENSE-ray-optics.txt](LICENSE-ray-optics.txt) - BSD 3-Clause, for the code derived from Michael Hayford's RayOptics.
* [LICENSE-Minpack.txt](LICENSE-Minpack.txt) - for the code derived from MINPACK.
* [LICENSE-ryu.txt](LICENSE-ryu.txt) - for the bundled Ryu, used for number formatting (see also `third_party/ryu`).

The overall license is GNU GPL v3 or later.