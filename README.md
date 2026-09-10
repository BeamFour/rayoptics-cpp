# RayOptics-Cpp

This is a C++ port of [Beam42](https://github.com/BeamFour/Beam42) and Michael Hayford's [RayOptics](https://github.com/mjhoptics/ray-optics).
The functionality of this port is almost identical to Beam42's RayOptics component. For full details please refer
to that project. Differences are noted below.

## Executables

* LensTool2 as documented in [LensTool2](https://github.com/BeamFour/Beam42/blob/main/Documentation/LENSTOOL2.md). See exceptions below.

## Differences from Beam42

* The LensTool2 utility does not support obtaining prescriptions directly from [PhotonsToPhotos Optical Bench](https://www.photonstophotos.net/GeneralTopics/Lenses/OpticalBench/OpticalBenchHub.htm).
* The Optimizer is included but not yet documented.
* Michael Lampton's BeamFour product is not included.
* Some utilities used during development and verification are not included.

## Porting Strategy

* Beam42 remains the primary development project. Development occurs there first and is then ported here.
* The C++ port was created using Claude and Codex and aims to be a faithful replica of the Java version, except for details such as memory management.
* The performance of both projects is comparable. See [Performance](PERFORMANCE.md) doc for details. Extreme performance is not a goal for either project.
  The portions that are derived from Michael Hayford's RayOptics aim to maintain the same overall structure as the original Python version. This
  is to ease verification and maintainability.
* I tried GraalVM Community Edition as a native code generator for the Java version - unfortunately it produces executables that run 5x slower than the JVM.

## License

The project includes code derived from several opensource projects. See the individual license notices in the source code and in LICENSE notices.
The overall license is GNU GPL v3 or later; see [LICENSE-GPL-3.0.txt](LICENSE-GPL-3.0.txt).