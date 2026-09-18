# solver

A compact C++ finite-element solver focused on clear element formulations and an extensible solver architecture.

## v0.1

- 3D two-node Euler-Bernoulli beam/frame element (`Beam3D`)
- 6 DOFs/node: `ux, uy, uz, rx, ry, rz`
- axial, torsional and biaxial bending stiffness
- robust local coordinate frame and 12x12 local/global transformation
- deterministic DOF numbering and global stiffness assembly
- nodal loads, fixed supports, linear static solution and reaction recovery
- regression tests against analytical cantilever solutions and a rotated member

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Eigen 3.4 is used for linear algebra. CMake uses a system Eigen when available and otherwise fetches it.

See `examples/cantilever.cpp` and `docs/architecture.md`.
