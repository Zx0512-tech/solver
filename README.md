# solver

A compact C++ finite-element solver focused on clear element formulations and extensible analysis architecture.

## Implemented

### Structural model
- nodes, constraints and nodal loads
- linear-elastic material and beam section
- deterministic global DOF numbering
- dense global assembly

### Beam3D
- 2-node 3D Euler-Bernoulli beam/frame element
- 6 DOFs per node: `ux, uy, uz, rx, ry, rz`
- axial, torsional and biaxial bending stiffness
- consistent 12x12 mass matrix including torsional rotary inertia
- robust local coordinate system
- local/global transformation for both stiffness and mass

### Analysis
- linear static solution and reaction recovery
- generalized eigenvalue modal analysis `K phi = omega^2 M phi`
- Rayleigh damping `C = alpha M + beta K`
- Newmark-beta transient integration for `M a + C v + K u = F(t)`
- arbitrary nodal time-history loads
- initial displacement and velocity support

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Eigen 3.4 is used for linear algebra.

Examples:
- `examples/cantilever.cpp`: linear static beam
- `examples/dynamic_cantilever.cpp`: modal + Rayleigh + Newmark transient analysis
