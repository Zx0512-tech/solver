# solver

A compact C++ finite-element solver focused on clear element formulations and extensible analysis architecture.

## Implemented

### Structural model
- nodes, constraints and nodal loads
- linear-elastic material and beam section
- deterministic global DOF numbering
- dense global stiffness, mass and load assembly
- element-load abstraction with element-owned equivalent nodal loading

### Beam3D
- 2-node 3D Euler-Bernoulli beam/frame element
- 6 DOFs per node: `ux, uy, uz, rx, ry, rz`
- axial, torsional and biaxial bending stiffness
- consistent 12x12 mass matrix including torsional rotary inertia
- robust local coordinate system
- local/global transformation for stiffness, mass, loads and forces
- local deformation and end-force response recovery

### Loads and response
- nodal forces and moments
- uniform Beam3D line load in local x/y/z directions
- consistent equivalent nodal loads
- element resisting/end forces include fixed-end load effects
- static element response recording
- transient element-force history recording

For a Beam3D, the local end-force order is:

`[N_i, Vy_i, Vz_i, T_i, My_i, Mz_i, N_j, Vy_j, Vz_j, T_j, My_j, Mz_j]`

### Analysis
- linear static solution and reaction recovery
- generalized eigenvalue modal analysis `K phi = omega^2 M phi`
- Rayleigh damping `C = alpha M + beta K`
- Newmark-beta transient integration for `M a + C v + K u = F(t)`
- arbitrary nodal time-history loads
- initial displacement and velocity support

## Response architecture

```text
ElementLoad
   |
   +--> equivalent element nodal load --> global F
   |
Solver --> global displacement U
                     |
                     v
                 Element
                     |
             ElementResponse
          local deformation / forces
                     |
                     v
              ElementRecorder
           static or time history
```

The element, not the recorder, owns the mechanics needed to recover its response.

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
