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
- nodal forces and nodal moments
- `BeamUniformLoad3D`: uniform line load
- `BeamLinearLoad3D`: triangular/trapezoidal line load with independent i/j intensities
- `BeamPointLoad3D`: concentrated force and/or moment at an arbitrary distance from node i
- `BeamPartialLinearLoad3D`: uniform/triangular/trapezoidal line load acting only on a selected beam interval
- `TimeDependentElementLoad`: reusable scalar time history for any spatial element load, supporting `P(t)` and `q(x,t)`
- beam loads may be defined in local or global coordinates
- consistent equivalent nodal loads based on beam interpolation functions
- element resisting/end forces retain fixed-end effects from element loads
- static element response recording
- transient element-force history recording
- section-force recovery at arbitrary beam positions: `N(x), Vy(x), Vz(x), T(x), My(x), Mz(x)`
- `BeamResponseSampler` for uniform sampling plus automatic load breakpoints
- concentrated-load locations preserve left/right limits so force jumps are not smoothed out
- CSV export with `x, side, N, Vy, Vz, T, My, Mz`
- dependency-free SVG diagrams for axial force, shear, bending moment and torsion

For a Beam3D, the local end-force order is:

`[N_i, Vy_i, Vz_i, T_i, My_i, Mz_i, N_j, Vy_j, Vz_j, T_j, My_j, Mz_j]`

Examples:

```cpp
// Uniform load in local axes.
model.addElementLoad<fem::BeamUniformLoad3D>(
    1, Eigen::Vector3d(0.0, -2000.0, 0.0));

// Global vertical load on an arbitrarily oriented beam.
model.addElementLoad<fem::BeamUniformLoad3D>(
    1,
    Eigen::Vector3d(0.0, 0.0, -5000.0),
    fem::BeamLoadCoordinateSystem::Global);

// Trapezoidal load from qi to qj.
model.addElementLoad<fem::BeamLinearLoad3D>(
    1,
    Eigen::Vector3d(0.0, -1000.0, 0.0),
    Eigen::Vector3d(0.0, -3000.0, 0.0));

// Concentrated force 1.5 m from node i.
model.addElementLoad<fem::BeamPointLoad3D>(
    1, 1.5, Eigen::Vector3d(0.0, -10000.0, 0.0));

// Partial trapezoidal load acting only from x=1.0 m to x=4.0 m.
model.addElementLoad<fem::BeamPartialLinearLoad3D>(
    1,
    1.0,
    4.0,
    Eigen::Vector3d(0.0, -1000.0, 0.0),
    Eigen::Vector3d(0.0, -3000.0, 0.0));

// Time-dependent point load P(t) = P0 sin(omega t).
model.addTimeDependentElementLoad<fem::BeamPointLoad3D>(
    [omega](double t) { return std::sin(omega * t); },
    1,
    1.5,
    Eigen::Vector3d(0.0, -10000.0, 0.0));

// Time-dependent partial distributed load q(x,t) = q0(x) * scale(t).
model.addTimeDependentElementLoad<fem::BeamPartialLinearLoad3D>(
    [](double t) { return 1.0 + 0.5 * std::sin(20.0 * t); },
    1,
    1.0,
    4.0,
    Eigen::Vector3d(0.0, -1000.0, 0.0),
    Eigen::Vector3d(0.0, -3000.0, 0.0));
```

### Beam force diagrams and CSV output

```cpp
const auto result = fem::LinearStaticSolver{}.solve(model);

const auto samples =
    fem::BeamResponseSampler{}.sample(
        model,
        beam_id,
        result.displacement,
        31);

fem::BeamForceCsvWriter{}.writeFile(
    "results/beam_1_forces.csv",
    samples);

fem::BeamForceDiagramSvgWriter{}.writeSet(
    "results",
    "beam_1",
    samples);
```

The diagram writer creates:

```text
beam_1_axial.svg
beam_1_shear.svg
beam_1_bending.svg
beam_1_torsion.svg
```

For 3D beams, the shear diagram contains both `Vy` and `Vz`, and the bending diagram contains both `My` and `Mz`. SVG output has no external plotting dependency and can be opened directly in a browser.

### Analysis
- linear static solution and reaction recovery
- generalized eigenvalue modal analysis `K phi = omega^2 M phi`
- Rayleigh damping `C = alpha M + beta K`
- Newmark-beta transient integration for `M a + C v + K u = F(t)`
- arbitrary nodal time-history loads
- arbitrary time-dependent element loads evaluated at every Newmark step
- time-aware element-force recovery, so recorder output subtracts the current equivalent element load
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

Run the force-diagram example with:

```bash
./build/beam_force_diagrams_example
```

It writes a CSV file and four SVG diagrams into `beam_force_output/`.
