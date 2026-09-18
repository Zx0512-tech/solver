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

### Section geometry and stress recovery

The section layer keeps the original property-only API and adds geometry-aware convenience types:

```cpp
fem::GeneralSection general(A, Iy, Iz, J);
fem::RectangleSection rectangle(size_y, size_z);
fem::CircularSection circle(radius);
```

`RectangleSection` derives `A`, `Iy`, `Iz` and an engineering Saint-Venant torsion-constant approximation from its dimensions. `CircularSection` derives the exact solid-circle properties.

Normal stress is available for every section:

```text
sigma_x = N/A - Mz*y/Iz + My*z/Iy
```

For geometry-aware sections the solver also validates that the requested `(y,z)` point lies inside the section.

Supported full point-stress recovery in v1.0:

- `RectangleSection`: axial/biaxial-bending normal stress plus classical rectangular `Vy/Vz` transverse shear.
- `CircularSection`: axial/biaxial-bending normal stress plus solid-circle Saint-Venant torsional shear.
- `GeneralSection`: axial/biaxial-bending normal stress only because boundary/shear geometry is unknown.

Unsupported combinations throw explicitly instead of dropping a stress contribution silently. In particular, rectangular torsional point stress and circular transverse-shear point stress are intentionally deferred.

Example:

```cpp
const double sigma_x =
    model.beamNormalStressAt(
        beam_id,
        x,
        y,
        z,
        result.displacement);

const auto stress =
    model.beamStressAt(
        beam_id,
        x,
        y,
        z,
        result.displacement);

const auto extrema =
    model.beamNormalStressExtrema(
        beam_id,
        x,
        result.displacement);
```

`beamNormalStressAt` always recovers the axial+bending normal stress, even when the section's full shear/torsion point-stress model is not implemented. `beamNormalStressExtrema` returns `sigma_min` / `sigma_max` and their section coordinates for rectangles and solid circles.

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

Run the section-stress example with:

```bash
./build/beam_stress_recovery_example
```

It reports the root-section point stress and minimum/maximum normal stress for a rectangular Beam3D.


### Result recorders and dynamic envelopes

The result layer provides three focused recorders:

```cpp
const auto node_history =
    fem::NodeRecorder{}.record(node_id, newmark_result);

const auto section_history =
    fem::SectionRecorder{}.record(
        model, beam_id, x, newmark_result);

const auto stress_history =
    fem::SectionRecorder{}.recordNormalStress(
        model, beam_id, x, y, z, newmark_result);
```

`EnvelopeRecorder` returns minimum, maximum and maximum-absolute response together with the occurrence time and step:

```cpp
const auto uy_envelope =
    fem::EnvelopeRecorder{}.record(
        node_history,
        fem::NodeResponseQuantity::Displacement,
        fem::Dof::UY);

const auto mz_envelope =
    fem::EnvelopeRecorder{}.record(
        section_history,
        fem::BeamSectionForceComponent::Mz);

const auto sigma_envelope =
    fem::EnvelopeRecorder{}.record(stress_history);
```

Each envelope contains:
- minimum value, time and step
- maximum value, time and step
- maximum absolute magnitude, signed value, time and step

Run:

```bash
./build/result_envelope_example
```

to see tip-displacement, root-bending-moment and root-normal-stress envelopes from a Newmark transient analysis.
