# v1.0 Verification Suite

The v1.0 verification suite is an end-to-end acceptance layer for the solver. It complements the focused unit/regression tests by checking complete analysis paths against analytical structural-mechanics results.

Run it with:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
./build/test_v1_verification
```

The suite is also registered with CTest as `v1_verification`.

## Acceptance cases

### 1. Combined 3D cantilever static response

A single Beam3D is loaded simultaneously by axial force, local-y force, local-z force and torsion.

The suite checks:

```text
ux = Fx L / EA
uy = Fy L^3 / (3 E Iz)
uz = Fz L^3 / (3 E Iy)
rx = Mx L / GJ
```

and the six support-reaction equilibrium terms relevant to those loads.

### 2. Uniform distributed load and section equilibrium

A cantilever under a full-span local-y UDL is checked against:

```text
u_tip = q L^4 / (8 E Iz)
theta_tip = q L^3 / (6 E Iz)

Vy(x) = q (L - x)
Mz(x) = q (L - x)^2 / 2
```

This verifies the chain:

```text
element load
  -> equivalent nodal load
  -> global solve
  -> section-force recovery
```

### 3. First Euler-Bernoulli cantilever bending mode

An 8-element planar cantilever is checked against the analytical first bending circular frequency:

```text
omega_1 = beta_1^2 sqrt(E Iz / (rho A L^4))
beta_1 = 1.875104068711961
```

The acceptance tolerance is intentionally tighter than a typical engineering reporting tolerance while still accounting for finite-element discretization.

### 4. Newmark-beta free vibration

A one-DOF axial cantilever with initial displacement and zero velocity is integrated for one analytical period.

The suite checks the half-period and full-period displacement against:

```text
u(T/2) = -u0
u(T)   =  u0
```

### 5. Recorder and envelope chain

The same transient result is passed through:

```text
NewmarkResult
  -> NodeRecorder
  -> EnvelopeRecorder

NewmarkResult
  -> SectionRecorder
  -> EnvelopeRecorder
```

The displacement and axial-force envelopes are checked against the analytical amplitude.

### 6. Rayleigh damping targets

`RayleighDamping::fromModalTargets` is checked at both requested target frequencies.

### 7. Rectangle-section normal stress

The chain:

```text
static solve
  -> section N/M
  -> BeamSectionStressRecovery
```

is checked against:

```text
sigma_x = N/A - Mz y/Iz + My z/Iy
```

## Relationship to the focused regression suite

The v1.0 suite does not replace the existing test groups. The focused tests continue to cover edge cases such as:

- local/global Beam3D transformations
- consistent mass
- partial and time-dependent loads
- concentrated-load section-force jumps
- CSV/SVG force output
- section geometry and unsupported stress guards
- result recorders
- sparse assembly
- sparse reduced-system solvers

The v1.0 suite is the release-level, end-to-end acceptance layer above them.
