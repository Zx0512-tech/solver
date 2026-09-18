# Architecture

The solver separates model data, element formulation, DOF numbering, assembly and solution. The design is conceptually influenced by OpenSees (Domain/Node/Element), deal.II (separate DOF management), and the local-element/global-assembly split used by MFEM and CalculiX. No upstream implementation is copied.

## Current flow

```text
Node + Material + Section
          |
        Beam3D
   local 12x12 K
          |
 coordinate transform
          |
  global element K
          |
   Model::assemble
          |
     DofManager
          |
      global K,F
          |
 LinearStaticSolver
          |
 displacement + reaction
```

`Beam3D` is a two-node 3D Euler-Bernoulli frame element with six DOFs per node: `[ux, uy, uz, rx, ry, rz]`. It includes axial extension, Saint-Venant torsion and bending about both local section axes. The local x axis is node i -> node j; a local-y hint controls section orientation. Global stiffness is `K_global = T^T K_local T`.

## Next extensions

1. Consistent/lumped mass matrices and global mass assembly.
2. Timoshenko shear deformation as a separate formulation.
3. Distributed element loads and equivalent nodal loads.
4. Sparse matrices and sparse direct/iterative solvers.
5. Modal analysis and Newmark-beta transient integration.
6. Geometric stiffness/corotational transformation.
7. Local end-force, section-force and stress recovery.

## Reference projects

- https://github.com/OpenSees/OpenSees
- https://github.com/dealii/dealii
- https://github.com/mfem/mfem
- http://www.calculix.de/
