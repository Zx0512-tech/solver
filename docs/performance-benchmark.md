# Performance Benchmark

The benchmark executable measures the sparse v1.0 analysis path using an axial Beam3D chain. It is intended for repeatable trend tracking, not as a hard real-time performance guarantee.

Build and run:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel

./build/solver_benchmark \
  --sizes 100,500,2000 \
  --steps 50 \
  --repeats 5 \
  --csv benchmark.csv
```

## Reported columns

The CSV contains:

```text
elements
global_dofs
free_dofs
k_nnz
m_nnz
k_fill_ratio
sparse_k_bytes
dense_k_bytes
dense_to_sparse_k_memory_ratio
assembly_ms
static_total_ms
newmark_total_ms
newmark_steps
static_tip_relative_error
```

### Timing definitions

- `assembly_ms`: median time for `Model::assemble()`.
- `static_total_ms`: median end-to-end `LinearStaticSolver::solve()` time, including assembly, sparse reduction, factorization, solve and reaction recovery.
- `newmark_total_ms`: median end-to-end `NewmarkBetaSolver::solve()` time, including assembly, sparse reduction, mass factorization, one-time effective-stiffness factorization, all requested time steps and result-history storage.

The median of repeated runs is used to reduce scheduler/noise sensitivity.

### Sparse-memory estimate

`sparse_k_bytes` is an approximate compressed-sparse storage estimate based on values, storage indices and outer pointers.

`dense_k_bytes` is the bytes required by an equivalent full `double` matrix:

```text
global_dofs * global_dofs * sizeof(double)
```

The ratio is a storage comparison for K only. It is not a process resident-memory measurement.

## CI benchmark

CI runs:

```text
100 elements
500 elements
2000 elements
50 Newmark steps
5 timing repeats
```

and uploads `benchmark.csv` as the `solver-benchmark` artifact.

Absolute timings from hosted CI runners should be treated as samples, because CPU scheduling and runner hardware can vary. Trend comparisons should use comparable runners/build types and several runs.

## Correctness guard

The benchmark also checks the static tip displacement of the axial chain against:

```text
u_tip = P L_total / EA
```

and fails if the relative error exceeds the benchmark correctness tolerance.

The modal solver is intentionally not included in the scale benchmark because the current v1.0 modal reduced eigensolve is dense. Modal correctness is covered by the v1.0 verification suite.
