# spmv-omp-phi

High-performance **Sparse Matrix–Vector Multiplication (SpMV)** in C with **OpenMP**, designed for highly parallel architectures such as Intel Xeon Phi and multi-core CPUs.

## Overview

This library implements \( y = A \cdot x \) for sparse matrices \( A \) stored in several formats. Each format offers different trade-offs between memory layout, cache behavior, and parallelism, so the best choice depends on the matrix structure and target hardware.

## Supported storage formats

| Format | Abbreviation | Description |
|--------|--------------|-------------|
| **CRS/CSR** | Compressed Row Storage | Row-oriented; one row per thread, good for row-parallel SpMV. |
| **CCS/CSC** | Compressed Column Storage | Column-oriented; one column per thread, scatter into result. |
| **BCRS/BCSR** | Block Compressed Row Storage | Row-blocked; groups consecutive row nonzeros for better locality. |
| **ELL/ELLPACK** | ELLPACK | Fixed-length rows (padded); regular memory access, good for vector/SIMD. |

## Building

Using the provided Makefile (recommended):

```bash
make          # build libspmv.a
make clean    # remove object files and library
```

Or compile manually with an OpenMP-enabled compiler (e.g. GCC, Clang, ICC):

```bash
gcc -fopenmp -O3 -I include -c src/mmio.c src/util.c src/dense.c src/sparse.c src/blocked.c src/reorder.c
ar rcs libspmv.a src/*.o
```

Include the `include/` directory when compiling your code (`-I include`). Link with `-fopenmp` and the math library if needed (`-lm`).

## Usage outline

1. **Read a matrix** (Matrix Market format) into a dense `MATRIX`:
   ```c
   FILE *f = fopen("matrix.mtx", "r");
   MATRIX *m = ReadMatrix(f);
   ```

2. **Convert** to a sparse format, e.g. CRS:
   ```c
   CRS *crs = CreateCRS(m);
   ```

3. **Allocate** input vector `x` and output vector `y` (length = number of rows for \( y \), length = number of columns for \( x \)).

4. **Compute** \( y = A \cdot x \):
   ```c
   MultiplyCRS(crs, x, y);
   ```

5. **Clean up**:
   ```c
   DestroyCRS(crs);
   DestroyMatrix(m);
   ```

Matrix Market (`.mtx`) files are supported via the bundled MM I/O code (see [Matrix Market](http://math.nist.gov/MatrixMarket/)).

## Project layout

```
.
├── include/
│   ├── sparse.h   # CRS, CCS types and SpMV
│   ├── dense.h    # Dense MATRIX and I/O
│   ├── blocked.h  # BCRS, ELL types and SpMV
│   ├── reorder.h  # Column reordering (e.g. Gray code)
│   ├── util.h     # Timing, quicksort, helpers
│   └── mmio.h     # Matrix Market I/O
├── src/
│   ├── sparse.c   # CRS/CCS create, multiply, destroy
│   ├── dense.c    # Read/write matrix, dense SpMV
│   ├── blocked.c  # BCRS/ELL create, multiply, destroy
│   ├── reorder.c  # Reordering for cache/conflict reduction
│   ├── util.c     # dtime(), QuickSort, findMax, etc.
│   └── mmio.c     # Matrix Market implementation
├── Makefile
├── README.md
└── LICENSE        # Public domain (Unlicense)
```

## Fixes and robustness (review)

The following issues were identified and fixed during review:

- **MultiplyCCS**: Removed incorrect use of `temp` and overwriting of `r[rowInd[j]]` after the inner loop (which also used `j` out of scope). CCS SpMV now correctly scatters per column.
- **MultiplyMatrix**: Fixed parallel reduction — each row now uses a private `result` and correct per-row accumulation.
- **DestroyELL**: Corrected to free `c->values[i]` and `c->indices[i]` (and use `c->nrows`) instead of wrongly using `m->mel`.
- **MultiplyELL**: Skip padding entries where `indices[i][j] == -1` to avoid out-of-bounds and wrong values; result vector size set to `nrows`.
- **util.c**: Added `#include <sys/time.h>` for `gettimeofday` and `struct timeval`.

## Implementation notes

- **CRS**: Parallel over rows; each thread accumulates a row dot product. No scatter, good for cache and OpenMP.
- **CCS**: Parallel over columns; each thread scatters contributions into `y`. Can cause write conflicts if the same row appears in multiple columns (current code uses parallel for over columns; for correctness with conflicting writes, use atomic adds or a different strategy).
- **BCRS**: Blocks are contiguous in memory; inner loop over block entries uses `#pragma ivdep` to hint vectorization.
- **ELL**: Rows padded to `max_entries_per_row`; padding indices are -1 and skipped in the multiply. Regular stride per row can help SIMD.

## Possible improvements and extensions

- **Build system**: Add a `Makefile` or CMake to build library and a small driver (e.g. read `.mtx`, run one format, report time).
- **COO (Coordinate)**: Store `(i, j, val)` and optionally sort by row for a simple SpMV; good for very irregular matrices and GPU-style kernels.
- **HYB (COO + ELL)**: Use ELL for “regular” rows and COO for the tail of long rows to limit padding.
- **Symmetric matrices**: For symmetric storage (e.g. only lower triangle), add a code path that uses both `A[i,j]` and `A[j,i]` when needed.
- **Vectorization**: Experiment with explicit SIMD (e.g. AVX2/AVX-512) for ELL and BCRS inner loops; keep OpenMP for thread parallelism.
- **Format auto-selection**: Heuristic (e.g. row-length variance, nnz) to choose CRS vs ELL vs BCRS.
- **Reordering**: `reorder.c` provides column intersection and Gray-code ordering; document how to apply the permutation to the matrix and how it improves locality or reduces conflict.
- **Testing**: Add unit tests (e.g. compare CRS/CCS/ELL/BCRS/dense result against a reference) and optionally a small benchmark harness.
- **Error handling**: Check `malloc` and file I/O returns; avoid `exit(1)` inside library code when possible (return error codes or let the caller decide).
- **Caller-allocated output**: Allow the caller to pass a pre-allocated `y` buffer to avoid internal allocation and to support repeated runs without reallocating.

## License

Public domain. See [LICENSE](LICENSE) (Unlicense).
