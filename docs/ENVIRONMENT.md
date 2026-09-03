# Environment

Full reproduction is a claim about a computation, and a computation is not separable from the arithmetic that runs it. The certificates record interval endpoints produced with directed rounding; a different arithmetic backend can move the last digits, and a different summation order can move a decision. The environment is therefore pinned, not recommended.

## Reference versions

| component | version |
| --- | --- |
| Python | 3.14.4 (CPython) |
| numpy | 2.5.1 |
| mpmath | 1.3.0 |
| sympy | 1.14.0 |
| scipy | 1.18.0 |
| matplotlib | 3.11.1 |
| BLAS | scipy-openblas 0.3.33.112.0 |
| LAPACK | scipy-openblas |

## Determinism

These variables are set by the container and should be set by hand for a local run. Each is here for a reason:

- `OMP_NUM_THREADS=1` — the summation order of a multithreaded BLAS is not guaranteed reproducible.
- `OPENBLAS_NUM_THREADS=1` — the same, for the OpenBLAS backend.
- `MKL_NUM_THREADS=1` — the same, should the backend be MKL.
- `PYTHONHASHSEED=0` — set iteration order influences any path that derives an ordering from it.
- `SOURCE_DATE_EPOCH=0` — a reproducible timestamp for any artefact that would carry one.

## What the certificates themselves record

| certificate | python | platform |
| --- | --- | --- |
| bridge_atlas_panel | 3.14.4 | Linux-7.0.0-30-generic-x86_64-with-glibc2.43 |
| bridge_metric_path | 3.14.4 | Linux-7.0.0-30-generic-x86_64-with-glibc2.43 |
| face_traversal_leaf | 3.14.4 | Linux-7.0.0-30-generic-x86_64-with-glibc2.43 |

Certificates that record no environment were produced by exact rational or interval arithmetic whose result does not depend on the backend beyond the pinned mpmath.

## Shipped inputs whose recorded platform is localised

`platform.platform()` asks the operating system for the processor name, and the operating system answers in its own language. Everything this repository rebuilds is produced under `LANG=C`, where that name comes back empty. The files below are INPUTS: no producer shipped here writes them, so they cannot be regenerated without shipping a chain this repository deliberately does not carry. Their provenance is published as recorded rather than edited, because a provenance corrected by hand is no longer one.

| shipped input | recorded platform |
| --- | --- |
| exact_gluing | Linux-7.0.0-28-generic-x86_64-inconnu-with-glibc2.43 |
| face_preliminary | Linux-7.0.0-28-generic-x86_64-inconnu-with-glibc2.43 |
| halo_metric | Linux-7.0.0-28-generic-x86_64-inconnu-with-glibc2.43 |
