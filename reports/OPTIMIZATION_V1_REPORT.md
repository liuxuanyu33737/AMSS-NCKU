# AMSS-NCKU CPU Optimization V1 Report

## Environment

- Ubuntu 22.04 WSL
- OpenMPI
- GCC 11.4

## Baseline

- MPI processes: 8
- Evolution steps: 20
- Total Evolve Time: 3852.32 seconds
- Total Running Time: 4009.23 seconds
- Real Time: 66m49.233s

## Optimization V1

- MPI processes: 8
- Evolution steps: 20
- Total Evolve Time: 821.391 seconds
- Total Running Time: 824.376 seconds
- Real Time: 57m35.947s

## Performance Improvement

- Evolution speedup: 3852.32 / 821.391 = approximately 4.69x
- Running speedup: 4009.23 / 824.376 = approximately 4.86x
- Real Time improvement: 66.82 min -> 57.60 min, approximately 16%

The core computation time was reduced substantially. However, wall-clock real
time is also affected by MPI communication, WSL scheduling, and other runtime
overheads.

## Optimization Summary

### Compiler flags

- Removed `-pg` profiling instrumentation overhead.
- Enabled native CPU instructions with `-march=native`.
- Enabled additional Fortran inlining with `-finline-functions`.

### Six-point `polint` specialization

- Added a fixed path for the production case `n == 6`.
- Expanded nearest-node selection and the Neville recurrence stages.
- Preserved the original generic algorithm for `n != 6`.
- Preserved the interface, interpolation order, mathematical algorithm, and
  physical computation flow.
