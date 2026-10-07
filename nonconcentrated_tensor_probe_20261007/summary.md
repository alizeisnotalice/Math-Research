# Implicit overlapping tensor mixture: three frozen numerical rounds

B41 shared-resolution mechanism on the full finite measure. No parameter or lambda selection; no actualgeom or FIRST gates. Values below are independent A/B cross-product finite-node energies. `cert` means the screened, count-bound-certified residual subset; `upper` retains every unresolved residual query and floating boundary ambiguity. These are numerical screens, not interval certificates.

| round | n | samples / side | nodes | Gamma1 cert | Theta hard cert | Gamma vstar cert | Gamma vstar raw cert | Gamma vstar upper |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 16 | 32 | 65 | 0.0402756 | 0.0296252 | 0.15996 | 0.15996 | 0.797644 |
| 1 | 32 | 32 | 65 | 0.0772514 | 0.0595685 | 0.423656 | 0.423656 | 0.670259 |
| 1 | 64 | 32 | 65 | 0.032322 | 0.0258653 | 0.175064 | 0.175064 | 0.254684 |
| 2 | 16 | 64 | 129 | 0.0661513 | 0.0494122 | 0.233893 | 0.233893 | 0.80432 |
| 2 | 32 | 64 | 129 | 0.0496042 | 0.0386504 | 0.277392 | 0.277392 | 0.467956 |
| 2 | 64 | 64 | 129 | 0.0345008 | 0.0280064 | 0.208507 | 0.208507 | 0.279888 |
| 3 | 16 | 128 | 257 | 0.0637015 | 0.0466642 | 0.223607 | 0.223607 | 0.749701 |
| 3 | 32 | 128 | 257 | 0.0606369 | 0.0478353 | 0.311164 | 0.311164 | 0.481081 |
| 3 | 64 | 128 | 257 | 0.0350376 | 0.0279542 | 0.202223 | 0.202223 | 0.280009 |

The final n=16/32/64 values are three separate fixed-dimensional inputs; no growth law is fitted. All rounds use the same source formula and frozen log(lambda)=-n log64, but different seeds, sample sizes and nested time grids. The upper column omits the microbox gate: it is a full hard-relaxation upper envelope that retains unknown queries, and does not prove a high residual energy. The certified residual is the exact complement of the box-existence event Eexists, conditional on computed winners/counts; no finite greedy-cover complement is used.

Approximate t3 intervals and bounded Hoeffding intervals for each finite grid statistic are stored in `round*.json`; four replicate pairs give little precision, and the bounded intervals are wide. Nested coarse/fine differences are discontinuous-integrand resolution diagnostics, without a deterministic continuum quadrature error bound. Atomic-to-L1 transfer and arbitrary allowed-scale or actualgeom qualification remain open.

All source/cone draws use the complete probability measure. Distinct atom counts are 256^n, with log10 counts approximately 38.53,77.06,154.13. Overlapping source sites carry summed component masses. No sparse finite subsample replaces the measure when computing maximal averages.
