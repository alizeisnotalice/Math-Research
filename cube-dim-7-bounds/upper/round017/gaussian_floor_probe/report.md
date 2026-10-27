# R017 Gaussian cube-uniform surrogate diagnostic

The 8192-shape cache gives empirical means of approximately **1.9350e-8** for the per-shape vertex surrogate and **1.4390e-7** for the symmetry-averaged trace surrogate. These are about 0.75% and 5.55% of the R016 conditional target g_required≈2.59312e-6. The observed values suggest this particular surrogate route is weak at the target scale; they do not certify either expectation or rule out a contribution missed by the finite sample.

For each frozen Dirichlet shape P, write B=A(P)^(-1), d=det(A)^(-1/2), and KD=Γ(65/2)/(π^(7/2)Γ32). Convexity of xᵀBx puts its maximum on the cube at a vertex, so

    L_vertex(P)=KD d/[max_(σ∈{−1,+1}⁷) σᵀBσ]^3

is a cube-uniform lower contribution of that shape. Averaging σᵀBσ over uniform vertex signs gives tr B, hence L_trace(P)=KD d/(tr B)^3≥L_vertex(P). After using the original shape law's full sign/permutation invariance and Jensen on Q⁻³, the expectation of L_trace is also a true cube-uniform lower surrogate. This mathematical statement is distinct from estimating its expectation by the cache mean.

| diagnostic | sample mean | nominal shape SE | eight-block SE | mean / R016 target |
|---|---:|---:|---:|---:|
| L_vertex | 1.93503e-8 | 4.00e-10 | 4.63e-10 | 0.00746 |
| L_trace | 1.43903e-7 | 2.13e-9 | 1.83e-9 | 0.05549 |

The eight original block means are in results.json. The reported standard errors describe the finite sample under the relevant sampling assumptions; they are not deterministic error bounds. The theoretical per-shape envelope, from ψ(u)≥32√2|u| and determinant AM–GM, is KD/(√2048·343)≈6.60534e-6. It bounds the surrogate but does not turn the observed mean into a strict expectation certificate.

The sole input was the original hash-verified cache, read without modification or reconstruction. The run evaluated 8192×128=1,048,576 vertex quadratics and 8192 traces. For stability, the eight cached shapes with largest estimated condition numbers were cross-checked on the same inputs using their cached scaled Cholesky factors: an additional 1024 vertex solves, eight traces and eight determinant weights. The largest relative vertex discrepancy was about 1.04e-8, below the registered 1e-6 stop threshold; no shape was dropped or replaced. All cached B matrices passed finite/SPD checks, and the sample hierarchy and theoretical envelope passed. These floating checks are not interval certificates.

Wall time was 0.100 s, CPU time 0.079 s and peak process RSS 85,508,096 bytes; the single process used one BLAS thread. The checker verifies cache and source hashes, all 8192 retained scalar observables, eight block summaries and the eight cross-check records without recomputing the vertex quadratics. No new shape sample, optimizer, Green grid, MATLAB, commit or push was used. No global g, strict expectation bound or Cn result is claimed.
