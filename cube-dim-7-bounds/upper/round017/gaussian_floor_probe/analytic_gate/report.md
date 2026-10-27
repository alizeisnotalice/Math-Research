# R017 Gaussian-proxy analytic materiality gate

The named trace and vertex surrogate expectations have a deterministic ceiling approximately **1.1688271e-6**, strictly below the **2.5931207e-6** floor required by the rounded R016 budget for its comparison expression to reach 128. This ceiling applies to the two surrogates; it is not an upper bound on the true cube floor or on G.

For B=A⁻¹, d=det(B)^(1/2), KD=C D and C=1/(2π^(7/2)), determinant AM–GM gives

    L_trace = KD d/(tr B)^3 ≤ (KD/343) det(B)^(1/14).

The established radial identity M=D E det(B)^(1/14) therefore gives E L_vertex≤E L_trace≤C M/343. The trace expression becomes a uniform lower surrogate only after using the shape law's sign/permutation symmetry and Jensen; it is not in general an unaveraged single-shape floor. R016's specific rounded budget J₀=4018421141119/10^18 comes from F M_RB^(7/6), so M₀=(J₀/F)^(6/7) is a legal upper bound on M. This provenance is essential: an arbitrary upper bound on I would not imply such a moment bound.

A stronger exact obstruction concerns the specified proof chain. If its numerator is F M_U^(7/6), M_U≥M, and the certified floor satisfies g_cert≤E L_trace, then its reported cost expression is at least

    K = 7^(14/3) π³ / 840 ≈ 324.3106312 > 128.

The frozen proof also establishes K>128 without numerical assumptions using π>3 and 7^(2/3)>3. This is a barrier for that selected surrogate-plus-Minkowski certificate, **not a lower bound on C7**, an upper bound on actual g, or a rejection of more refined x-dependent/symmetrized methods or another numerator inequality.

The four registered 512-bit enclosures are K≈324.3106311571, U_proxy≈1.16882708483e-6, g_required≈2.59312073254e-6 and the per-shape envelope U_shape≈6.60533646443e-6. Both strict comparisons and the ratio identity U_proxy/g_required=(128/K)^(6/7) were checked by outward Arb intervals. Exact dyadic endpoints are in result.json; the R016 rational source and proof files are hash-locked.

The separate 8192-shape cache diagnostic had already finished before this analytic gate was requested; it was preserved and was not rerun. This scalar run added no shape, vertex or phase evaluation, and its strict conclusion does not use the empirical means. Wall time was 0.00425 s, CPU 0.00198 s, with peak process RSS 30,556,160 bytes. The checker repeats the same four scalar inputs and source comparisons as verification. No MATLAB, commit or push was used.

The initial checker invocation encountered a missing parenthesis. The corrected checker passes; its preregistered original bytes and a repair note are retained. The scientific scalar run was not repeated.
