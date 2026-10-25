# Existing work that must not be recomputed as a new pilot

All paths below are relative to the curvature clone's cube-dim-6-bounds,
unless an absolute path is given. The inventory hashes are in
reuse_inventory.json. No source was executed in A027.

|Existing data|Actual computation and evidence status|Reuse / limitation|
|---|---|---|
|upper/numerics/body32_conditional_pilot_results.json|32face+32cube points; conditional expectation over26 remaining Cauchy variables, two512 scrambled Sobol replicates,80-node truncated logarithmic time rule|This was NOT a full fan kernel at each point. Formal costs171–1434 had material replica spread and caused a stop. They do not exclude the model.|
|upper/round003/direct_kernel/direct_kernel_pilot_results.json|Early denominator257 generic-only code encountered zero/repeated nodes and reported STRUCTURAL_FAILURE|Invalid values, including negatives, must not be recycled as G evidence.|
|upper/round003/direct_kernel/direct_kernel_generic_results.json and direct_kernel_diagnostic_summary.json|32 exact-rational inputs with denominator1000000007, all967680 terms at each point; binary64 arithmetic despite long-double spelling|True complete fan evaluations, but only numerical diagnostics. Mean(G^q)~7.73539e-6 is not an I bound. Index0 has separate exact fmpq+Arb reference.|
|upper/round003/direct_kernel/local_minimum/local_direct_results.json|121 complete967680-term generic evaluations; index24 independently exact|Index24 c=(1000000007,998836847,892738545,807464540,467291939,96968636)/1000000007 has G~2.25013386e-5. It was an endpoint, not a certified local minimum.|
|upper/round005/integral_probe/full512_results.json|Four128-point scrambled Sobol groups, every H_.001 point uses the full21 cones x46080 group; binary64 sum|Formal norm5.76520181480976e-5; no integration/error certificate. Do not redo this as a new32/128/512-point pilot.|
|upper/round005/integral_probe/arb_max_reference_results.json|Index90 is a full967680-term Arb384 recomputation; positive value~2.65642527970578e-4|One implementation calibration only. It neither certifies the other511 floating values nor the integral.|
|upper/round006/integral/gauss3_relative_results.json|The full243-node tensor Gauss3 functional,21 permutation orbits, each orbit point full-group Arb384 with exact Q(sqrt15) nodes|Q~1.02503438853374e-5 and Q^(5/6)~6.95475811139387e-5 are certified FUNCTIONAL values. No signed or absolute quadrature remainder was closed.|
|upper/round006/g and upper/round007/g_time_moments|Eight local boxes, then one time-band/Taylor box and120 permuted images|G>1/60000 on the listed regions only. A007 total covered face volume~1.25869e-6; not globalg.|
|upper/round004/primary/spectral_moment.json|Accurate spectral second moment gives G>=~6.35027e-8 throughout cube|A genuine but noncompetitive global lower bound.|
|upper/round006/integral/global_budget_results.json|True whole-face Minkowski envelope B_full~0.00139999122 and B<7/5000|A genuine but noncompetitive norm upper bound. It is not a norm lower bound.|

The H512 input denominator is1073741824; its stored seeds are6005001–6005004.
The preliminary8 indices were0,1,128,129,256,257,384,385, all included in
those512 and not additional new inputs. Index90 is replicate0/index90,
with numerators(1073741824,1009044139,794370107,872619305,961019777,807530479).
Exact point files and existing raw results are retained at their original
paths, not copied into this report.

Additional completed exact-G points in the upper-astra clone:
`/Users/zhengzhihao/Desktop/T/Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round010/short_minimum_diagnostic/`:
point1.json/point1_results.json and point2.json/batch2_results.json contain
two complete fmpq G evaluations. Point2's denominator is
50000000350000000000000; it cannot be rounded into an int64 interface.
The B023/F023 x0 value is already exact and must not be recomputed.

B015 z=(100000000000000,99999998999999,99999997999988,90268521812170,
53714088833376,157107163435)/10^14 has a complete K_.048/p/gradient result
in the upper-astra round015/spatial_minimum_search/checkpoint.json and
native_outputs/02.stdout.json. That batch explicitly performed no G or H
passes. Its gradients belong to K_.048, not G. It is a potential new
G-gradient seed, subject to a duplicate check against later completed work.

C013 is already available in the main clone under upper/round013/generic_green_engine/.
It evaluated G and all6 derivatives at the old local index24/A007 center,
and independently certified its x2 derivative by a full fmpq DD6 pass.
Its native_result.json stores exact dyadic intervals. This is an existing
G-gradient baseline, unlike the later K/H gradients; do not repeat it.
Its node-nonzero and int64 input restrictions remain mandatory.

B008 upper-astra/upper/round008/off_strata_integration/report.md also
contains completed singular-tube and off-stratum integral budget work;
it supplies no competitive I and does not exclude the true model.
