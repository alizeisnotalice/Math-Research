# R009 exact ordered-cone arrangement rays

The finite enumeration produced **46 distinct primitive rays** in the ordered cone. The catalog is a list of arrangement rays, not a triangulated fan, a Green-kernel evaluation, or a bound on g, I or C7.

For w_i=ξ_i−ξ_(i+1)≥0 and ξ₈=0, the input consists of seven coordinate-wall normals followed by the 29 mixed-sign prefix-walk normals in the frozen theory order. The remaining 35 of the 64 original body forms have nonnegative prefix normals; their boundary zeros are already represented by coordinate walls, and all 64 signs are retained. The 29-row first-hit template was regenerated and checked literally against the frozen input. The section uses Σξ_i=Σi w_i=1, with exact fractions.

Every one of the 1,947,792 lexicographic six-row subsets was processed exactly once in 1903 terminal chunks. Counts are:

| category | six-row subsets |
|---|---:|
| rank below 6 | 552,979 |
| rank 6, mixed-sign kernel | 1,066,827 |
| rank 6, nonnegative kernel after orientation | 327,986 |
| total | 1,947,792 |

Those 327,986 hits collapse to 46 distinct nonnegative gcd-primitive w vectors; zeros are retained. The numbers of rays with 1, 2, 3, 4 and 5 positive w coordinates are respectively 7, 20, 14, 4 and 1. All seven coordinate axes occur. The largest component among these primitive w vectors is 5, and the largest Σξ denominator is 24. These counts describe this arrangement and should not be conflated with the earlier spectral-vertex catalog.

The first and last 4096 lexicographic subsets were included as pilots: their FLINT integer-rank/nullspace results were checked one by one against seven independently evaluated signed 6×6 cofactors, then reused. All seven wall-row strata agree with the exact binomial counts. Each retained record contains w, ξ, normalized rational ξ, 36 incidence bits, signs of all 64 body dot products, six active witness rows and an explicit nonzero integer 6×6 minor. The separate theory proves every chamber-section vertex has such a witness and that every accepted ray is a chamber-closure vertex; enumeration alone does not prove a chamber triangulation.

The run took 2.820 s wall and 14.096 s summed chunk CPU, about 5.00 core-equivalents over wall; the largest individual process peak RSS was 35,274,752 bytes. The five-second sampler captured only the starting process because the run completed before its next sample; no sustained worker occupancy or aggregate memory peak is claimed. The 15 worker processes used single-thread settings after the serial pilot.

check.py uses independent standard-library integer arithmetic, including a Bareiss determinant, to verify every stored ray, rank witness, full incidence and body-sign list, exact section normalization, terminal spans, pilot counts and wall strata. It also verifies the pre-run and final hashes. It does not re-enumerate the 1,947,792 subsets; completeness relies on the frozen enumeration generator and the independently stated active-rank proof. No MATLAB, probabilistic sample, commit or push was used.
