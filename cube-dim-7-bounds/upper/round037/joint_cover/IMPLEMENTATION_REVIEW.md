# R037 source-level implementation review

This is a read-only source review, not a production rerun, new LP solve, or independent reconstruction of the 48-cube atom table. The separate frozen-instance audit owns that reconstruction. The all-dimensional four-cube proof in `SPLITTING_PROOF.md` is independent of the stochastic implementation.

The reviewed production directory is `/Users/zhengzhihao/Desktop/T/cube-r030-compute/cube-dim-7-bounds/upper/round037/compute/`.

## Mathematical implementation interfaces

The following core interfaces match the intended exact model at source level.

* `joint_atoms.py` partitions the two active coordinates and inactive L-infinity radial slabs, and multiplies exact rational cell area by `(2*rho_hi)^(n-2)-(2*rho_lo)^(n-2)`. It aggregates Q/E masks only where the original union is nonempty, which is sufficient for increments and rejected volume. Strict boundaries are null sets.
* `initial_flow` processes a nonincreasing-radius order once and permanently skips ineligible candidates. This is equivalent to largest-eligible selection: adding nonnegative fixed p terms cannot restore an ineligible candidate. Increment weights p are fixed when selected; the final combined load q includes the proportional rejected-volume charges.
* The exact proportional rejected flow uses P_sigma>tau and payments d_sigma*p_j/P_sigma. The routine checks total payment equals the exact union volume.
* The witness graph uses strict open-box intersection. Complete maximal-clique enumeration and nonnegative box Helly give the true load. Hitting enumeration caps is a censored case, not a complete-clique certificate.
* The mass LP has the correct min-mass primal and max-packing dual directions. Its floating outputs are proposals; the implemented rational primal repair adds uncovered residual to a containing maximal clique, and the dual repair scales to enforce every complete clique constraint.
* `joint_compress.py` uses full allowed-graph true-twin classes, exact payment-preserving maximal-budget anchors, the correct current-q boundary load, and exact lexicographic safe-candidate selection. It verifies the resulting true load does not increase. The absence of a separately stored L_rest does not invalidate that direct complete-clique load check, but it leaves the requested detailed bill output incomplete.

These observations validate the mathematical form of the implementation. They do not by themselves certify every emitted case or every complete-clique claim; exact frozen case data and independent verification are still required.

## Registration and output deviations

The actual batch must be identified by its frozen coordinates, histories and code hashes. It is not identical to the written SPEC case set.

1. Tie orders use per-vertex hash domains `tie-v-i`, whereas the SPEC describes radius-class domain permutations. Equal L0 orders retain the first encountered order instead of the lexicographically smallest order.
2. Template A's direction iteration is `e1,-e1,e2,-e2`, whereas the SPEC lists `e1,e2,-e1,-e2`. The full core set is unchanged, but the initial probe-slot ordering before seed shuffling differs, so a selected probe subset can differ geometrically. B's core direction order also changes IDs and subsequent tie histories.
3. C rejects a whole family candidate at the first duplicated sign vector, rather than allowing up to 256 redraws of that vertex. Its actual retry seed domain is `retry-k` rather than the written `retry/k`.
4. The stopped-pool diagnostic for unsafe moves which nevertheless decrease the true load is absent from the initially reviewed code.
5. The output combines negligible-gap and unresolved-gap statuses and does not retain all case-level primal/dual vectors and charge data. These limitations require explicit ledger interpretation or incremental archival; they must not be repaired by pretending the original run used a different registration.

The geometry/history deviations still produce lawful rational finite instances and lawful first-rejection payments. They invalidate a claim that raw labels are rigorous statistics for an identical pre-registered finite case set. Only an audited actual rational instance supports an exact counterexample claim. H failure is a scientific result, not an internal invariant failure; resource censors are likewise not scientific counterexamples.

## Structural diagnosis of the frozen failure

The frozen A/S0/case6/partial witness has positive q only at vertices 0,2,3, all with the same largest budget H. Proposition 1 of `SPLITTING_PROOF.md` therefore makes the concentration saving s zero for every clique, including candidates outside the registered pool. Enlarging the pool cannot repair this instance of the strict-saving rule.

Its leading four witness centers are `(4,0),(-4,0),(0,4),(0,-4)` with radius 4, the same C4 contact pattern as the standalone exact proof. Zero outside charge means fixed-q independent splitting of {0,1} is a legal alternative even if smaller zero-charge witnesses distinguish their full graph neighborhoods. If those outside variables may later be reactivated in an LP, the full allowed-support restriction applies.

## Reviewed identities

The following hashes were read after the source inspection. No file in the production directory was edited by this review.

    joint_geometry.py 858104b4e0c7f2789c69f5c442ef284d62d16eebe912748d944c16ed1718eb87
    joint_atoms.py 6e12c2da2837e2d565b3a1e903702540237fb2d1efbf2382ca579d149c48d061
    joint_core.py df8e193bb1970215045bd358e4561a9cc8f84a717fca67e7098fe8098115ba0c
    joint_compress.py 77a64b0deecb9f5591ce6e191053fea75cfc68be89fd40f970f805ea8b71f5ce
    r037_joint.py 11d6a69f4f8dbd29e66bbd13d4a7abfc171d1849941937e754046b82aeaa07fc
    r037_failure_witness.json f82c98c9fcfc8a9ef769aa645442e2f8ef0130b4b9e97b0d97066013498f2ac2
    r037_failure_geometry.json 6f77aca883bc6d48faa3d7f5896ec0195cc8b516d3fc1b909a9c442729a307e2

No new production computation or previous baseline verification was run for this review.
