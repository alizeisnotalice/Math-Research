# Issue list and decision

Frozen source: `79f8acb20c418e60c5196a3c830103253e34c900`.

**Analytic theorem: ACCEPT. Core unresolved issues: none.** The independent
full proof is `proof_review.md`; historical audit status is not a premise.

| Required interface | Review result | Exact location in reconstructed proof |
|---|---|---|
| Finite maximal eligible radius selection | Valid; weights freeze and eligibility only decreases | Section 2 |
| w=small increment volume / star volume | Valid; increments, not stars, are disjoint | Sections 2 and 4 |
| Rejected cube covered using b+1 | Valid using the first rejection prefix with every radius >= rejected radius | Section 2 |
| All other positive terms earlier than k | Valid only with the last selected index among minimum-radius contributors | Section 3 |
| Captured terms <=tau | Actual selection inequality, nonnegative sub-sum | Section 3 |
| Self term | a^(-n) retained until the scalar limit | Sections 3 and 8 |
| Unpaid radii | 1<=r<t from the triangle inequality and b=a(1+1/t) | Section 3 |
| Cap containment and weighted Jacobian | Each actual unpaid I_j maps inside D(u); factor P and 2^(-n) explicit | Section 4 |
| Cap coordinate monotonicity | Valid face decomposition; no moving-volume derivative substitution | Section 5 |
| Corner shell coefficient | n/2 and finite ell_a exactly accounted for | Section 6 |
| General L1, full R^n, all positive common radii | Valid via openness, finite compact subcover, inner regularity | Section 7 |
| Finite a>2t and scalar a->infinity | Valid; no convergence of families or functions assumed | Sections 1, 6, 8 |
| B_n(t)>0 and tau choice | Valid for every integer n>=1 and t>1 | Section 8 |
| R028 rational-plus-log formula | Matches independent binomial integration and bounded Fraction check | Section 9 |
| Decimal upper bound | CONDITIONAL in this operator review on independent log/sqrt enclosure verification | Section 9 |

Paper presentation requests (expository, not mathematical blockers): spell out
the first rejection prefix, the last-index tie rule, and retention of finite-a
self and b+1 terms. Do not replace them by a reference to a PASS receipt.

Do not state any global optimizer, sharp constant, public novelty/priority,
arbitrary-rectangle theorem, or new lower bound. No such assertion was reviewed.

Operational observations: SymPy was unavailable in the attempted python3 import;
existing Python Fraction sufficed for the bounded exact check. Checkout emitted
an unrelated missing-object warning outside sparse scope. Every actually read
source was separately verified byte-identical to its frozen commit, with SHA256
and full read range recorded in `source_manifest.json`. These observations do
not create a gap in the analytic proof.
