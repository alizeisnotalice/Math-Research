# A033: exact low-dimensional projection peaks

All nine registered cases completed. The boundary-trace proposal is valid
with the full marginal continuity argument in `theorem.md`. The result is
an exact projection-density peak, not a new bound for C6. No spatial G/H
kernel or singular-tube budget was evaluated in this batch.

For normal rows L, put d=6-r and D=det(LL^T). The table gives the exact
rational R such that

    q_d(1,0)=sqrt(D) R/(2pi)^d,
    kappa_L=R/(2pi)^d,
    sup_y q_d(t,y)=t^-d q_d(1,0),  t>0.

The decimal columns below only display the Arb512 results saved in
`results.json`; the exact rational and formula are the certificate.

| Case | d | D | Exact R | q_d(1,0), display | B008 old peak / new peak |
|---|---:|---:|---|---:|---:|
| r1 | 5 | 6 | 67320683/3251298631680 | 5.179262651578715e-9 | 5.78583024 |
| r2, k=1 | 4 | 20 | 14455/254803968 | 1.627823948529488e-7 | 4.00679092 |
| r2, k=2 | 4 | 32 | 121723/2645913600 | 1.669754046626909e-7 | 3.90617422 |
| r2, k=3 | 4 | 36 | 2051/47038464 | 1.678588739979709e-7 | 3.88561537 |
| r3, 4110 | 3 | 64 | 3/16384 | 5.905432330297758e-6 | 2.72069905 |
| r3, 3210 | 3 | 96 | 17/110592 | 6.071852608296218e-6 | 2.64612881 |
| r3, 2220 | 3 | 128 | 23/168960 | 6.208823510083999e-6 | 2.58775339 |
| r3, 3111 | 3 | 160 | 612943/5291827200 | 5.906562742831061e-6 | 2.72017835 |
| r3, 2211 | 3 | 192 | 707/6614784 | 5.970555721256136e-6 | 2.69102322 |

All nine strict comparisons against the stored B008 peak balls passed.
Those older bounds use c=32/sqrt(3), with
Q5=2/(pi^3 c^5), Q4=3/(4 pi^2 c^4), Q3=1/(pi^2 c^3).
Only the nine peak comparisons were made: the old tube sums and global
face-integral budgets were not recomputed.

## Why the boundary calculation is legitimate

The map L pushes the full six-dimensional measure exp(-psi(xi)) dxi to
a continuous density F. Coarea gives integral_W exp(-psi)=sqrt(D)F(0).
Each fan cone has row-wise constant signs for the selected body forms.
After flipping normal rows, its nonzero normal columns lie in a positive
orthant; their nonnegative fibers are bounded. Along a generic ray
z=epsilon v, a cone with m>r nonzero normal columns therefore contributes
O(epsilon^(m-r)), hence zero in the limit. A cone with m=r contributes

    |det Rcone| / (|det(L Rcone_J)| product_(j notin J) q_j)

exactly when (L Rcone_J)^-1 v is strictly positive. Cone densities are
continuous near the chosen ray, so their almost-everywhere pushforward
identity holds there; continuity of F permits taking epsilon down to0.
No factor1/2 is assigned to fan faces, and all64 signs are enumerated.

The directions (1), (1009,1), (1009^2,1009,1) avoid every normal facet.
This follows from the integer bounds12 and288 and was also checked for
every actual facet in the pass. Primitive facet normals are retained in
`results.json`.

## Exact finite coverage and audit

The original exact fan cover is reused from round002. Every integer ray,
rate and determinant in the round003 fan was matched to that geometric
certificate. The integer fan SHA256 is
`7a57caad9b7713d9ce0896e51ab1dd402566feb1cba50829104c593a6d12a812`.
Every body form has constant weak sign on every base cone.

All4960 unordered triples were independently classified into the five
projective count types. Their counts, in table order, are480,1920,480,
640,1440. Signed permutations of normal rows induce the complete S4 on
projective classes, so these cases exhaust the triple geometry. Normal
rows for each representative are explicitly stored in `results.json`.

The sole full pass used701568 unique pattern/cone items and8709120 logical
items (967680 per case). It enumerated unique coordinate patterns, all64
signs and21 cones, restoring repeated coordinate permutations by the
class-factorial multiplicity. Its elapsed time was11.459174 seconds and
maximum recorded RSS33128448 bytes. No failure or partial case occurred.

The260 surviving terms are preserved in `selected_terms.json.gz`, including
case, permutation index, sign mask, cone, selected columns, determinant,
Cramer numerators and exact weight. The weight ALREADY includes repeated
coordinate permutations; it does not include a hidden global-sign factor.
The selected-term counts in table order are56,12,40,64,4,8,12,24,40.

`check.py` independently reconstructs those260 matrices with SymPy and
their rational sums with fmpq. It passed in0.054611 seconds. This minimal
checker is deliberately not a second exhaustive enumeration: omitted-term
completeness still depends on `compute.py`, its count ledger and the trace
proof. This limitation is explicit in `checked_results.json`.

## Files, scope and next acceptance object

- `registration.md` and `input_manifest.json`: scope, time budget, source hashes.
- `theorem.md`: coarea, trace, regularity, genericity and orbit proofs.
- `compute.py` and `preflight.json`: sole full pass and exact structural audit.
- `results.json`: nine rational peaks, Arb balls, facet data and count ledger.
- `selected_terms.json.gz`: complete surviving-term certificate.
- `check.py` and `checked_results.json`: separate selected-term arithmetic audit.
- `final_receipt.json`: frozen output hashes and final resource accounting.

There were nine scalar comparisons, zero native G/H calls, zero spatial
grid points, zero optimization calls, zero MATLAB calls and no edits to old
sources. The full enumeration must not be silently rerun: `compute.py`
refuses to overwrite an existing `results.json`.

The next bounded acceptance object would be substitution of these exact
kappa values into the existing B008 coarea estimates, with orbit weights
and existing analytic tube-integration errors preserved. Whether that
revised budget meets a useful global target has not been computed here.
No new g, I, or C6 certificate follows from this peak improvement alone.
