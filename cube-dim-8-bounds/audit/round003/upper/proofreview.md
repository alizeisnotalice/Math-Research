# R003 selector-coupling audit: orthant energy and a feasible sharp family

2026-10-27. **Verdict: PASS for the stated selector-model calibration, its sharp-coefficient obstruction, and the comparison of the two coarse load formulas.** The conclusion supplies neither a lower bound for \(C_8\) nor a new accepted operator upper bound. The frozen author package and this audit's actual execution are identified in `manifest.json`.

The tested model statement is that, with \(d=2^n-1\),

\[
\frac{S^2}{2d}-\frac12\sum_{j\in J}w_j^2
\le E_J\le\tau S,
\quad 0\le w_j\le M=a^{-n}.
\]

The coefficient \(1/(2d)\) cannot be increased uniformly while retaining that exact diagonal correction, or the weaker correction \(MS/2\), even for families produced by the actual finite selector. The proof is geometric and uses the real increment weights; it is not a claim about an arbitrarily weighted abstract graph. Nonuniform inequalities depending on a large-load regime, actual radii or cap occupancy remain possible.

All newly authored files are confined to `audit/round003/upper/`. Previously frozen R001 upper and Green audit files are unchanged. This role did not execute the author program, repeat its 32,640-pair loop, optimize, scan other candidates, invoke MATLAB, spawn agents, commit or push.

## 1. Energy calibration, including the directed factor and diagonal

Normalize the observation point to zero and the selected last minimum-radius positive contributor to radius 1, with index \(k\). The earlier R001 rule gives all other positive contributors indices less than \(k\), nonincreasing radii, and the unpaid set \(J\) has \(1\le r_j<t\). Here \(a>2t,t>1,b=a(1+1/t)>a\), and \(|y_j|_\infty<ar_j\). The actual weights are

\[
w_j=|Q_j\setminus\bigcup_{i<j}Q_i|/(2ar_j)^n\le a^{-n}=M.
\]

Define

\[
S=\sum_{j\in J}w_j,\quad D_2=\sum_{j\in J}w_j^2,\quad
E_J=\sum_{i,j\in J,\ i<j}w_iw_j1_{bQ_i}(y_j).
\]

At each selected \(j\), the sum over earlier \(J\)-indices is a nonnegative sub-sum of the actual eligibility test, hence at most \(\tau\). Multiplying by \(w_j\) and summing gives \(E_J\le\tau S\). The use of \(i<j\) counts each captured ordered-in-time pair once; there is no second reversed contribution.

Assign zero coordinates to the positive side when partitioning centers by sign. If \(i<j\) lie in the same sign class, every coordinate satisfies

\[
|y_{i,\ell}-y_{j,\ell}|\le
\max(|y_{i,\ell}|,|y_{j,\ell}|)<ar_i<br_i,
\]

because \(r_i\ge r_j\). Thus every earlier same-class center captures every later one, with a strict margin inside the **open** \(bQ_i\).

Reflect the full configuration so that every anchor coordinate is nonnegative. This preserves the selector, weights and all distances. An unpaid center cannot be coordinatewise nonnegative: otherwise \(r_j\ge1\), \(0\le y_{k,\ell}<a\), and \(0\le y_{j,\ell}<ar_j\) imply \(|y_{k,\ell}-y_{j,\ell}|<ar_j<br_j\) in every coordinate. It would capture the anchor, a contradiction. Therefore at least one coordinate is strictly negative, leaving at most \(d=2^n-1\) sign classes.

Writing their masses as \(S_\sigma\), actual captured pairs give

\[
E_J\ge\frac12\left(\sum_\sigma S_\sigma^2-D_2\right)
\ge\frac{S^2}{2d}-\frac{D_2}2
\ge\frac{S^2}{2d}-\frac{MS}2.
\]

The first factor \(1/2\) removes the doubled cross terms in each square; the negative diagonal is exactly \(-D_2/2\). Since \(D_2\le MS\), replacing it by \(MS\) makes the lower bound weaker, which is the required direction. The same proof without excluding the anchor's sign class gives the requested \(2^n\)-orthant calibration. With \(E_J\le\tau S\), division only for \(S>0\) gives \(S\le2d\tau+dM\); the zero case is immediate.

## 2. The stronger direct linear consequence must also be compared

For each nonempty sign class, take its last selected positive-weight \(J\)-term, denoted \(j_\sigma\). Every other member of that class is earlier and captures this term. Consequently the real eligibility test gives

\[
S_\sigma-w_{j_\sigma}\le\tau,
\quad S\le d\tau+\sum_\sigma w_{j_\sigma}\le d(\tau+M).
\]

This is stronger than the mass bound inferred through the aggregate energy by a factor of two in its \(\tau\) term. The author includes it in §3.3 and evaluates its competitiveness. Thus sharpness of the quadratic coefficient is not being misrepresented as sharpness of every possible selector mass estimate. The construction below saturates this direct linear estimate for \(m\ge2\), but it does not imply that the complete R001 load bound, including captured terms and the anchor, is sharp.

## 3. The fixed probe is a real selector run, including the last anchor

For the registered \(n=8,t=2,a=8\) probe, \(b=12\), and every supplied cube has radius 1 with center in \(\{-7,+7\}^8\). Put the all-positive center last. All radii are equal, so any prescribed tie order is allowed by maximum-radius selection. Each star contains zero strictly, since \(7<8\).

Any two distinct sign vectors differ in at least one coordinate; their sup-distance is exactly 14. Hence they are farther apart than 2, so the small open cubes are disjoint, and farther apart than 12, so no distinct center captures another. Every frozen increment is its entire cube and every weight is exactly

\[
M=8^{-8}=1/16777216.
\]

Every eligibility test is zero, including the final positive anchor test. The registered positive threshold \(\tau=2^{-30}\), although smaller than \(M\), is sufficient because the test is zero; no lower bound \(\tau\ge M\) is required. All unselected cubes remain eligible throughout, and the selector can select the entire family in the prescribed order. The anchor is a minimum-radius positive contributor and is the last such index, exactly as required by the R001 tie convention.

All other 255 cubes genuinely lie in \(J\): they contribute at zero, precede the anchor, have radius \(1<t\), and do not capture its center. Consequently

\[
S=255M,\quad D_2=255M^2=MS,\quad E_J=0.
\]

The calibrated lower bound is exactly zero. Replacing \(1/(2\cdot255)\) by \(1/(2\cdot254)\) would instead require

\[
E_J\ge\frac{S^2}{508}-\frac{D_2}2
=\frac{255}{142989288169013248}>0,
\]

which is false. More generally any strict increase of the coefficient gives \((\kappa-1/510)S^2>0\), so the one calibrated equality excludes all such uniform increases, not just the displayed rational test value. The weakened diagonal version fails as well because \(D_2=MS\).

The source reports checking \({256\choose2}=32640\) distinct center pairs, including pairs with the anchor. The energy sums only pairs within \(J\); there are \({255\choose2}=32385\) potential such pairs, and every one is uncaptured. There is no pair-count or diagonal mismatch. This audit does not enumerate them: the single sign-distance argument establishes all of them.

For compatibility with a true \(L^1\) witness, take \(f=2\,1_{(-16,16)^8}\) and \(\lambda=1\). Each star has coordinate extent at most 15 in absolute value and is strictly inside the input support. Its average is exactly 2. Supplying the associated bad small cubes is therefore a permitted input to the actual finite covering selector. This is feasibility of the model family, not a computation of a weak ratio or a \(C_8\) lower bound.

## 4. Sharpness is not confined to a vanishing-mass one-point family

For the concrete eight-dimensional sequence, keep \(t=2\), take integer \(m\ge1\), \(a=8m\), \(b=12m\). Each non-all-positive sign class has \(N=m^8\) radius-one centers with coordinate magnitudes

\[
6m+1,6m+3,\ldots,8m-1.
\]

Add the all-positive anchor \((8m-1)\mathbf1\) and select it last. Endpoint bounds establish all required geometry:

- Every star contains zero with margin at least 1, since every center magnitude is at most \(a-1\).
- Distinct centers in a class differ by at least 2 in some coordinate. The corresponding radius-one **open** cubes are disjoint even at the equality case: the adjacent endpoints belong to neither open interval. There is no positive-volume overlap to reduce the frozen increments.
- The within-class diameter is at most \(2m-2<12m\), so all earlier terms in a class capture later terms.
- Centers in different sign classes have an opposite-sign coordinate with separation at least \(12m+2>12m\). No cross-class capture occurs.
- The anchor is at distance at least \((8m-1)+(6m+1)=14m>12m\) from every unpaid center in at least one negative coordinate. Its incoming test is therefore exactly zero.

All small cubes, including the anchor, are disjoint. Thus \(M=(8m)^{-8}\) is the actual common increment weight. At any intermediate selection, a class member receives only its previously selected class mass, at most \((N-1)M\). Taking \(\tau=NM=8^{-8}>0\) always works. For \(m\ge2\), the smaller \(\tau=(N-1)M>0\) also works, including its final in-class equality, because original eligibility uses **\(\le\tau\)** rather than a strict inequality. There are no nonselected larger-radius cubes to interfere; all radii equal 1. The prescribed anchor-last order is a legal maximum-radius order at every stage.

This gives exact, rather than asymptotic, identities

\[
S=255NM=255/8^8>0,\quad D_2=255NM^2=MS,
\]
\[
E_J=255\frac{N(N-1)}2M^2
=\frac{S^2}{510}-\frac{MS}2.
\]

As \(m\to\infty\), \(M\to0\), the total unpaid mass remains fixed and positive, and \(E_J/S^2\to1/510\). For the exact threshold \((N-1)M\), one also has \(S=255(\tau+M)\), saturating the direct linear class estimate. This does not make the energy upper bound \(E_J\le\tau S\) sharp: here \(E_J=\tau S/2\), a distinction relevant to any further selector analysis.

The same construction works in general dimension and for arbitrary fixed \(t>1\) with \(q=1+1/t\), provided \(L>2t\), \(L(2-q)\ge4\), \(a=Lm\), and magnitudes \((L-2)m+2h-1\). Then cross-class separation is at least \(2(L-2)m+2\ge qLm+2\), and the anchor separation is \((2L-2)m>qLm\). The fixed mass is \(d/L^n\). The author's additional fixed-\(a,t\) singleton construction, with \(s=(2+q)/4\), is also valid: \(q/2<s<1\) makes \(2sa>b\), and \(a>2t>2,s>3/4\) ensure small-cube disjointness. That family proves the coefficient obstruction even at each fixed admissible parameter pair.

The larger family has a compatible bounded compactly supported input \(f_m=2\,1_{(-2Lm,2Lm)^n}\), since its centers have magnitudes below \(Lm\) and its stars have radius \(Lm\). Thus it is a genuine bad-average selector family for each finite \(m\). No limit of these input functions is needed for the coefficient sharpness statement.

## 5. Competitive comparison is correct and narrowly scoped

The old load decomposition has \(W\le\tau+M+S\), while its finite cover factor is \(1+A/\tau\), \(A=(1+1/t+1/a)^8>1\). The quadratic energy estimate alone gives

\[
W\le511\tau+256M.
\]

Multiplication by the cover factor yields a scalar expression at least \(511A>511\). Its infimum over admissible parameters is 511: first let finite \(a\to\infty\), then \(\tau\downarrow0\), then \(t\to\infty\). Each step is a scalar comparison in the already established covering skeleton; no new operator limiting assertion is needed.

The stronger last-index bound instead gives \(W\le256(\tau+M)\). Its analogous scalar-family infimum is 256, by the same ordered limits; every individual expression is at least \(256A>256\). These numbers are minima of the two **coarse bound families**, not lower bounds for the true operator constant or for every possible selection coupling.

A purely analytic existing-cap comparison suffices: \(B_8(2)<1+8\log2<9\) and \((1+1/2)^4=81/16\), hence

\[
F_{\rm cap}(2)<(81/16+3)^2=16641/256<256<511.
\]

The strict \(\log2<1\) follows directly from its integral over \([1,2]\). Thus no old decimal optimum, new scalar grid, or numerical search is involved. The author correctly concludes that these uniform orthant formulas cannot compete with the established cap family. The construction's fixed positive mass is still small; it does not claim to fill the cap comparison region, and cannot exclude stronger inequalities restricted to large \(S\) or detailed radius/geometry regimes.

## 6. Actual checks, provenance, and final scope

The independently authored `independent_check.py` actually exited 0. It uses only Fraction arithmetic and endpoint identities, checks the registered probe values and false stronger residual exactly, verifies the sign-separation margins, source masses, energy polynomial coefficients, and rough \(16641/256<256<511\) comparison. It neither imports nor executes author code and does not add a new geometric search. Its actual UTC/runtime/source-hash receipt is `actualreceipt.json`; output is `result.json`.

The author program implements the written one-candidate/one-probe registration and the saved verification agrees with that code. Its recorded plan hash matches the actual plan. The final source manifest reports a successful exact run and stopped-as-registered status, and hashes the four associated source/result files. These records are mutually consistent. The saved package does not independently timestamp immutable preregistration or supply runtime timing; no stronger historical chronology claim is made. Mathematical feasibility and the coefficient refutation are independently proved here regardless of that provenance limit.

No substantive unresolved gap was found in the stated model result. The accepted scope is the exact orthant-energy calibration, a feasible uniform-coefficient obstruction with nonvanishing mass, and the coarse-family comparison. It neither advances the current numerical \(C_8\) record nor restricts all more detailed coupling strategies.

Reproduce this audit's small exact check from the clone root with:

```text
python3 cube-dim-8-bounds/audit/round003/upper/independent_check.py
```

Frozen author hashes, external-source identity checks, audit deliverable hashes, and the unchanged earlier audit manifests are recorded in `manifest.json`. The author base-source and repository state do not imply that unrelated historical blobs were repaired or validated.
