# R032 multiscale-theory review

2026-10-28. **ACCEPT with the source's stated scope restrictions.** The finite
two-region guard, exact one-dimensional shell interfaces, single-box
obstructions, and shared-deficit fractional-charging theorem are correct.
They are useful interfaces and scoped obstructions, not an explicit new
uniform operator upper bound.

Source read in full:
`/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round032/multiscale_theory/REPORT.md`,
SHA256 `ee433a4037b672f11fd27a32bb3460f6f3836647985be6b513c0e5d4d9d67420`.
The already accepted actual-increment and finite covering framework remains
an explicit dependency. This review ran no search or LP solver and does not
claim to have re-executed the author's source checks.

## 1. The actual two-region implication

The measure mu has exactly the original uncaptured weight as total mass;
the common scaling Jacobian, actual increment disjointness and r_j>=m(z)
give mu<=lambda_a. For A=(-s,-r)^q times (-h,h)^(n-q), the assumptions
1<r<s<t, h<=s and r>beta*s-1 imply A subset the limiting corner cap and
also the finite-a corner envelope. No actual item radius is replaced by
the spatial radius m(w).

If actual mu(A)>0, choose one actual item j and w in its increment intersect A.
Its uncapturedness gives c*r_j<=||u-w||<=alpha+s, hence r_j<=R under the
source's finite guard. Every actual item i represented in the common P has
r_i>=m(z)>R>=r_j; strict radius priority makes i<j. The center displacement
is strictly smaller than [kappa+2/(a+1)]r_i<=alpha*beta*r_i. Thus all such
i belong to this same j's eligibility test, giving mu(P)<=tau. The alternate
mu(A)=0 branch loses all lambda_a(A). Applying both branches to the single
nonnegative deficit measure d=lambda_a-mu proves (2.2).

For fixed R>(1+s)/beta and 1<kappa<beta, both finite guards hold for every
sufficiently large a, by their strict limiting margins. If R>=t, P is empty;
the source correctly excludes deriving a positive cut from that case.

## 2. Shell formulas and boundary care

The cumulative volume

```
V_A(rho)=[(min(rho,s)-r)_+]^q [2 min(rho,h)]^(n-q)
```

is exact. Because it is continuous and piecewise polynomial, its Stieltjes
integral against rho^(-n) has no endpoint atoms. Since r>1 and A lies in the
cap, this gives precisely the stated weighted mass.

The common-capture region in coordinate i, with A interval (l_i,u_i), is
(u_i-kappa*rho,l_i+kappa*rho) for the source's strict supremum condition.
The lengths L_i are its intersections with the full cube, and C_i intersect
it additionally with the corner-captured interval. For 1<=rho<t, the
corner-captured box has all positive maximal faces and no negative ones,
since rho<alpha/(c-1). Therefore a face decomposition gives exactly

```
sum_i [(F_i^-+F_i^+) product_{j!=i}L_j
       -F_i^+ product_{j!=i}C_j].
```

The normalization is 2^(-n)rho^(-n) d rho, not Euclidean surface measure.
Ties of maximal coordinates have zero n-volume; a face-indicator equality
may be a full face at one radius and must not be called zero face area.
Its exceptional radius has zero d rho measure. All coefficient changes
occur at finitely many linear breakpoints, as claimed.

Strict capture of every individual w in an open A can include a boundary
where the supremum over w equals kappa*rho. The source's P explicitly uses
the strict supremum condition, whose interval formula above is exact;
the potentially larger every-w set is not required. These boundary
distinctions do not affect the stated weighted integral or its safe guard.

## 3. Single-box obstruction is correctly restricted

If h<=1 and q<n, the radial integral extends upward to infinity to yield

```
lambda(A)<=2^(-q)(h/r)^(n-q)/binom(n-1,q).
```

Indeed the substitution z=(rho-r)/rho gives the integer Beta integral
B(q,n-q), and multiplying by q gives 1/binom(n-1,q). The inequalities
binom(n-1,q)>= (n-1)/q and 2^q>=2q imply the further bound 1/[2(n-1)].
The q=n case instead has the full-negative-orthant bound
2^(-n)n log(s/r); the source does not incorrectly apply the former Beta
integral at its divergent endpoint n=q.

For h>1, no transverse coordinate can provide a common captured maximal
face: it would require (kappa-1)rho>h, while rho<t and kappa<=beta imply
(kappa-1)rho<1. Nor can a positive face occur in a negative-coordinate
direction, since that would require (kappa-1)rho>s>1. Thus only q negative
faces remain, giving lambda(P)<=q log(t/R)/2. Also A occupies at most the
2^(-q) sign fraction of its full radial shell. If lambda(P)>tau, then
q>2tau/log t, proving the advertised exponential bound on the cut budget.

The no-growth conclusion is consequently correct only for this fixed-axis
single-box family on the target scale t comparable to n/log n and tau
comparable to n/log n. It does not rule out arbitrary t,tau, all outer-shell
sets, multiple compensated cuts, or altered covering products. In
particular the later R033 good shell set is not one of these boxes, so
its growth deficit is compatible with this obstruction.

## 4. Fractional charging and missing global improvement

For each valid cut, either d(A_i)=lambda_a(A_i)=a_i or d(P_i)>=p_i. Thus
d(A_i)/a_i+d(P_i)/p_i>=1. Multiply by nonnegative gamma_i, sum, and integrate
the pointwise capacity inequality against the *same* d<=lambda_a. This
proves sum gamma_i<=d(total), exactly the source's bound direction. Positive
lower mass bounds in the denominators remain valid: they increase the
required pointwise charge, and in each branch the corresponding deficit
ratio is still at least one.

All shared-P and overlapping-A costs are included in this one pointwise
capacity. If one actual j catches several guarded blocks, their union is
paid by one eligibility threshold; separate thresholds may not simply be
subtracted or added. The source correctly requires one joint union rather
than summing overlapping masses.

No certified geometric collection with growing sum gamma_i is supplied.
The theorem is therefore a legal certificate interface, not a proof of a
new C_n value or a leading-coefficient improvement. Any later numerical
candidate must still verify guards, mass directions and pointwise capacity
over every geometric cell.

## 5. Actual bounded checks

`check_small.py` in this audit directory actually passed exact integer
versions of the Beta/combinatorial inequalities at five fixed pairs,
including q=n-1. These are finite transcription checks; the general result
is established analytically above. No parameter search, n=7 computation,
MATLAB, Git mutation or change to frozen R033 was performed.
