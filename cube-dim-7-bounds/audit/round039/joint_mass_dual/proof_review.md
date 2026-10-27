# R039 common-witness mass dual and clusters: independent audit

2026-10-28. **ACCEPT with the stated interface scope.** The common-star
payment theorem holds for every common witness, the fractional cluster
certificate pays representative reuse correctly, and the finite
arrangement-cell LP preserves all certificates of that form. Both
counterexample diagnostics are valid in their stated domains. There is
no new numerical or uniform C_n bound until a uniform representative-
overlap certificate is constructed for all relevant histories.

## 1. Locked source and review scope

Source read independently in full and checked against its actual hash:
`/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round039/joint_mass_dual/PROOF.md`,
SHA256 `3fa6117e186343466f3a1cbeb18ed4084a12f4b303e966f14480687a1edd24f6`.

The review checks the entire finite mapping, mass dual, geometric payment,
cover factor, arrangement-column reduction, general-input/all-scale bridge
and counterexamples. It performs no production computation, numerical
search, n=7 calculation, MATLAB or Git operation. Only
`audit/round039/joint_mass_dual/` is written; frozen earlier rounds remain
unchanged. A finite family is understood to be nonempty with positive
radii. An empty selected family has zero selected volume and needs no
mass ratio or normalization by ||W||.

## 2. Exact witness-mass LP and the old particular dual

For a fixed original selected history, the actual increments I_i are
disjoint and satisfy I_i subset Q_i. Their sum U is the volume of the
selected union. Every normalized original common witness f/lambda is
feasible for the mass problem, so its cost mfrak is at most ||f||/lambda.
Since U>=0 and mfrak>0,

```
lambda U/||f||<=U/mfrak.
```

This is an upper payment ratio for selected volume, not a C_n lower bound
or a claim that the entire selected union lies in a strict maximal level
set at the same level.

For the finite LP, any mass outside the union of stars can first be
removed without affecting a constraint or increasing cost. The positive-
volume coordinate arrangement cells then carry every relevant L1 mass.
Null faces cannot carry positive L1 mass. Conversely constant densities
on these cells realize every feasible cell-mass vector. The LP is feasible,
its nonnegative bounded objective sublevels are compact, and its optimum
is positive because of any one positive star-volume constraint.
Finite strong duality gives exactly the stated packing dual, with an
almost-everywhere overlap constraint. Strict-witness infima have the same
value by multiplication by 1+epsilon.

For W=sum_i(v_i/V_i)1_{S_i}, the vector z_i=(v_i/V_i)/||W||_infty is
feasible and has value U/||W||_infty. The essential supremum is positive
for a nonempty history because its first increment has positive volume.
This normalization is only one feasible dual; the optimal common-mass
dual can be strictly stronger. No independently assigned witness masses
have been introduced.

## 3. Whole common-star payment, for every common f

All stars are open. A finite nonempty intersection is therefore an open
set of positive volume. If x lies there and R is a largest small radius,
every point of every Q_i lies in Q(x,(a+1)R). Actual increment disjointness
then gives

```
U<= [2(a+1)R]^n=c_a V_j, c_a=((a+1)/a)^n.
```

The same representative j is a true witness for the original f, so
lambda U<=c_a lambda V_j<=c_a integral_{S_j}f<=c_a||f||.
This payment is valid for *every* common f and does not require its mass
to be concentrated in the intersection.

The mass optimum equals V_j: its j-th constraint gives the lower bound,
and a density of total mass V_j supported in the common open intersection
is feasible for every i since V_i<=V_j. Equivalently the single packing
dual z_j=1 already supplies this optimum.

Thus the entire R037 common-star high-load family has selected-volume
payment at most c_a, despite large W at a common point. This is stronger
than the earlier particular-spike diagnostic. It still concerns the
selected union of a family whose stars have a common intersection; it
does not claim that arbitrary global selected families have such a pivot
or that the general maximal constant is one.

## 4. Fractional clusters pay overlap with the same common input

For any cluster J, its original increments remain disjoint even if J
is not a prefix or contiguous subhistory. Their union lies in the same
pivot-centered box, so sum_{i in J}v_i<=c_a V_{j(J)}. Recomputing increments
inside J is neither necessary nor allowed by this argument.

Nonnegative theta_J with item coverage at least one gives

```
lambda U
 <=lambda sum_J theta_J sum_{i in J}v_i
 <=c_a sum_J theta_J lambda V_{j(J)}
 <=c_a integral f sum_J theta_J 1_{S_j(J)}
 <=c_a K ||f||.
```

Every inequality uses the same original common f. The representative
overlap constraint is required almost everywhere, which is sufficient
for arbitrary L1 f; it cannot be dropped merely because each cluster
has a good individual payment.

For the original threshold algorithm, the accepted first-rejection
prefix cover step bounds all input-union volume by
[1+(beta+1/a)^n/tau] times selected-union volume. Applying the above
payment yields the source's conditional operator bound, provided K_*
bounds a valid cluster certificate for *every* selected history from
the finite true-witness covering problem. The ordinary compact-subset
and inner-regularity argument then covers general L1 and all positive
common radii. No bounded-radius or compact-support assumption is added.

This bridge inherits the original algorithm's allowed parameters; it
is not an assertion that an unconstructed global K_* has been bounded.
A finite collection of solved instances would not control the supremum
over all histories, nor justify a scalar a-to-infinity optimization
without a corresponding bound on K_*(a).

## 5. Arrangement columns are complete while preserving the representative

For finite open axis-parallel cubes, pairwise intersection is equivalent
to total intersection: in each coordinate, the largest lower endpoint
is strictly below the smallest upper endpoint. This proves the stated
Helly-number-two test with strict inequalities.

More importantly, take any admissible cluster J with chosen maximum-
radius member j. Its common intersection is open, so choose a pivot
off every star-face hyperplane. The positive-volume arrangement cell C
containing that pivot lies in every star of J. Hence j belongs to T_C
and

```
J subset {i in T_C:r_i<=r_j}=J_{C,j}.
```

Enlarging the cluster this way preserves the representative exactly,
keeps its paid star S_j unchanged, and only increases item coverage.
Several upgraded columns can be merged by adding their coefficients.
This proves that the finite columns lose no certificate in Theorem 2.
In contrast, extending to arbitrary maximal cliques might introduce a
larger representative or a different paid star and has no such guarantee.

The representative-overlap constraint for cell D is correctly indexed
by j in T_D, independent of the pivot cell C. This is the load of S_j,
not of the pivot cell. Every feasible nonnegative LP primal thus gives
the right direction for an upper payment bound. No dual lower bound or
floating-point optimizer value has been substituted for a feasible
certificate. Positive-volume cells suffice because all paid indicators
are constant almost everywhere.

The continuous pivot-measure notation is consistent with this finite
form. Any pivot induces a cluster of stars containing it whose radii
are no larger than the representative. Boundary pivot atoms cause no
loss: each such finite cluster still has an open common intersection
and can be enlarged using an off-face cell as above. The finite form
does not depend on a measurable choice over uncountably many pivots.

## 6. The two diagnostics and their precise quantifiers

For equal small radii r and centers +2r e1 and -2r e1, a>2 gives a
positive-volume star intersection. A spike supported there has integral
equal to its whole mass over each star, so separately adding the two
representative payments charges 2||f||. This exactly exhibits the missing
overlap factor. One combined cluster legitimately has K=1; the example
is not an operator lower bound.

For the original-center greedy rule, fix any N and then choose a large
enough that a>max(3,2t) and (N-1)a^(-n)<=tau. With R_i=3^(-i),
y_i=(3/4)R_i e1 and r_i=R_i/a:

- Neighboring small cubes are disjoint for a>8/3, which is implied.
  Their actual increments equal themselves and have weights a^(-n).
- All stars contain 0 in their interiors. Every preceding bQ_i captures
  each later center, but its whole previous load is at most
  (N-1)a^(-n)<=tau. Largest-radius selection therefore legally selects
  every item in radius order.
- The largest remaining star's own center belongs to no later star:
  its coordinate (3/4)R_i exceeds every later right endpoint, at most
  (7/12)R_i. The proposed rule removes only itself each time and produces
  N singleton representatives.
- Their overlap is N on a *positive-volume neighborhood* of 0, not only
  at one point. The genuine free pivot 0 yields one cluster and K=1.

This gives arbitrarily large loss over the allowed parameter families.
It does not assert that N can be unbounded with n,a,tau all fixed:
the displayed eligibility condition is deliberately preserved by the
choice of a. A common strict L1 witness can also be supplied on the
finite open star intersection, so witness consistency does not rescue
this particular greedy rule.

## 7. Final interpretation

All finite common-mass mappings, payment inequalities, arrangement LP
directions and diagnostics are accepted. They preserve the actual
increments and original witnesses, and the conditional cover bridge
preserves the original general-input/all-scale scope. The unresolved
step is a uniform small K_* construction, using geometry/history while
paying representative overlap. This source gives no new numerical C_n,
no C_n lower bound and no global optimality claim for the cluster method.
