# R032 joint witness-mass and covering-flow review

2026-10-28. **ACCEPT as finite exact interfaces with the source's limitations.**
The witness-mass primal/packing dual, variational sandwich, adaptive
selection, covering-flow primal, probability/knapsack dual, and margin-only
obstruction are mathematically sound. None alone is an improved explicit
uniform bound for the centered-cube operator.

Source read in full:
`/Users/zhengzhihao/Desktop/T/cube-r030-soft/cube-dim-7-bounds/upper/round032/joint_cover_theory/proof.md`,
SHA256 `6bfd95576ef401f27c23be1eaf3e7b9b19293883f9b2b524a329246ee7c3fdd2`.
Finite families are understood to be nonempty with positive cube radii.
An empty family can be handled separately with zero covered volume, without
forming a 0/0 ratio. The author's separate numerical LP examples were not
rerun or independently certified by this review.

## 1. Witness-mass LP and exact reformulation

Partition the finite union of inflated cubes V_i=aQ_i by all positive-volume
incidence patterns E_sigma. A nonnegative L1 function gives cell masses
m_sigma=lambda^(-1) integral_Esigma |f| with total at most ||f||_1/lambda
and every witness constraint satisfied. Conversely every feasible vector is
realized by its constant cell densities m_sigma/|E_sigma|. Multiplying by
1+epsilon makes all constraints strict, including those originally at equality.
Null cells cannot hide positive mass from an L1 function, and they need not
enter the finite LP.

The LP is feasible, and has positive finite attained optimum: its nonnegative
objective has compact bounded sublevel sets, and each positive witness
constraint bounds total mass away from zero. Finite strong duality gives
the displayed packing dual with sum_{i in sigma}z_i<=1 on every positive
cell. A feasible dual lower-bounds required input mass and therefore has
the correct direction for bounding |U|/M from above. A feasible primal
upper-bound on M alone has the opposite direction and cannot prove that
upper bound.

For K_{n,a}, the left sandwich C_n<=K_{n,a} follows by covering an arbitrary
compact part of a strict maximal level set with finitely many true-witness
cubes shrunk by a. The original f supplies all constraints simultaneously.
No compact support or radius cutoff is needed. For the right sandwich,
the optimal realized g has integral over every V_i at least |V_i|. For each
x in Q_i, V_i lies in Q(x,(a+1)r_i), giving maximal average at least
(a/(a+1))^n. Applying the weak bound at a strictly smaller threshold and
taking that threshold upward gives

```
C_n<=K_{n,a}<=(1+1/a)^n C_n.
```

The known finiteness of C_n ensures a squeeze as a tends to infinity and
inf_{a>1}K_{n,a}=C_n. No monotonicity of K in a is needed or proved. This
is exact but does not evaluate the supremum over every finite geometric
family; a finite sample of families cannot certify it from above.

## 2. Adaptive budgets preserve the real selection prefix

The h_i are positive and fixed before selection, with |V_i|<=h_i<=M_i/lambda.
After an item is selected, its actual increment I_j and p_j=|I_j|/h_j stay
fixed. Nonnegative additions can only destroy eligibility, so largest-radius
selection has nonincreasing radii. This argument remains valid with p_j=0.

For a rejected cube, consider the first prefix that destroys its eligibility.
It was still eligible before each selection in that prefix, hence every
selected radius there is at least its radius. Each prefix term that captures
the rejected center through bQ_j captures the entire rejected cube through
(b+1)Q_j. Their prefix weight exceeds tau, proving P(x)>tau on rejected
volume. This remains true for arbitrary b>0; no t-dependent cap estimate
has been silently used.

The rejected-volume flow is feasible: on each positive cover-incidence cell,
z_sigma,j=d_sigma p_j/P_sigma sums to d_sigma and is at most
d_sigma p_j/tau. The real rejected volume is counted once, including when
the geometric E_j overlap.

## 3. The combined load and the optimization direction

Selected increments sum to selected volume; the flow's c_j sum to rejected
volume. For every feasible flow, the same original witness masses give

```
integral |f| G_z=sum_j [(v_j+c_j)/h_j]M_j
 >=lambda sum_j(v_j+c_j)=lambda |U|.
```

The essential maximum of G_z is precisely the maximum over positive-volume
witness incidence cells. It follows that lambda |U|/||f||_1 is at most every
feasible flow's load bound, hence at most the attained minimum L_*.
The flow polytope is nonempty and compact, and its maximum-of-finitely-many-
linear-functions objective attains its minimum. Thus a feasible flow is a
valid *upper* certificate; a feasible dual lower-bounds the best load
certificate and cannot alone establish an operator upper bound.

The old-product recovery has the correct factors: capacity gives
c_j<=p_j |E_j intersect rejected|/tau, so
(v_j+c_j)/h_j<=p_j[1+|E_j|/(tau h_j)] and
|E_j|/h_j<=((b+1)/a)^n. It is explicitly conditional on a pointwise load
bound for this same adaptive selection. It does not assert that the
original old load proof automatically applies after arbitrary adaptive
budget changes. The coupled flow LP avoids that lossy recovery when a
sharper actual flow certificate is available.

## 4. Exact dual and the margin obstruction

Writing the maximum witness-cell load as a maximum over the finite probability
simplex is exact. The bilinear form on that simplex and the compact feasible
flow polytope permits finite minimax/LP strong duality. For fixed probability,
each rejected cell separates into a fractional-knapsack problem with unit
total allocation and capacities p_j/tau. These capacities sum above one, so
the local problem is feasible; increasing-price filling is its exact solution.
The baseline selected load is sum_j kappa_j v_j, giving precisely (CFD).
The probability is a dual probe of pointwise load and does not replace the
original f or its simultaneous budget constraints.

For any finite geometry, f_epsilon=(1+epsilon)lambda 1_{union V_i} is L1 and
gives all margins exactly 1+epsilon. Thus witness validity cannot provide
a fixed strictly positive uniform margin gap. At h=(1+epsilon)|V|,
p and the threshold both scale by 1/(1+epsilon); deterministic tie choices,
selected increments, flow capacities and cover cells are unchanged, and
the load objective scales by that same factor. This proves the stated
scaling identity even for equality ties.

This only blocks a margin-gap-based proof; it does not block incidence
constraints, flow savings or the later actual-shell deficit method. It is
not a lower bound on C_n or on all possible covering certificates.

## 5. What is and is not proved

All finite interface and bound-direction claims are accepted. A new uniform
operator estimate would still require a packing dual or feasible covering
flow with a suitable load bound for *every* relevant geometric family and
original input, followed by the existing compact-level-set bridge. That
missing uniform theorem is correctly acknowledged by the source.

An actual small Fraction check in `check_small.py` verifies one abstract
two-witness cell-mass primal/dual equality without an optimizer. It is
only an algebra check; it is not a rerun of the author's numerical
examples or a finite substitute for the general LP proof. No search,
n=7 work, MATLAB, Git mutation or R033 change was performed.
