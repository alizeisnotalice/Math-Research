# Round 019: one-axis alternating-mass first-order diagnostic

No new C6 bound is obtained. Translation symmetry rigorously removes a
**differentiable** first-order improvement at epsilon=0, but it does not
exclude a positive absolute-value cusp. Concavity of each mass-to-radius
map does not imply concavity of the finite witness union. The sign of the
actual round008 family's one-sided union derivative remains pending an
exact exposed-boundary calculation. No such union or search was run here.

## Fixed, precisely defined inherited subfamily

Keep the accepted round008 spacings, with positive u^6=2,v^6=3:

```
a1=1+u^4/2,
a2=u^2*v^2/3+u^2*v^4/6,
a3=a4=a5=a6=(1+u)/2.
```

The fixed label set L is **all 127 raw block-length vectors** exported in
the accepted `lower/round008/spacing_chamber/chamber_interface.json`, plus
their ordinary coordinate reflections and receiver-cell translations.
No 96-box dominance pruning is inherited: changing masses could change
containment. No additional block labels or previously empty boxes are
added. Thus the following conclusions apply to this exact finite subfamily,
not to every cube witnessing the alternating-mass source. Length-5/6 and
old baselines are not re-evaluated.

At site a*k put positive mass

```
w_epsilon(k)=1+epsilon*(-1)^k1,    |epsilon|<1.
```

The true periods are 2a1 in axis one and aj in the other five axes. A
half-open period cell contains two sites with masses 1+epsilon and
1-epsilon, hence **total source mass exactly two**. This mass identity
alone says nothing about the level-set volume.

Use two adjacent receiver cells centered at 0 and a1*e1, each with local
coordinates in C=product[-aj/2,aj/2]. Together they form one true double
cell, up to null boundaries. Let H=product[0,aj/2]. In a receiver cell
whose first-coordinate index parity is s in {0,1}, the positive-orthant
block of length ell has local starts

```
qj=-floor(ellj/2)    if ellj is odd,
qj=1-ellj/2         if ellj is even.
```

Write P=product_{j>1}ellj and

```
S(q,l)=sum_{k=q}^{q+l-1}(-1)^k
      =0 (l even), or (-1)^q (l odd).
M_s(ell,epsilon)=P*(ell1+epsilon*(-1)^s*S(q1,ell1)).
```

This counts the selected atoms directly, including blocks crossing a
period boundary. Reflection about the receiver-cell center sends local
k1 to -k1. It preserves (-1)^(s+k1), and therefore the selected mass.
Equivalently the reflected block start is -q1-ell1+1 and has the same
parity sum. Reflections in the other axes also preserve mass. Thus each
cell still has 64 admissible sign-reflected half-cell witness copies.

## Threshold and common-radius interface

For any fixed lambda>0, a block uses the **one common radius**

```
r_s=(M_s/lambda)^(1/6)/2,
```

in every coordinate. With local qj as above, its center box is

```
B_s(ell,epsilon,lambda)
 = product_j([(qj+ellj-1)*aj-r_s, qj*aj+r_s]) intersect H.
```

Every such cube contains the selected block and has average at least
lambda. Null/reversed intervals give empty boxes. For strict weak levels,
use a level below lambda or inward boxes and a smaller threshold, as in
the accepted finite-L1 transfer. This diagnostic fixes **lambda=1** and
varies only epsilon; varying lambda would be an additional parameter.
Keeping the old radii while assuming every old block still witnesses
level one would be incorrect when its selected mass decreases.

Let sigma=S(q1,ell1), m=P*ell1, and r0=m^(1/6)/2. Then

```
r_s(epsilon)=r0*(1+epsilon*(-1)^s*sigma/ell1)^(1/6),
r_s'(0)=r0*(-1)^s*sigma/(6*ell1),
r_s''(0)=-5*r0*sigma^2/(36*ell1^2).
```

In particular even ell1 has unchanged radius, ell1=1 expands in the even
receiver cell, and ell1=3 shrinks there. The odd receiver cell interchanges
these changes. Radius perturbation acts on **all six slab coordinates**;
it is not a one-dimensional interval-length perturbation only.

## Exact symmetry and the derivative criterion

Define V(epsilon)=measure of the union of the 127 B_0 boxes in H at
lambda=1. The neighboring cell has half-cell union V(-epsilon). Reflection
copies occupy disjoint half-cell interiors within each receiver cell,
so the true double-cell witness ratio is exactly

```
R(epsilon)=[64*V(epsilon)+64*V(-epsilon)]/2
          =32*[V(epsilon)+V(-epsilon)].
```

At zero this is the accepted ordinary-lattice witness ratio 64V(0).
Translation by a1 in axis one changes epsilon to -epsilon, agreeing with
the direct double-cell formula. Hence **R is even**.

The finite boxes have algebraic/analytic radius functions near zero,
followed by finitely many min/max, interval-clipping and finite-union
operations. Endpoint differences either agree identically or have a
definite sign on each sufficiently small punctured side of zero. The
union is therefore piecewise analytic near zero and has finite one-sided
derivatives. Writing V'_+(0),V'_-(0), the exact criterion is

```
R'_+(0)=32*[V'_+(0)-V'_-(0)].
```

If V is differentiable at zero, the right side is zero and there is no
first-order gain. More generally:

* an upward slope jump of V gives a positive `|epsilon|` gain;
* a downward slope jump gives a negative `|epsilon|` loss;
* equal one-sided slopes leave the second-order sign unresolved.

Evenness alone cannot choose among these possibilities. If one could
prove V concave on a neighborhood, it would imply R(epsilon)<=R(0),
but that concavity has not been established for this subfamily.

## Why radius concavity is insufficient

For a simple geometric counterexample to the proposed inference, let
I_+(epsilon)=[0,(1+epsilon)^(1/6)] and
I_-(epsilon)=[0,(1-epsilon)^(1/6)]. Both individual lengths are concave,
but their union length is

```
max((1+epsilon)^(1/6),(1-epsilon)^(1/6))
 =(1+|epsilon|)^(1/6),
```

which has a positive first-order cusp at zero. Multiplying by five fixed
positive intervals gives the same counterexample in six dimensions.
This is only a counterexample to the **general logical inference** from
individual radius concavity to union concavity; it is not claimed to be
an admissible mass assignment within the fixed round008 block family.

Single-box volumes also need not be concave in mass. For an unclipped
box with six equal positive source-span offsets c, its volume is
`(r-c)^6`. As a function of m=(2r)^6 in the positive-width region, this is
`2^-6*(m^(1/6)-2c)^6`, whose derivative is
`2^-6*(1-2c*m^(-1/6))^5` and whose second derivative is positive for c>0.
Thus neither the union operation nor even the mass-to-box-volume step
allows a general Jensen exclusion based on r''<0.

## A specific current degeneracy does not supply a positive mechanism

The accepted spacing event has an axis-two tie between the blocks

```
A=(2,3,2,3,3,3),  m_A=324,  r_A=u^2*v^4/2,
B=(3,4,3,4,4,4),  m_B=2304, r_B=u^2*v^2,
3*a2=r_A+r_B.
```

For A, ell1=2, so its radius is fixed. For B, ell1=3 and q1=-1, so
`r_B(epsilon)=r_B*(1-epsilon/3)^(1/6)` in the even receiver cell. The
axis-two gap between A's upper face and B's lower face is

```
g(epsilon)=3*a2-r_A-r_B(epsilon),
g(0)=0,  g'(0)=r_B/18>0.
```

This confirms a real first-order raw-endpoint splitting, rather than
arguing solely from fixed total mass. However, A and B have strictly
disjoint **axis-one** projections at zero: A starts at a1-r_A, while B
ends at r_B-a1, and

```
2*a1 > r_A+r_B.
```

A short exact proof writes x=u^2=2^(1/3), y=v^2=3^(1/3). Cubing positive
rationals gives x<4/3, y<3/2, and x^2>3/2. Consequently
`r_A+r_B=x*y*(1+y/2)<7/2<2+x^2=2*a1`.
The separation persists for sufficiently small epsilon by continuity.
The corresponding tail permutations have the same axis-one separation.
Therefore this known A/B ordering tie cannot itself create an overlap
cusp between those blocks. It is not evidence of a positive union gain.
This argument does not inventory every endpoint coincidence or exposed
face of the full 127-box union.

## Remaining gap and one minimal future candidate

There is no certified sign for `V'_+(0)-V'_-(0)` in this exact subfamily.
The missing object is the exposed union-boundary/overlap incidence at
epsilon=0, including coincident endpoint functions that split when their
block parity responses differ. Raw endpoint order or total mass alone
does not determine that object. The accepted spacing chamber's grouping
of identical affine functions cannot automatically be reused for a new
mass direction: functions equal along the spacing direction can acquire
different epsilon velocities.

The smallest useful next acceptance object would be **one** exact pair
of one-sided union first derivatives on the unchanged 127 labels, with
lexicographic endpoint value/velocity ordering and an exposed-coverage
certificate. It would decide the cusp without a new nonlinear parameter
scan; if the slopes agree it would still leave the quadratic term open.

If main instead authorizes a single finite nonzero candidate, nominate
axis one at the same exact spacings, **epsilon=1/1024**, lambda=1, and
the same 127 labels in both parity phases. Its acceptance object is the
joint exact ratio `32*(V(1/1024)+V(-1/1024))`, with common radii from
the rational block masses above, compared against the **complete**
accepted round008 ideal upper endpoint. Roots may be bounded by rational
integer-power certificates; no assumption that the new radii remain in
Q(u,v) is warranted. Either improvement or exclusion of that one candidate
would require a new accurate union/derivative calculation. Neither has
been authorized or run in this diagnostic. No automatic certification or
new batch is started.

## Finite resource and evidence closure

Registration was written and hashed before checks. Its SHA-256 is
`ed16731637533491ddb2903c441e3dd77f06d72f0b760c18e1d2fefcb485fffe`;
it records exact hashes of the locked ledger, state, interface and static
candidate result. `parity_counts.py` performs exactly **12** Fraction
cases: starts -1,0,1 and lengths 1,2,3,4, epsilon=1/32, other-axis mass 6.
It verifies the direct weighted sum, formula, reflection invariance,
neighbor-phase/negative-epsilon equivalence, and constant paired mass.
This is a count-interface check, not a new geometric candidate.

No full union, optimization scan, old lower baseline replay, length
extension, MATLAB call, agent expansion, staging, commit or push occurs.
Only Python standard-library integer/Fraction arithmetic is required.
Actual PID, CPU time, source hash and exit status are recorded with the
check output. All writes are in this owned directory. The current
certified C6 lower bound and its complete interval remain unchanged.

Actual count run: PID **34698**, exec chunk **221f32**, exit **0**,
CPU **0.000278 seconds**. The source SHA-256 is
`646c52fab28979ce62756bf5b0277173861e6b21d70eb29daf791563031e0874`.
The actual exit, output hashes, resource counts and final stop are also
saved in `execution_receipt.json`. The five retained files are this report,
registration, the small checker, its exact results, and that receipt.

**Stopped after this one low-resource diagnostic; await main's next batch.**
