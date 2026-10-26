# 080: explicit finite-six allocation in the Stein–Strömberg shell argument

This is a quantitative diagnosis of one specified family of covering
certificates, not an exclusion of covering arguments or a new bound for C6.
The original source is Stein–Strömberg (1983), Lemma and Theorem 1, print
pp.261–263. The proof below replaces its expansion n by beta>1 and uses exactly
m equal logarithmic radius shells, m a positive integer. No constants are
inferred from a displayed O(n log n) assertion.

For n=6 the resulting allocated weak constant is

    K(beta,m)=[1+((beta+2)/beta)^6]
              [1+m((beta+1)/beta)^6 beta^(6/m)].                 (1)

All members of this explicitly allocated family satisfy

    K(beta,m)>1602/25=64.08>64.                                 (2)

Thus none of these allocations improves the existing 64 certificate. Equation
(2) is not a lower bound on C6, the true best covering constant, or any improved
version of the original selection scheme.

## Complete selection and multiplicity argument

Let B be an open bounded convex symmetric body, and use its induced norm.
For a finite family of translates of dilates, write B_j^*, B_j^**, B_j^*** for
the concentric expansions by beta, beta+1, beta+2. Greedily select a remaining
ball of maximal radius among those whose center y_k satisfies

    sum_(j<k) [|I_j|/|B_j^*|] 1_(B_j^**)(y_k)<=1,

where I_j=B_j minus the previously selected balls. Selected radii are
nonincreasing: an ineligible ball cannot become eligible after adding
nonnegative terms. Each selected ball is removed from the finite list, so the
procedure stops. The I_j are disjoint.

For an unselected ball of radius r, the sum of weights of earlier selected
balls of radii at least r, evaluated at its center using **, exceeds 1.
Otherwise it would have been eligible before any smaller radius was chosen.
If x is in this unselected ball and its center is in B_j^**, then x is in
B_j^***, since r<=r_j. Integrating this covering inequality over unselected
balls gives

    |union original balls|
      <=[1+((beta+2)/beta)^6] sum_j |I_j|.                       (3)

Now fix x in at least one B_j^*. Choose a selected containing ball of smallest
radius, translate x to zero and normalize that radius to 1. Every contributing
radius is at least 1. Let its selected center be y_k, so |y_k|_B<beta.
For every contributing radius r_j>=beta, triangle inequality gives
|y_j-y_k|_B<beta*r_j+beta<=(beta+1)r_j. These j occur before k because beta>1.
Their total weight is at most 1 by the selection condition at y_k.

For contributing radii in a shell [a,b), 1<=a<b<=beta, the relevant disjoint
increments lie in the ball of radius (beta+1)b about zero, whereas each
denominator |B_j^*| is at least |B|(beta*a)^6. Their sum is therefore at most
((beta+1)/beta)^6 (b/a)^6. Taking m shells with ratio beta^(1/m), covering
[1,beta), proves the pointwise multiplicity bound

    sum_j [|I_j|/|B_j^*|] 1_(B_j^*)(x)
       <=1+m((beta+1)/beta)^6 beta^(6/m).                       (4)

Open boundaries or equal radii cause no omitted intervals: assign shell
endpoints consistently; the outer part includes r_j=beta. The argument uses
only volume dilation and disjoint increments, and remains valid for cubes.

## General L1 input and every scale

Fix nonnegative f in L1 and lambda>0. For every x in a compact subset K of
{M_B f>lambda}, choose a witnessing radius R_x and let the small ball have
radius R_x/beta centered at x. Its * expansion is the original common-scale
witness. A finite subfamily of these open small balls covers K. Apply (3),
use |B_j^*|<lambda^(-1) integral_(B_j^*) f, and then (4). This gives
lambda|K|<=K(beta,m)||f||1. Supremizing over compact subsets gives the full
weak inequality. Use |f| for signed inputs. Cube boundaries have zero Lebesgue
measure, and every positive witness radius was allowed. Thus this is a genuine
upper certificate for the requested operator, albeit an uncompetitive one.

## Uniform exact exclusion of the allocated formula

For a=6 log beta>0, m exp(a/m)>=e*a: exp(u)/u has its minimum e at u=1.
Consequently the second factor in (1) is at least
1+6e log beta ((beta+1)/beta)^6. The following elementary rational strict
inequalities suffice: e>8/3, log 2>2/3, log 7>97/50. The first follows from
the first four terms of exp(1) plus its positive tail; the second from the
first term of 2 atanh(1/3) plus its positive tail. For the third, the first
eight terms of 2 atanh(3/4) already exceed 97/50. All series converge with
positive terms.

* If 1<beta<=2, the first factor is at least 65 and the second is greater
  than 1; hence the product exceeds 65>1602/25.
* If 2<beta<=7, the first factor is at least 1+(9/7)^6 and the second exceeds
  1+6(8/3)(2/3)=35/3. Their rational product exceeds 1602/25.
* If beta>=7, the first factor exceeds 2 and the second exceeds
  1+6(8/3)(97/50)=801/25. Their product exceeds 1602/25.

These overlapping cases cover all beta>1 and every integer m>=1. The code
checks only the finite exact rational facts; it does not replace the analytical
selection or minimization proof, and no beta/m search or old baseline is run.

## Status

The main derivation and exact checks are new internal work. The public source
provides the selection argument, not the numerical family-wide statement (2).
No publication or exhaustive novelty claim is made. Public and internal C6
thresholds are unchanged. An independent review is still required before
accepting this diagnostic in the final ledger. Better overlaps, a different
selection rule, or cube-specific geometry are outside the family (1).
