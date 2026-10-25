# A true two-coordinate L1 bound for the same middle tube

Keep precisely the original body32 Cauchy frame and the043 domain

    A={1,2,3}, B={4,5,6}, U0=span(1_A,1_B),
    D={x=(1,y), y in[0,1]^5: dist(x,U0)<rho,
                                  mean_B(x) in[1/4,3/4]},
    rho=1/64, epsilon=1/1000.

Use s for semigroup time, u for the mean_B spatial coordinate, and choose
one fixed index j=4. There is no sum over the three possible indices.
All densities are actual Lebesgue densities of the original common-time
Cauchy parameters. The only new numerical interface is an L1 bound;
the pointwise043 constants are read from their complete saved intervals.

## 1. Exact two-coordinate marginal

Partition the32 body directions according to their j-th sign. There
are2^4=16 directions in each block. Let U and V be the sums of their
independent Cauchy(s) parameters. Then U,V are independent Cauchy(16s),
and

    (X1,Xj)=(U+V,U-V).

The absolute inverse Jacobian is1/2. With f_a(t)=a/[pi(a^2+t^2)], the
two-coordinate marginal is therefore

    p2_s(1,z)=(1/2)f_(16s)((1+z)/2)f_(16s)((1-z)/2)
      =2048s^2/[pi^2(((1+z)^2+1024s^2)((1-z)^2+1024s^2))]. (1)

The marginal integral identity holds at every fixed(1,z), not merely at
almost every pair of retained coordinates. To verify this directly,
choose a six-column body basis whose first two columns have j signs+1
and-1. Their2x2 matrix in coordinates1,j has determinant-2. The remaining
four basis parameters are free. In the positive basis formula for p_s,
solve the first two parameters using the retained coordinates, and change
the other four physical coordinates into those four free basis parameters.
The Jacobian is |det(basis)|/2; it cancels the original1/|det(basis)|,
leaving1/2. Tonelli over the other30 Cauchy parameters factors into the
two16-fold convolutions in(1). This is an identity of explicit positive
integrals for every retained coordinate pair.

One explicit basis is ones, followed by its sign flips at coordinates
4,2,3,5,6. Its full determinant has absolute value32, the retained pivot
determinant has absolute value2, and the four-variable Schur Jacobian has
absolute value16. The exact structural check verifies these identities.

## 2. The positive-coordinate slice factor is exactly1/16

The full symbol, and hence p_s, is invariant under each of the four
coordinate reflections outside{1,j}. Such reflections permute the body
directions and leave the two retained coordinates fixed. The16 orthants
of those four coordinates have equal integrals; their boundaries have
four-dimensional Lebesgue measure zero. Thus for every z,

    integral_(remaining4 coordinates>=0) p_s(1,z,rest) drest
       =p2_s(1,z)/16.                                    (2)

There is no additional factor for coordinate1, whose value is fixed+1,
nor for coordinatej, whose value is fixed z. No reflection symmetry is
claimed for a parameter-restricted branch density.

## 3. The true projection range of D

For the B block write x_B=u1_B+q, sum q_i=0. Orthogonal projection to U0
gives dist(x,U0)^2=|e_A|^2+|q|^2, so |q|<rho on D. The norm of the
coordinate functional on the three-dimensional zero-sum space is
sqrt(2/3): q_j=(e_j-(1/3)1_B).q. Therefore

    z=xj in[L,U],
    L=1/4-sqrt(2/3)rho,
    U=3/4+sqrt(2/3)rho.                                  (3)

Here0<L<U<1 since sqrt(2/3)rho<1/64<1/4. This bound uses the actual
B-block mean and residual, not a sup-norm proxy for the tube. The actual
D is contained in the set with this z range and the remaining four
coordinates nonnegative. Positivity permits extending those four box
constraints to their entire positive orthant. Each point is integrated
once in the chosen fixed coordinate; no overlapping tube union occurs.

Combining(1)-(3) and dropping positive time terms in the two denominators,

    integral_D p_s dy <=C_joint s^2,  for every s>0,

    C_joint=(128/pi^2) integral_L^U (1-z^2)^-2 dz.         (4)

The constant128 comes from2048/16. It does not require a1/32 symmetry
factor for five coordinates after one of them has already been retained.
The elementary primitive is

    F(z)=z/[2(1-z^2)]+(1/4)log((1+z)/(1-z)),
    F'(z)=(1-z^2)^-2 on(0,1),
    C_joint=(128/pi^2)(F(U)-F(L)).                        (5)

## 4. Transfer to the actual043 branches

Use the exact043 event partition, with v=ones and w=(+++---):

    E={|a_v|>=1/2,|a_w|>=1/16},
    E_v={|a_v|<1/2},
    E_w={|a_v|>=1/2,|a_w|<1/16}.

Let P_s be the density of E and B_s the density of the actual disjoint
complement E_v union E_w. Then P_s+B_s=p_s and0<=P_s,B_s<=p_s. The043
proof supplies pointwise bounds on this same D,

    P_s<=D_big s^-2, B_s<=M_bad,

where D_big is the saved043 C2 and M_bad is its saved M_combined. Their
source expressions, normalizations, complete complementary events and
common-time laws are unchanged. The enlarged events used to bound
M_bad are not asserted disjoint here; only the actual B_s is used in
the L1 domination. From(4), separately,

    ||P_s||_L1(D)<=C_joint s^2,
    ||B_s||_L1(D)<=C_joint s^2.                            (6)

These are domination bounds from the full density. We have not applied
the1/16 orthant symmetry to either branch.

For q=6/5, interpolation has exponents1/6 and5/6. Thus

    ||P_s||_Lq(D)<=D_big^(1/6)C_joint^(5/6)s^(4/3),
    ||B_s||_Lq(D)<=M_bad^(1/6)C_joint^(5/6)s^(5/3).

Minkowski on positive time truncations and monotone convergence yield

    ||integral_0^epsilon p_s ds||_Lq(D)
      <=(3/7)D_big^(1/6)C_joint^(5/6)epsilon^(7/3)
        +(3/8)M_bad^(1/6)C_joint^(5/6)epsilon^(8/3).        (7)

All actual event mass is included, including every cancellation and
other large parameter already paid by043. The only new improvement is
the true O(s^2) joint marginal L1 interface.

## 5. Scope of the numerical certificate

Read D_big,M_bad and the old043 short-budget total from their full lower
and upper dyadic endpoints. Do not reconstruct any value from a printed
decimal and do not reevaluate the043 pointwise or short-budget formulas.
Evaluate only the fixed L,U,C_joint and the new formula(7). Comparison
with the saved old budget is a comparison of certified upper-budget
formulas, not a lower bound on the actual norm or an actual improvement
factor of that norm. No larger time range, different tube, other plane,
endpoint, whole-face integral, global g or C6 result is certified here.
