# upper004 spatial continuous interfaces audit

## Verdict

Accepted. No sign or normalization error was found in the real coarea,
log-semiconvexity, vertex-integration, or coordinate ODE interfaces. The
complex-tube assertion is locally holomorphic under the stated componentwise
tube condition: it keeps the real coarea plane and has an integrable
majorant, so no unsupported contour displacement is involved.

These are interfaces, not a completed grid certificate. No new C6 bound
follows until certified vertex intervals and a finite box cover are supplied.

## Actual symbolic check

The independent checker ran successfully:

    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python
      /Users/zhengzhihao/Desktop/T/Math-Research-cube6-literature/cube-dim-6-bounds/audit/upper004/spatial/check_spatial_symbolics.py

It verifies A*A^T=32I, A*(A^T/32)=I, sum(c_v^2)=1/32, and
sum(abs(c_v))=1 for c=A^T e_i/32. The two real log-Cauchy inequalities
reduce to the nonnegative rational functions

    2*y*(y+3)/(1+y)^2
    (y-3)^2/(4*(1+y)^2).

The squared complex denominator difference is

    s*y*(6+2*y-s*(2+y)),

which is nonnegative for 0<=s<1 and y>=0. These are exactly the signs
needed by the source.

## Real coarea and vertex certificates

With R=A^T/32, Rx is perpendicular to w in ker A and norm(Rx)^2 is
norm(x)^2/32. The lower log-second-derivative inequality therefore makes
p_t(x)*exp(norm(x)^2/(32*t^2)) log-convex: after completing the square in
the coarea integral, the remaining w weight is fixed and positive. Holder
preserves log-convexity under integration. For t>=epsilon, an additional
nonnegative quadratic factor and positive time mixture prove the stated
log-convexity of H_epsilon(x)*exp(norm(x)^2/(32*epsilon^2)).

Barycentric interpolation from the 32 face-box vertices has quadratic gap
sum_i(r_i^2-(x_i-c_i)^2), at most sum_i r_i^2. Weighted AM-GM, followed by
integrating barycentric weights, gives the claimed vertex upper bound for
the face integral. The factor and the 1/32 vertex average are correct.

For the ODE side, the exact derivative identity is

    F_ii/F = (sum c_v*(log f_t)')^2 + sum c_v^2*(log f_t)''.

The checked c norms and the upper scalar inequality give
F_ii/F <= 1/t^2 + 1/(128*t^2). Thus
H_epsilon,ii <= (129/(128*epsilon^2))*H_epsilon. Differentiation under
coarea and time integration is justified by the absolute bound
(1+1/16)*t^-2*F and t>=epsilon. The comparison solution for
f''-k^2 f<=0 has the source's lower-bound direction. Positive coordinate
operators commute, so equal-radius box vertices at least L give
H_epsilon >= L/cosh(k*r)^5.

## Complex tube

For x+i*h, w remains real and the integrand uses Rx+i*R*h+w. If every
component abs((R*h)_v) is below t, the checked scalar inequality bounds
each factor by the corresponding real Cauchy factor times a finite
constant. On a compact set strictly inside this tube this is a common
integrable majorant, proving holomorphy by dominated differentiation (or
Morera). For H_epsilon, t>=epsilon makes the factor uniform and the real
H_epsilon integral provides the time majorant.

For one axis, abs((R*s*e_i)_v)=abs(s)/32, giving the power-32 estimate.
For several coordinates the true condition remains
max_v abs((R*h)_v)<epsilon. The source correctly refuses to use the
single-axis 32*epsilon width as a five-coordinate max-norm tube.

## Face reduction for g

A six-dimensional H_epsilon grid is unnecessary for the cube minimum.
Stable-potential homogeneity is G(rho*y)=rho^-5*G(y). For nonzero x in the
unit cube, rho=norm(x)_infty is at most one and y=x/rho lies on its
boundary, hence G(x)>=G(y). Sign and permutation symmetry reduce that
boundary to the positive face. At zero, G(0) is infinite because
p_t(0)=C0*t^-6. Thus a face lower certificate for G>=H_epsilon is enough
for g. It remains separate from the face integral upper certificate.

## Hashes and scope

| artifact | SHA-256 |
|---|---|
| spatial second-order source | e0785ef2226eac8d31bdcb4be2df7b58865fdce7f0d8dce73af1619a149f9489 |
| spatial multiplicative source | b8211c42fbebf036019d786e1bdbff7414fbb464fff8d1ee0e5f3a2a09283d2d |
| independent symbolic check | 3514c6a049c82243f43d08fb156c3349c86e740c485c295a4494f586eabde846 |

This audit creates no numerical vertex interval, box cover, face I bound,
cost bound, or C6 conclusion.
