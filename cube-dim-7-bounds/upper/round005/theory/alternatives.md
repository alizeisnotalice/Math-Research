# R005 alternatives and complexity audit

This is a supplementary feasibility derivation. The frozen numerical source files are `merge4_proof.md` and `registration.md`; this appendix changes none of their parameters. No numerical run is performed by this owner.

## 1. Exact complex-pair nesting is valid, but not yet the cheapest certificate

For Re z>0 define q(z)=E U^(-z), U=T1+T2, and H(w)=Gamma(w)m(w). Absolute beta/Mellin inversion gives

    B(z):=Gamma(z)q(z)
       =1/pi integral_0^infinity H(z/2+iu)H(z/2-iu)du.               (A1)

Unlike a real q, the inner integrand for complex z is not a modulus square, although it is even in u. The outer four-block identity, with s=1/4 and a=s/2=1/8, is

    q_four(s)=1/[pi Gamma(s)] integral_0^infinity
                    B(a+it)B(a-it)dt.                             (A2)

The outer real integrand is |B(a+it)|²>=0. One should evaluate B directly via(A1), avoiding the huge reciprocal Gamma(z) and its later cancellation. Gamma has no zeros; reciprocal Gamma has zeros at nonpositive integers, not extra poles. In the required Re z>0 half-plane, all factors in(A1) are holomorphic. There is no Gamma-zero obstruction.

This is two levels of scalar integration in an implementation. If one expands the squared inner integral, it is a three-real-variable integral, not a flat positive two-dimensional integrand. The inner complex phases and their propagated errors cannot be discarded. The one-dimensional Radau route avoids this nesting rather than pretending it is free.

## 2. Valid domains, majorants and tails for a possible future nested method

For an outer complex t disk with |Im t|<=rho_o<a, z=a+it has real part v in[a-rho_o,a+rho_o]. For an inner u disk with |Im u|<=rho_i, require

    rho_i<(a-rho_o)/2.                                             (A3)

Then both z/2+iu and z/2-iu have positive real part. Their real parts are(v/2-y_i,v/2+y_i). Since H(v)=Gamma(v)m(v) is log convex, a bound for the inner integrand on the whole outer/inner product disk is the maximum over the four endpoint pairs

    H((a-y_o)/2-y_i) H((a-y_o)/2+y_i),
    y_o in{−rho_o,+rho_o}, y_i in{−rho_i,+rho_i}.                    (A4)

The product is symmetric in y_i, so only two distinct values remain. This follows from convexity of the sum of the two log-H functions on the rectangle; its maximum occurs at a vertex. No unproved Gamma monotonicity is needed.

For the outer true integrand, |q(v+it)|<=q(v), and pair AM–GM gives q(v)<=2^(-v)m(v/2)^2. Hence its whole-disk bound is

    Gamma(a-rho_o)Gamma(a+rho_o)
           ·2^(-2a)m((a-rho_o)/2)^2m((a+rho_o)/2)^2.               (A5)

The outer tail uses the same recurrence as R004:

    tail_qfour(T)<=q(a)^2Gamma(a+N)^2/
                              [pi Gamma(s)(2N-1)T^(2N-1)].        (A6)

Here q(a)=q(1/8) is already certified in R004, so this tail is particularly straightforward.

The inner tail is also explicit uniformly on an outer disk. Suppose |Im z|<=W, Re z=v in[v-,v+], U0>W/2, N>=1. Gamma recurrence gives

    |H(z/2+iu)H(z/2-iu)|
      <=C(v)^2(u-W/2)^(-2N),
    C(v)=Gamma(v/2+N)m(v/2).

C is log convex, so C(v)<=max(C(v-),C(v+)). Thus the uncomputed B tail is at most

    max(C(v-),C(v+))²/
                    [pi(2N-1)(U0-W/2)^(2N-1)].                    (A7)

This is a finite bound, with no Gamma asymptotic. It may be too coarse; its actual quality must be evaluated in a future pilot. It must not be copied from the real-parameter R004 tail because the shifts Im z/2 matter.

For a tensor Taylor polynomial of the inner kernel with degrees N_o,N_i and Cauchy ratios eta_o,eta_i<1, the omitted-series supremum is bounded by

    M_inner[eta_o^(N_o+1)+eta_i^(N_i+1)]/
                             [(1-eta_o)(1-eta_i)].                 (A8)

This follows by the Cauchy coefficient bound and a union of the two omitted index tails. Integrate against the full rectangle area and sum all rectangles. Alternatively, an inner approximation uniform on an outer Cauchy circle with error e_B yields kth outer coefficient uncertainty e_B/rho_o^k; when squaring B, the product error must additionally pay2B_max e_B+e_B². An inner value error at only the real centre does not bound the outer Taylor coefficients.

As a concrete conservative cost example, take outer T=32, rho_o=3/32, h_o=1/16 (256 panels). The narrowest inner real part is1/64, so rho_i=3/256 and h_i=1/128 are valid. Inner cutoff64 would require4096 panels. A naïve tensor method has over one million rectangles, each potentially carrying two series orders129, and the existing runtime supports univariate rather than automatically nested Gamma series. Replacing outer series by many interpolation nodes also multiplies a complete inner integral for every outer node. This is a genuine implementation/error cost, not merely an instruction to increase a2D mesh.

Adaptive phase-aware contours or a better Laplace representation might improve that cost, but they need a separate finite derivation/registration. R005's1792 real1D panels use an already exercised kernel and directly certify a strict four-block correction, so they are the first acceptance gate.

## 3. A possible real-moment hierarchy after the present batch

The present Radau correction uses the first two moments of W=cosh(tau Delta)-1. More generally, choose integer d>=1, tau=a/(2d), and g(w)=T_(1/(2tau))(1+w)^(-s) when1/(2tau) is an integer. Its roots again lie to the left of the nonnegative w axis, so g is completely monotone.

A Gauss–Radau rule with a fixed node0 and d positive nodes matches moments through degree2d. Its Hermite remainder has sign g^(2d+1)(xi)w product(w-node_j)^2<=0, yielding a global polynomial majorant. Moments of W through2d use only real q(a+j tau), j from−2d to2d, and q0=1. This replaces complex q evaluations by a small univariate moment matrix.

For d=2, tau=1/32 and g=T16^(-1/4). The required real pair exponents are0,1/32,...,1/4. After R005, only four new exponents1/32,3/32,5/32,7/32 would be missing. However the W distribution and g differ from the present d=1 choice; a stronger certified numerical endpoint must be tested rather than asserted from a formal increase in moment count. No such next batch is authorised or executed by this appendix.

## 4. Joint two-block matrix determinants: what is and is not gained

Pair merging at the scalar T level still pays

    det(A+B)^(1/7)>=det(A)^(1/7)+det(B)^(1/7).

For SPD A,B let C=A^(-1/2)BA^(-1/2), rho=det(C)^(1/7). Replacing only elementary-symmetric terms1 and6 by their exact values while retaining AM–GM for the other terms gives the sharper pointwise expression

    det(A+B)>=det A [(1+rho)^7
                     +(tr C-7rho)
                     +(rho^7tr C^(-1)-7rho^6)].                    (A9)

Both corrections are nonnegative, and they vanish together only when C is scalar. This is a concrete way to quantify matrix-shape mismatch if paired diagnostics show determinant Minkowski is the main remaining slack. Its expectation involves a joint16-variable block law; it is not supplied by the scalar moments m(z) alone.

One tempting shortcut does not work: two Hadamard blocks together form a tight16-column frame. Applying the generic16-column determinant AM–GM and then four-group AM–GM reproduces the old coarse R001 Gamma bound exactly. In detail

    E det(Sigma_pair)^(-z/7)
      <=2^(-2z)[Gamma(1/2+z/16)/Gamma(1/2)]^16.

For z=1/8, the four-group bound is

    (1/2)[E det(Sigma_pair)^(-1/56)]^4
      <=(1/4)[Gamma(65/128)/Gamma(1/2)]^64=M_Gamma.

It therefore loses the improvement already obtained, and should not be run as an alleged stronger route. A useful matrix-pair method must retain mixed-minor/trace information such as(A9), or another genuinely stronger joint distribution statement.

The16-column Cauchy–Binet polynomial has binom(16,7) finite terms and positive coefficients; enumerating them under a separate finite registration could identify useful mixed-minor groups. That is a distinct theoretical/computational project, not permission to grow the current1792-panel batch or to revert to63-dimensional microscopic boxes.

The present owner has not evaluated any of these alternatives numerically. They explain why the frozen real-moment Radau computation is the preferred immediate object and identify a separate matrix-shape route if scalar block merging proves insufficient.
