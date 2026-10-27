# Pair-Mellin certificate: frozen analytical specification

2026-10-27. Sol derivation only, no numerical run by this owner. This file is for Terra's separate `upper/round004/pair_mellin/` computation. The original equal-weight body64/common-time Cauchy model is fixed. No g or C7 claim follows from this computation alone.

## 1. Exact moment and expectation bound

Partition64 directions into8 disjoint7x8 tight Hadamard blocks, as proved and checked in R002. For a block Sigma_g=sum_(i in g)S_i v_iv_i^T, S_i=1/(4Y_i), Y_i iidGamma(1/2,1), let T_g=det(Sigma_g)^(1/7). The8 T_g are iid positive nondegenerate variables. For Re z>0,

    m(z)=E T^(-z)
      =2^(-4z/7) Gamma(4+z)/Gamma(4+8z/7)
                      [Gamma(1/2+z/7)/Gamma(1/2)]^8.                (1)

Indeed the block determinant is16(sumY)/(productY). Write Y_i=R p_i, R~Gamma4 independent of p~Dirichlet(1/2)^8. Then detSigma_g=16R^(-7)/(product p), and direct Gamma/beta integration proves(1), including complex z by absolute integration for Re z>0. In particular |m(v+it)|<=m(v), v>0; m is analytic in that half-plane and log m(v) is convex by Holder.

SPD determinant Minkowski, whose self-contained R002 proof applies, gives

    M=E det(Sigma)^(-1/14)<=E(sum_(g=1)^8 T_g)^(-1/2).

For U_j=T_(2j-1)+T_(2j), four independent iid pair sums, AM–GM gives

    M<=M_pair:=(1/2)[q_s]^4, s=1/8,
    q_s=E(T1+T2)^(-s).                                             (2)

The factor is4^(-1/2)=1/2 and each pair exponent is1/8, not1/16.

## 2. Positive one-dimensional integral, including its constant

Put a=s/2=1/16. The beta integral under y=e^r gives

    Gamma(a+it)Gamma(a-it)/Gamma(2a)
      =integral_R e^(itr)[2cosh(r/2)]^(-2a)dr.

Fourier inversion and

    (x+y)^(-s)=(xy)^(-s/2)
                    [2cosh((log x-log y)/2)]^(-s)

therefore give, after independent expectation,

    q_s=1/[2pi Gamma(s)] integral_R
             Gamma(a+it)Gamma(a-it)m(a+it)m(a-it)dt
       =1/[pi Gamma(s)] integral_0^infinity f(t)dt,                 (3)
    f(t)=|Gamma(a+it)|² |m(a+it)|²>=0 on real t.

Fubini is absolute: E(T1T2)^(-a)=m(a)^2 and the Gamma-product integral is finite by recurrence below. The factor1/(pi Gamma(s)) in the half-line formula is thus exact.

Pair AM–GM implies q_s<2^(-s)m(s/2)^2, strict since T is nondegenerate continuous. Hence M_pair<8^(-1/2)m(1/16)^8=M_H analytically. This does not supply a numerical size of the improvement; q_s and the resulting budget still require the finite certificate.

## 3. Holomorphic integrand and full Cauchy-envelope proof

Use the analytic continuation in the **integration variable** w,

    f(w)=Gamma(a+iw)Gamma(a-iw)m(a+iw)m(a-iw).                       (4)

It is holomorphic for |Im w|<a. If w=t+iy, the real parts of a+iw and a-iw are a-y and a+y. The Gamma integral implies |Gamma(v+iu)|<=Gamma(v), v>0; the moment integral implies |m(v+iu)|<=m(v). Thus

    |f(t+iy)|<=h(a-y)h(a+y), h(v)=Gamma(v)m(v).

Both Gamma and m are positive log-convex functions on v>0, so h is log-convex. For |y|<=r<a, write theta=(y+r)/(2r) in[0,1]. Then

    h(a+y)<=h(a+r)^theta h(a-r)^(1-theta),
    h(a-y)<=h(a+r)^(1-theta)h(a-r)^theta.

Multiplying proves the exact endpoint envelope

    |f(t+iy)|<=M_*:=h(a-r)h(a+r)                                   (5)

for every real t and |y|<=r. No unproved monotonicity of Gamma is used. Choose r=3/64, so a-r=1/64>0 and a+r=7/64. Every relevant Cauchy circle remains strictly within the analyticity strip.

## 4. Fixed finite segment and Taylor remainder constant

Take T=32. Split[0,32] into512 panels

    panel_j=[j/16,(j+1)/16], centre t_j=(2j+1)/32,
    halfwidth h=1/32, j=0,...,511.

Compute the Taylor coefficients b_(j,k) of(4) about the **exact rational** centre t_j through degree128, with outward512-bit complex ball power series and series order at least129. The Taylor tail for |u|<=h is bounded, by Cauchy on radius r, by

    M_* sum_(k=129)^infinity(h/r)^k
      =M_* (2/3)^129/[1-2/3].

Each panel has length2h=1/16. Summing512 integrated tails yields

    E_raw_quad=32 M_* (2/3)^129/[1-2/3].                            (6)

The factor32 is the **total interval length**, not512 or64. Retain degree128, hence first omitted index129.

Integrate every polynomial term, retaining possible negative coefficients and all interval radii:

    P_j=sum_(k=0)^64 [2h^(2k+1)/(2k+1)] b_(j,2k).                  (7)

The exact coefficients are real since f is real on the real axis. Certified imaginary intervals should contain0; otherwise stop as an implementation inconsistency. Do not drop a negative term, take absolute coefficients, or treat a node sum as a whole-panel integral. The true finite integral lies in sum_j Re(P_j)+[-E_raw_quad,E_raw_quad].

## 5. Full positive tail without asymptotic Gamma estimates

For integer N>=1, recurrence and the defining Gamma integral imply

    |Gamma(a+it)|=|Gamma(a+N+it)|/
                         product_(j=0)^(N-1)|a+j+it|
       <=Gamma(a+N)|t|^(-N).

Therefore f(t)<=m(a)^2Gamma(a+N)^2t^(-2N). Choose N=16. The omitted tail in q_s is bounded by

    E_q_tail=m(a)^2Gamma(a+16)^2/
                      [pi Gamma(s)·31·32^31].                     (8)

This pays all t>=32. There is no numerical Gamma asymptotic or unrecorded tail cutoff. With the exact scalar factor c_q=1/[pi Gamma(s)], assemble

    q_lower=max(0,c_q·(sum Re P_j -E_raw_quad)),
    q_upper=c_q·(sum Re P_j +E_raw_quad)+E_q_tail.                   (9)

All interval operations use outward endpoints; the formulas denote the appropriate lower/upper endpoints of the entire expressions. The positive tail need not be subtracted from the lower bound.

## 6. q to M to complete I, and comparison scope

The original full-boundary Gaussian/Minkowski identity, proved in R001 and reconstructed in R002, is

    I=integral_[0,1]^6 G(1,y)^(7/6)dy
      <=M^(7/6)/(105·2^(25/6)pi^(13/12)).

Thus a certified q_upper gives

    M<=M_pair_upper=(1/2)q_upper^4,
    I<=I_pair_upper=((1/2)q_upper^4)^(7/6)/
                             [105·2^(25/6)pi^(13/12)].             (10)

This is the **complete** original-G/time/full-boundary budget. No short-time correction is added and no per-shape positive-face symmetry is assumed. Compare its outward upper endpoint against a freshly enclosed lower endpoint of the exact R002 M_H/I_H **budget expressions**. Their lower endpoints are not lower bounds on true M/I. If separated, round the new I upper endpoint outward to a declared rational grid and record the strict gap.

No accurate q has been evaluated by this owner. There is no numerical improvement claim before Terra's certificate. An improved I bound alone supplies no g and no C7 upper bound.

## 7. Computation contract and implementation checks

Terra owns all runtime probing and numerical files in its separate `pair_mellin/` directory. A gamma operation on a local complex Taylor series is mathematically valid here: all Gamma arguments have positive real part on each Cauchy disk, Gamma has no zeros, and m's reciprocal Gamma denominators are safe analytically. This does not assert the installed Python/FLINT API spelling or default series cap; the actual interface must be inspected. Direct Gamma series avoids unnecessary log-Gamma branch issues. Implement2^(-4z/7) as exp[(-4log2/7)z], with positive real log2 interval.

A frozen four-panel pilot j=0,1,2,3 should check exact rational centres/order129, f(0)=Gamma(a)^2m(a)^2, coefficient realness, and local Gamma recurrence Gamma(z+1)=zGamma(z) at series level. Overlapping interval identities are implementation checks, not a second arithmetic kernel. Reuse the four accepted panel records in the512-panel batch; do not silently re-run them.

After pilot acceptance,15 processes split the remaining fixed panel indices, with BLAS/OpenMP1 and immutable code/config. No extra panel, t cutoff, degree, epsilon or MC sample is authorised by unused CPU capacity. Any inability to compute an interval coefficient or separate a comparison is an explicit finite-batch failure/uncertainty, not a license to publish the midpoint as certified. The main owner may explicitly register a repair if necessary.

## 8. Independent Laplace/Plancherel check of the coefficient

For an independent confirmation of(3), let L(u)=E exp(-uT). The elementary Laplace identity x^(-s)=Gamma(s)^(-1)integral_0^infinity u^(s-1)e^(-ux)du gives

    q_s=Gamma(s)^(-1)integral_0^infinity u^(s-1)L(u)^2du.

Set g(r)=exp(ar)L(exp r), a=s/2. It belongs to L1 because its integral is Gamma(a)m(a), and to L2 because the last q integral is finite by pair AM–GM. Under the Fourier convention ghat(t)=integral exp(-itr)g(r)dr, the substitution u=exp r and absolute Tonelli/Fubini give

    ghat(t)=Gamma(a-it)m(a-it).

Plancherel then gives

    q_s=1/[2pi Gamma(s)]integral_R
                        |Gamma(a+it)m(a+it)|²dt,

exactly(3), including its half-line coefficient. This also identifies the positive integrand as the logarithmic Laplace-transform energy of one block, rather than an arbitrary contour formula. Both derivations are purely analytical; no numerical check is asserted here.
