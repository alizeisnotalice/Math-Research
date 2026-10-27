# R043 · Exact fractional moments and a scoped truncation barrier

Status: new proof submitted for independent review. R042's accepted proof is unchanged. The result here is a precise limitation of mean-normalized, independent-coordinate column products and uniform column cutoffs. It is **not** a lower bound for Cn and does not rule out a new kernel, actual-μ geometry, or row-dependent payment. No production scan, n=7 work, MATLAB, or Git operation.

## 1. The exact geometric one-coordinate law

Fix r,ρ,g>0 with

\[
g\le\min(r/2,\rho/2).
\]

Let X be Uniform[0,ρ] and define the normalized target side length

\[
\ell(X)=1-\frac{(X-(\rho-g))_+}{2r}.
\]

This is exactly the free-coordinate source-face law in R042, with g=r_-−q+ζρ. Its values lie in [1−g/(2r),1]⊂[3/4,1]. For every real p>1,

\[
\boxed{\quad
\mathbb E\ell^p
=1-\frac g\rho+
\frac{2r}{\rho(p+1)}
\left[1-\left(1-\frac g{2r}\right)^{p+1}\right].\quad} \tag{1}
\]

**Proof.** The atom at ℓ=1 has probability 1−g/ρ. The excess U=(X−(ρ−g))_+ has density 1/ρ on (0,g). Integrate (1−u/(2r))^p there. ∎

In particular, writing m=Eℓ,

\[
m=1-\frac{g^2}{4r\rho},\qquad
\mathbb E\ell^2
=1-\frac{g^2}{2r\rho}+\frac{g^3}{12r^2\rho},
\]

\[
\operatorname{Var}\ell
=\frac{g^3}{12r^2\rho}-\frac{g^4}{16r^2\rho^2}
\ge\frac{g^3}{24r^2\rho}. \tag{2}
\]

The variance bound uses g/ρ≤1/2. It is a strict positive fluctuation term; p→1 does not remove it from a normalized product law.

For each fixed p>1, Taylor expansion of (1), with g/r and g/ρ tending to zero, gives

\[
\log\frac{\mathbb E\ell^p}{m^p}
=\frac{p(p-1)g^3}{24r^2\rho}
+O_p\left(\frac{g^4}{r^3\rho}+\frac{g^4}{r^2\rho^2}\right). \tag{3}
\]

Equation (3) is a fixed-p expansion. The uniform p↓1 argument below does not divide an unspecified O_p remainder by p−1.

## 2. A uniform entropy lower bound for every p>1

Let U=ℓ/m, so EU=1. Since ℓ∈[3/4,1] and 3/4≤m≤1, U∈[3/4,4/3]. The function f(u)=u logu has f''(u)=1/u≥3/4 on this interval. Taylor's lower inequality about u=1 yields

\[
\mathbb E[U\log U]
\ge\frac38\operatorname{Var}U
\ge\frac38\operatorname{Var}\ell
\ge\frac{g^3}{64r^2\rho}. \tag{4}
\]

An exact entropy formula is available. Set a=1−g/(2r); then

\[
\mathbb E[U\log U]
=\frac{r(a^2-1)/(2\rho)-(r/\rho)a^2\log a}{m}
-\log m. \tag{5}
\]

It follows by integrating y logy after the change y=1−u/(2r). Its small-gap leading term is g³/(24r²ρ), consistent with (3).

Now take N=n−1 independent free coordinates and their product Z=∏iℓi. A fixed maximal-face coordinate and any deterministic prefactor in the column envelope multiply Z by a positive deterministic number, and cancel after division by its expectation. Therefore the normalized column product is

\[
\widehat Z=Z/\mathbb EZ=\prod_{i=1}^N U_i.
\]

The convex function ψ(p)=logEU^p has ψ(1)=0 and ψ'(1)=EU logU. Thus for every p>1,

\[
\boxed{\quad
\frac1{p-1}\log\mathbb E\widehat Z^p
=\frac{N}{p-1}\log\mathbb EU^p
\ge\frac{Ng^3}{64r^2\rho}.\quad} \tag{6}
\]

This bound covers arbitrarily p close to 1 and all larger p, without a compact-p restriction. In particular the mean-normalized moment cost does not vanish with p−1. The fixed-face law has N=n−1 free coordinates, not n; no independent law is assigned to the pinned face coordinate.

## 3. What a fractional-moment Markov certificate actually costs

Let Y be a positive deterministic multiple of Z, M=EY, and suppose a proof pays its removed column energy only by

\[
\mathbb E[Y\mathbf1_{Y>HM}]
\le\frac{\mathbb EY^p}{(HM)^{p-1}}
=M\frac{\mathbb E\widehat Z^p}{H^{p-1}}. \tag{7}
\]

To make this particular upper bound at most θM for a fixed 0<θ<1, it is necessary that

\[
H^{p-1}\ge\mathbb E\widehat Z^p/\theta.
\]

Using (6), every such certificate must pay

\[
\boxed{\quad
\log H\ge\frac{Ng^3}{64r^2\rho}
+\frac{\log(1/\theta)}{p-1}.\quad} \tag{8}
\]

The cubic first term persists as p↓1. The second term diverges when θ is fixed and p↓1. For fixed p the exact expansion (3) refines the first term to approximately pNg³/(24r²ρ); it can affect constants but not the cubic scale.

### A certificate artifact must be distinguished from geometry

If Y is deterministic, its true removed energy is zero as soon as H>1. Nevertheless (7) still requires H≥θ^(−1/(p−1)), which diverges as p↓1. Thus the second term in (8) is a cost of that Markov certificate, not an intrinsic lower bound for every tail estimate. By contrast the first term comes from a strictly positive entropy of the actual geometric product law.

## 4. The cubic cost in the frozen product-envelope energy law

The following diagnostic does not use Markov's p-moment bound, but its object is the frozen-slice product envelope Y itself. In R042 only C_A≤Y is proved for the actual kernel column C_A; equality or a positive reverse bound is not established. Thus a retention necessity for Y-energy cannot be transferred to C_A-energy. Define a probability Q by dQ=widehatZ dP. Under Q the N free coordinates remain independent, each tilted by density U. Its logarithmic mean is

\[
\mathbb E_Q\log\widehat Z
=N\mathbb E[U\log U]\ge A/64,
\qquad A=\frac{Ng^3}{r^2\rho}. \tag{9}
\]

The logarithmic variance is at most A. Indeed |logℓ|≤(4/3)(1−ℓ), ℓ/m≤4/3, and

\[
\operatorname{Var}_Q(\log\ell)
\le\mathbb E_Q(\log\ell)^2
\le\frac{64}{27}\mathbb E(1-\ell)^2
=\frac{16g^3}{81r^2\rho}
\le\frac{g^3}{r^2\rho}.
\]

Independence adds these variances. Chebyshev then gives

\[
\boxed{\quad
Q\{\widehat Z\le e^{A/128}\}
\le\frac{16384}{A}.\quad} \tag{10}
\]

The fraction of Y-envelope energy retained by Y≤HM equals

\[
\frac{\mathbb E[Y\mathbf1_{Y\le HM}]}{\mathbb EY}
=Q\{\widehat Z\le H\}.
\]

For instance, if A≥65536 and H≤exp(A/128), this fraction is at most 1/4. To retain at least half the conditional Y-envelope energy, a uniform cap must therefore exceed exp(A/128). This is an exact fixed-slice envelope statement, not merely an upper-bound artifact of (7).

It is not a lower bound for Cn, nor a retention necessity for the actual C_A, a cap-conditioned or mixed-radius column law, or actual μ. The normalization is the exact conditional mean of the specified Y law at that slice. Any such transfer would require an independently proved reverse energy comparison or a newly established law. None is supplied here. An argument retaining different fractions at different radii, changing row weights, or exploiting actual μ can lie outside this statement.

## 5. A rigorous order barrier for the specified uniform-cut proof template

Suppose an attempted extension keeps the independent-coordinate column law above with g comparable to d, r comparable to 1 and ρ comparable to n/L, with comparison constants independent of n. Then

\[
A=Ng^3/(r^2\rho)\ge a_0Ld^3
\]

for a fixed a0>0. Suppose further that its certified energy-deficit lower estimate has the same uniform-cap structure

\[
\text{reported deficit}\le Cnd/H,
\]

with fixed C, and that its tail payment uses the specific product-moment upper certificate (7). Equation (8) then forces a bound of the form below. A second, separate interpretation is possible if the retained-energy target is explicitly the Y-envelope energy itself; then (10) supplies the analogous large-A envelope restriction. Neither interpretation proves a necessity for retaining actual C_A-energy merely from C_A≤Y:

\[
\text{reported deficit}\le C'nd\exp(-a_1Ld^3), \tag{11}
\]

for a fixed a1>0, apart from the bounded-A regime where d=O(L^(-1/3)) already. Calculus gives the exact scale bound

\[
\sup_{d>0}d e^{-a_1Ld^3}
=(3ea_1L)^{-1/3}. \tag{12}
\]

Thus this specified certificate structure cannot improve the order nL^(-1/3) merely by choosing p>1 or sending p↓1. If d=L^-α with α<1/3, its cubic column cost grows like L^(1−3α); a uniform cap then loses exponentially more than the polynomial gain in nd.

The phrase “reported deficit” is essential: this upper bound limits what the displayed moment-payment template certifies, or the separately specified Y-energy normalization, not the actual missing mass Δ and not Cn. For the bounded-A envelope case, retaining a fixed fraction η>0 forces H≥η, since Q{widehatZ≤H}=E[widehatZ 1{widehatZ≤H}]≤H; this closes the bounded-A contribution to the stated order. No conclusion is drawn for arbitrary geometric kernels, actual C_A without a reverse comparison, variable cap thresholds, nonconstant retained-energy fractions, joint witness-mass payment, or additional actual-history constraints.

## 6. Concrete next routes left open

The exact column law suggests that a genuine order improvement must alter a cost in the template, rather than merely substitute a fractional moment into R042. Examples are:

- Use row-dependent weighting or a new transport dual that balances high columns without dividing the whole energy by one common H.
- Use actual-μ geometry or the common-witness mass dual, rather than λ-product moments alone, to pay high columns differently.
- Redesign the radial/capture kernel so that the tilted logarithmic variance and entropy cease to have size Ld³.

These are not proved positive results here. The proved result is the exact formula (1), the uniform all-p moment-certificate trade (8), and the scoped fixed-slice Y-envelope diagnostic (10).

## 7. Small symbolic checks

Default Python has no SymPy in this environment, so `check_moment_identities.py` uses standard-library exact Fraction polynomials instead. It checks the closed moment formula against direct polynomial integration for p=1,...,8, checks two fractional exponents by a rational-power substitution, and checks the variance identity and conservative constants. The finite checks diagnose algebra only; the all-real-p and all-n conclusions follow from the analytic proof above. No production scan is run.

