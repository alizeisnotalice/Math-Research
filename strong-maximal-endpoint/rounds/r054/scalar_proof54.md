# Explicit finite-dimensional strong-maximal lower bound

Let n=64, c=31/25, r=exp(-c), p=r/(1+r), and let H be the exact monotone staircase specified by scalar54.json: zero below b0>0, heights 1+j/512 between the 513 boundaries, and height 2 after the last boundary. All boundaries are exact finite-decimal rationals. Its origin as an optimized profile is irrelevant to feasibility.

Let T=G+cK with independent G~Gamma(64,1), K~Binomial(64,p). The directed verifier certifies E H(T)>1.000001 and, with L=10^8 and alpha=1-63/L,
\[
J_\alpha=\int_0^\infty e^{-\alpha t}\Phi_{64}(H(t))dt
<6.583764\,10^{-32}<e^{-71.79}.
\]
The integral is exactly a finite sum of step-height costs times exponential differences divided by alpha. Each survival probability used for the mean is the positive finite sum
\[
\Pr(T\ge t)=\frac{1}{(1+e^{-c})^{64}}\left[
\sum_{ck\ge t}{64\choose k}e^{-ck}
+e^{-t}\sum_{ck<t}{64\choose k}\sum_{j=0}^{63}\frac{(t-ck)^j}{j!}\right].
\]
The j=0 term is included. No quantile/root residual is needed: the exact supplied boundaries define the input directly.

Define, almost everywhere,
\[
S(x)=\sum_{i=1}^{64}-\log|x_i|,\qquad
f_L(x)=H(S(x)-L)1_{[-1,1]^{64}}(x).
\]
Assign any value in [0,2] on the coordinate axes. For x in the cube with S(x)>=L use the actual rectangle whose ith interval is [-rx_i,x_i] for x_i>0 and [x_i,-rx_i] otherwise. It stays in the cube. A uniform point in the interval has logarithmic increment Exp(1)+c Bernoulli(p), independently across coordinates. Monotonicity gives average E H(S(x)-L+T)>=E H(T)>1. Thus the actual strong-maximal halo contains {S>=L} in the cube up to a null set. If rectangles are defined with the point in their interior, continuity of their integrals under a sufficiently small enlargement preserves the strict margin.

Set K_L=2^64 e^{-L} L^63/63!. Logarithmic polar coordinates give
\[
|\{S\ge L\}\cap[-1,1]^{64}|=K_L\int_0^\infty e^{-t}(1+t/L)^{63}dt\ge K_L,
\]
\[
\int\Phi_{64}(f_L)=K_L\int_0^\infty e^{-t}(1+t/L)^{63}\Phi_{64}(H(t))dt\le K_LJ_\alpha.
\]
Here the integrand below b0 is zero and (1+t/L)^63<=exp(63t/L) for t>=0. Consequently this one finite, bounded, compactly supported input satisfies
\[
\boxed{\frac{|\{M_S f_L>1\}|}{\int\Phi_{64}(f_L)}\ge J_\alpha^{-1}>e^{71.79}.}
\]
No limiting exchange, omitted orthant factor, or optimality claim is involved.

## Arithmetic contract

The Python verifier uses 110-digit directed Decimal intervals for +,-,*,/; exp and ln are correctly rounded nearest by the documented Decimal contract and are enlarged by one adjacent representable number at each endpoint. Integer powers use repeated interval multiplication. Finite-decimal preprocessing uses a 500-digit Inexact-trapped context; overflow or unsupported inputs fail. Parameter guards are explicit exceptions and survive python -O. The JSON's claimed intervals are ignored on replay and recomputed from its exact parameter strings. Independent Wolfram high-precision replay is numerical cross-validation, not the source of rigor.

Source: [Python Decimal documentation](https://docs.python.org/3/library/decimal.html), rounding contract and next_minus/next_plus. Geometric transfer follows the project's r002/r009/r028 constructions, with an explicit finite-L domination here. No public priority claim is made.
