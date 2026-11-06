# High-dimensional Orlicz price of the fixed-arity obstruction family

This is an exact cost/point-regret theorem for the explicit r057 family. It does not assert that all high-arity inputs obey this tradeoff, or that this family is near-extremal for C_N.

## Actual input and exact cost

Fix integers N>=2, M>=2 and 0<delta<=1/2. On [0,M], let
\[
g_i(x)=1-\delta 1_{[i-1,i]}(x)+\delta 1_{\{i<M\}}1_{[i,i+1]}(x),\quad i=1,\ldots,M,
\]
and set g_i=0 outside [0,M]. On outer slab i take f(x,y)=g_i(x), y in [i-1,i]. Embed this in N dimensions by multiplying by 1_[0,1] in each remaining coordinate. Then f is nonnegative, compactly supported and bounded by 3/2; its support volume is M². The modular is the actual Phi_N(t)=t(1+log^+t)^(N-1) modular, not a sum of Phi_2 fiber budgets.

At inner receiving coordinate zero and fixed outer interval J=[0,M], every fiber maximum is one. The average of fiber i on [0,j] is one except at i=j, where it is 1−delta/j. Every proper subset has a common maximizing prefix. The common average of all fibers is maximized at [0,M] and equals 1−delta/M². Thus
\[
R_J=\delta/M^2>0,\qquad \Delta_k=0\quad(k<M). \tag{1}
\]
Intervals extending to the left of zero only add zero mass and cannot improve any average. Between integer endpoints each prefix average is fractional linear, so the endpoint computation controls all intervals, not just a grid.

There are M downshifted cells of height 1−delta and M−1 upshifted cells of height 1+delta. Consequently the exact N-dimensional modular and its normalized cost are
\[
J_N=M^2-\delta+(M-1)(1+\delta)
 \big[(1+\log(1+\delta))^{N-1}-1\big],
\]
\[
\mathcal B_N:=J_N/M^2
=1-R_J+{M-1\over M^2}D_N(\delta),\quad
D_N(\delta)=\Phi_N(1+\delta)-(1+\delta). \tag{2}
\]
This specializes to the r057 Phi_2 cost, but keeps the essential N-dependent exponential price.

## Budget-regret inequality

For 0<=s<=1/2, log(1+s)<=s, and therefore
\[
{d\over ds}\log(1+\log(1+s))
={1\over(1+s)(1+\log(1+s))}
\ge{1\over(1+s)^2}\ge4/9.
\]
Integration from zero gives
\[
D_N(\delta)\ge\exp\{\tfrac49(N-1)\delta\}-1. \tag{3}
\]
Since 1−R_J>0, equation (2) implies
D_N<=M² B_N/(M−1)<=2M B_N. Combining with (1)--(3),
\[
\boxed{R_J\le {9\over4(N-1)M^2}
 \log(1+2M\mathcal B_N).} \tag{4}
\]
This is the quantitative cost of this explicit high-arity obstruction.

If B_N is bounded by a constant independent of N, then R_J=O(1/N) uniformly in M>=2: log(1+2M B_N)/M² is bounded. More generally B_N<=N^A and M<=N^C give R_J=O_{A,C}(logN/(N M²)). The polynomial-M restriction is automatic if one also requests R_J>=kappa(logN)²/N and delta<=1/2, because then M²<=N/[2kappa(logN)²]. These facts concern this family only.

In particular, suppose for some fixed kappa>0 that
\[
R_J\ge\kappa{(\log N)^2\over N}.
\]
For N>=2, (2)--(3) give directly
\[
\mathcal B_N\ge {1\over2M}
 \left[\exp\{\tfrac{2\kappa}{9}M^2(\log N)^2\}-1\right]. \tag{5}
\]
Once the exponent exceeds log2,
\[
\boxed{\log\mathcal B_N\ge
\tfrac{2\kappa}{9}M^2(\log N)^2-\log(4M).} \tag{6}
\]
The required log-cost is order at least (logN)², larger than the unresolved logN endpoint coefficient. Thus the r057 family cannot become a feedback-scale point-regret obstruction while keeping a polynomial normalized Orlicz budget. This does not exclude other high-rank constructions with a better budget, and it does not turn a fixed-outer point calculation into an endpoint theorem.

## Local computation and next gap

Native Wolfram 15 exactly integrated the actual cell matrices for N=2,4,8,16; M=2,3,5; delta=1/2,1/3, and checked (2) plus the elementary derivative inequalities. There are 29 exact checks. The same script evaluates 15 exact-symbolic parameter choices at 40- and 80-digit numerical precision for N=64,...,16384 and M=2,4,8. The two precisions agree within 1e-30 in the logged cost on every reported row. These are numerical diagnostics, not directed certificates; rows hitting delta=1/2 do not reach the selected critical regret target.

For example M=2 and delta=4(logN)²/N gives the critical regret when delta<=1/2. The log normalized costs are approximately 83.4089 at N256, 161.1471 at N1024 and 258.0185 at N4096. High cost is an output of the actual Phi_N formula, not a probability-tail or separate-fiber surrogate.

To improve the general upper bound, one still needs a theorem for arbitrary dangerous inputs: either controlled incompatibility rank follows from the global budget, or another family realizes large rank at a much lower cost. Equations (4)--(6) settle only this particular counterexample route.
