# Three-cell witness family: an exact exponential-rate obstacle

This statement concerns only a specified union of grid cells certified by the upper-corner rectangles. It does not bound the complete halo of the source, and does not determine any coefficient of the general endpoint constant.

Let E_{4,k} be the subset of [-4,4]^n having at least k coordinates of absolute value at most 1. Partition the positive cube [0,4]^n into coordinate cells with intervals [0,1], [1,2], [2,4]. The union of all positive-orthant cells of counts (a,b,c), a+b+c=n, has total volume
\[
w_{a,b,c}=\frac{n!}{a!b!c!}2^c.
\]
At its upper corner use the actual rectangle with a core coordinate entirely inside [-1,1], and remaining sides [-1,x_i]. Its average is exactly
\[
P\{a+\operatorname{Bin}(b,2/3)+\operatorname{Bin}(c,2/5)\ge k\}.
\]
The two binomials are independent. Let W_{n,k} be the union of cells where this probability is strictly greater than 1/2, with all sign reflections. The upper-corner rectangle contains the entire grid cell, so every point of that cell has a rectangle with this same strictly successful average. Thus W is contained in the complete half halo. Its volume and the source volume are exact integer sums, with common sign factor 2^n:
\[
|E_{4,k}|/2^n=\sum_{j=k}^n\binom nj3^{n-j},\qquad
|W_{n,k}|/2^n=\sum_{a,b,c:\text{test}>1/2}w_{a,b,c}.
\]

We prove the exact leading rate
\[
\lim_{n\to\infty}\frac1n\max_{1\le k\le n}\log\frac{|W_{n,k}|}{|E_{4,k}|}
=\kappa_{\rm grid}
=\log\frac{t_*^{15}+t_*^{10}+2t_*^6}{t_*^{15}+3},
\]
where t_* is the unique positive root of
\[
-5t^{19}-18t^{15}+45t^9+30t^4+36=0.
\]
In particular t_*>1 and kappa_grid is strictly less than 0.414, while mu>1. This grid family therefore loses a fixed positive amount at the exponential scale; it cannot probe the desired log n correction.

Write u=(a,b,c)/n and m(u)=u_a+2u_b/3+2u_c/5. Hoeffding's inequality implies that successful cells obey k/n <= m(u)+O(n^{-1/2}), uniformly over all counts. Conversely, whenever k/n <= m(u)-n^{-1/4}, the cell succeeds for sufficiently large n. This handles the strict half threshold without assuming a median identity.

Uniform multinomial Stirling estimates give
\[
n^{-1}\log w_{a,b,c}=H(u)+u_c\log2+O(\log n/n).
\]
For alpha=k/n the source exponent is log4-D_+(alpha||1/4), where D_+=0 for alpha<=1/4 and D_+=D otherwise, uniformly up to O(log n/n). The source estimate follows by comparing the binomial tail to its largest term and using at most n+1 terms. There are only O(n^2) cell types. Uniform continuity on the compact simplex, the Hoeffding sandwich, and rational approximation of each fixed u now show that the limiting optimal rate equals
\[
\max_u\left[H(u)+u_c\log2-\log4+D(m(u)\Vert1/4)\right].
\]
Indeed m(u)>=2/5>1/4, and for fixed u the increasing function D_+(alpha||1/4) is maximized under alpha<=m(u) at alpha=m(u). For the matching lower bound first choose alpha=m(u)-delta with fixed delta>0, approximate u by types, and then let delta decrease to zero. For the upper bound use the uniform Hoeffding error and continuity of D_+, including at the endpoints alpha=0,1. Thus the order of maximization and limiting types does not hide a boundary case.

Use the elementary Bernoulli entropy identity
\[
D(m\Vert1/4)=\sup_{\theta\ge0}\{\theta m-\log(3/4+e^\theta/4)\},\qquad m\ge2/5.
\]
Both optimizations are suprema over a product domain, so they commute. The finite Gibbs identity for the three types gives
\[
\kappa_{\rm grid}=\sup_{\theta\ge0}
\log\frac{e^\theta+e^{2\theta/3}+2e^{2\theta/5}}{3+e^\theta}.
\]
Set t=e^{theta/15}. The derivative of the ratio is
\[
\frac{t^5[-5t^{19}-18t^{15}+45t^9+30t^4+36]}{(t^{15}+3)^2}.
\]
Descending nonzero coefficients have exactly one sign change, so Descartes' rule bounds the number of positive roots by one. The polynomial is positive at 1 (value 88) and negative for large t, giving exactly one root, above 1. The ratio rises to this root and then falls to 1. This proves the displayed maximum formula. Exact rational bisection and an 80-term atanh logarithm enclosure certify the reported numerical interval; Mathematica separately verifies the derivative identity and 80/120-digit root stability.

The code's floating logarithms of exact integer volume ratios are diagnostic output. Strict cell membership and all volumes are computed using integers. The rate theorem is established by the argument above, not by extrapolating those finite ratios.
