# Independent review of the r081 exterior model

## Verdict

**PASS, with the stated sector/model scope.** The exterior-volume formula,
the Orlicz lower estimate, the flat-cutoff loss, and the dyadic-strip
accounting are consistent. No new general strong-maximal upper bound follows from these model identities.

A local Python sanity check evaluated the first terms at \(N=8,32\), finding
\(V/V_0-1=0.01263,0.00359\) and
\(V'/V-2N/\log2=0.0362,0.01034\), respectively; it also verified the
displayed flat-cutoff tail inequality at dyadic exponents \(8,9,10\).
These computations are consistency checks only, not used as proof.

## 1. Exterior-coordinate formula

For \(x_i>1\), write \(u(x_i)=e^{-v_i}\). Then
\(x_i=\cosh v_i\), \(dx_i=\sinh(v_i)\,dv_i\), and for
\(F_N=2\prod_i s_i{\bf1}_{[0,1]^N}\),
\[
 M_{\rm S}F_N(x)=2e^{-\sum_i v_i}.
\]
Consequently the exterior part of the level-one halo is \(\sum v_i<T\),
\(T=\log2\), and expansion of each \(\sinh v_i\) gives exactly
\[
 V_N(T)=\sum_{k\ge0}{N+k-1\choose k}{T^{2N+2k}\over(2N+2k)!}.\tag{1}
\]
The successive-term ratio is
\[
 {a_{k+1}\over a_k}
 ={N+k\over k+1}{T^2\over(2N+2k+2)(2N+2k+1)}
 \le {T^2\over2(2N+1)}.
\]
Thus
\[
 V_N(T)={T^{2N}\over(2N)!}\bigl(1+O(N^{-1})\bigr).     \tag{2}
\]
This is an exterior-sector volume only; it is not the entire strong halo.

For the modular, \(\int F_N=2^{1-N}\). On the exceptional region
\(F_N>1\), put \(z=-\sum_i\log s_i\in[0,T]\). The excess modular is
bounded, after division by \(2^{1-N}\), by
\[
 {2^N\over(N-1)!}\int_0^T
 z^{N-1}(1+T-z)^{N-1}\,dz
 \le{(2T)^N\over(N-1)!},                                \tag{3}
\]
because \(z(1+T-z)\le T\) on \([0,T]\). Hence the stated exterior
halo/modular ratio is extremely small; it cannot imply divergence of an
endpoint constant.

## 2. Flat-cutoff exterior loss

Let \(V_{\rm soft}\) mean the smooth integral over the **same exterior
sector**, with factor \(\chi((M_{\rm S}F_N-1)/\eta)\). The loss is
\[
 V_N-V_{\rm soft}=
 \int_0^\eta(1-\chi(h/\eta))
 {f_N(T-\log(1+h))\over1+h}\,dh,                       \tag{4}
\]
where \(f_N=V_N'\). Differentiating (2), or differentiating the positive
series (1), gives at fixed \(T=\log2\)
\[
 {f_N(T)\over V_N(T)}={2N\over T}+O(N^{-1}).            \tag{5}
\]
For \(0\le h\le\eta=N^{-2}\), the ratio in (4) relative to \(f_N(T)\)
is \(1+O(N\eta)\): the first term has logarithmic derivative \(O(N)\),
and the series correction has smaller derivative. The symmetry of \(\chi\)
gives \(\int_0^1(1-\chi(r))dr=1/2\). Equations (4)--(5) therefore yield
\[
 1-{V_{\rm soft}\over V_N}={1\over N\log2}+O(N^{-2}).  \tag{6}
\]
The qualification “exterior sector” is essential: (6) does not by itself
describe the smooth norm of the complete strong-maximal halo.

## 3. Dyadic strip accounting

With \(\delta_j=(96\,2^j)^{-1}\) and
\(\varepsilon_j=\delta_{j+1}\), the mesh prescription in `models80.md`
becomes
\[
 M_j=\left\lceil\sqrt{2^{j+2}/3}\right\rceil.
\]
It covers the strip \(\delta_{j+1}\le\delta\le\delta_j\). On that strip
the previous estimate \(M_{\rm S}F-1\le6\delta\le6\delta_j\) shows that
giving this family weight \(\chi(6\delta_j/\eta)\) dominates the desired
soft density pointwise. This is a fractional smooth-density cover, not an
exact unit-weight cover.

For \(J=k+4\), \(\eta=2^{-k}\), the first tail argument is \(2^{-8}\).
For \(m\ge1\),
\[
 \chi(2^{-m})\le e^2e^{-2^m}\le8\,2^{-2^m}.            \tag{7}
\]
Together with \(M_{j+1}+1\le2(M_j+1)\), the resulting tail is bounded by
\(16(M_J+1)2^{-256}\), as asserted. This only controls the displayed
dyadic construction near this two-dimensional boundary; it is not a
global column-budget estimate.
