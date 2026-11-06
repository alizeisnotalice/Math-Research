# r015 independent audit — 2026-11-06

Scope: read-only audit of `rounds/r015/report.md`, its two checking scripts,
and the claimed negative-Laplace mechanism. No repository or Git changes.

## Result

Accepted as stated in its deliberately limited scope. The report correctly does
**not** turn the exponential-moment statement into a CDF/quantile/endpoint
claim, and explicitly leaves the proposed `O(log n)` quantile coefficient
unproved.

For `a=1+lambda`, direct differentiation gives
\[
 \partial_c\log M_c=\frac{t^a-at-(a-1)}{(1+t)(1+t^a)},\qquad t=e^c.
\]
The numerator has value `-2(a-1)` at one and derivative
`a(t^(a-1)-1)>0` for `t>1`, hence has one root above one. At that root,
`t^a+1=a(t+1)`, so `M_{c_lambda}=t^{-(a-1)}=exp(-lambda*c_lambda)`.
The conditional one-step bound iterates and unrolling gives (3) exactly.
Its nonnegative summands prove the a.s. uniqueness assertion.

With `p=1+1/lambda`, `A=exp(lambda*c_lambda)=t^lambda`, multiplying the
stationary equation by `p-1` gives
\[
 (p-1)A^p-pA^{p-1}-1=0.
\]
Thus `A` is the (>1) GMS root and `k_lambda=A^{-1}`; this is an algebraic
identification, not an input to proof (1).

The removable stationary equation has `F(c,0)=e^c(c-1)-1` and
`F_c(mu,0)=mu*e^mu>0`. Coefficient matching gives
\[
c_\lambda=\mu-\frac\mu2\lambda+
\frac{\mu(9-\mu)}{24}\lambda^2+O(\lambda^3),
\]
and, using `m''(mu)=-(mu-1)/mu`, `v'(mu)=mu-1`,
\[
K(c_\lambda,\lambda)-K(\mu,\lambda)
=-\frac{\mu(\mu-1)}8\lambda^3+O(\lambda^4).
\]
The sign and local `lambda(c-c_lambda)^2` coercivity are correct.

## Actual tool checks

* MATLAB R2026a independently recovered the action coefficients: at
  `lambda=10^-4`, finite-difference values were `-0.639191142231` and
  `0.411291495489`, converging to `-mu/2=-0.639232271381` and
  `mu(9-mu)/24=0.411321220740`.
* Wolfram 15 exact series returned `a1=-m/2`, `a2=m(9-m)/24`; at
  `m=1+ProductLog[e^-1]`, its cubic value agrees with
  `-m(m-1)/8=-0.0445008805420259614`.
* MATLAB Code Analyzer reported no issues. The MATLAB `station` is the
  stationary numerator divided by `lambda`, and `derivative` is exactly
  `K_c`.

## Minor wording only

If still present in §2, “排除了端点极值” is more precise as “排除了端点成为
全局极小点”: the two endpoints equal `1/(1+lambda)` and are larger than the
interior minimum. This changes no result.
