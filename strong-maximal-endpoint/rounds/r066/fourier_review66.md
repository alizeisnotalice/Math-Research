# Independent audit of `fourier_proof66.md` and the three Wolfram scripts

## Verdict

The unconditional, projected Fourier statement is **correct after a small
proof repair**: the first \(c\)-derivative is uniformly bounded on the
zero-sum hyperplane, and at \(c=\mu\) it is uniformly bounded without the
zero-sum restriction.  The claimed second-derivative obstruction is also
correct, with the exact complex limits recorded below.  Neither statement
transfers to the conditional simplex density/kernel.  `conditional66.wl`
does implement the appropriate *fixed tilt, fixed observed sum* score at its
mean-matched evaluation point; it must not be read as a theorem about the
conditional kernel or as a derivative of its displayed centered CF alone.

Here \(q=e^{-c}\), \(p=q/(1+q)\), and
\[
 \ell_c(t)=\partial_c\log\psi_c(t)
 =p+{q e^{ict}(-1+it)\over1+q e^{ict}}.
\]
On every compact \(I\Subset(0,\infty)\), denominators in this formula are
uniformly separated from zero.  With \(m(c)=1+cp(c)\), Taylor expansion at
zero and a direct global estimate give
\[
 \ell_c(t)=i m'(c)t+O_I(h(t)),\qquad
 |\ell_c(t)-i m'(c)t|\le C_Ih(t),\qquad
 h(t)=\min(t^2,|t|). \tag{A}
\]
The global form of (A) follows by using the Taylor estimate for

\(|t|\le1\), and the elementary \(O_I(1+|t|)=O_I(|t|)\) bound otherwise.

For \(\sum_i t_i=0\), the linear terms in (A) cancel, so, with
\(w(t)=(1+t^2)^{-1/2}\),
\[
 |\partial_c P_c(t)|\le C_I\Big(\prod_iw(t_i)\Big)\sum_i h(t_i). \tag{B}
\]
This is dimension-free.  For completeness, split the sum into
\(|t_i|\le1\) and \(|t_i|>1\).  If \(S=\sum_{|t_i|\le1}t_i^2\), then
\[
 \prod_{|t_i|\le1}w(t_i)
 \le e^{-S/4},\qquad S e^{-S/4}\le4/e.
\]
If there are \(k\) large coordinates, their contribution is at most
\[
 \sum_{|t_i|>1}|t_i|\prod_jw(t_j)
 \le k\,2^{-(k-1)/2},
\]
which is uniformly bounded.  This supplies the missing elementary product
damping argument.  The display in `fourier_proof66.md` containing a comma after the
sum is not a valid inequality and should be replaced by this split proof.

The defining equation for \(\mu\) gives
\[
 p(\mu)=1-\mu^{-1},\quad m(\mu)=\mu,\quad
 m'(\mu)=p(\mu)-\mu p(\mu)(1-p(\mu))=0.
\]
Thus (A)--(B) apply at \(c=\mu\) without imposing \(\sum_i t_i=0\).

## Exact second-derivative counterexample

Put \(u=-1/M\),
\[
 S_M(c)=\psi_c(u)^{M^2},\qquad G_M(c)=\psi_c(M),\qquad P_M=S_MG_M.
\]
At \(c=\mu\), the CF expansion and the variance
\(v(\mu)=1+\mu^2p(\mu)(1-p(\mu))=\mu\) give
\[
 e^{i\mu M}S_M(\mu)\longrightarrow e^{-\mu/2}. \tag{C}
\]
Moreover \(S'_M(\mu)=O(S_M(\mu))\) and
\(S''_M(\mu)=O(MS_M(\mu))\).  The latter is harmless because
\(G_M(\mu)=O(M^{-1})\).  Directly differentiating the one large factor
gives the particularly clean limits
\[
 e^{-i\mu M}G'_M(\mu)\to-p(\mu),\qquad
 {e^{-i\mu M}\over M}G''_M(\mu)\to-i p(\mu). \tag{D}
\]
For example, the leading term in \(G'_M\) is exactly
\(q e^{icM}(-1+iM)/[(1+q)(1-iM)]\), and
\((-1+iM)/(1-iM)=-1\).  Applying the product rule to (C)--(D) proves the
exact assertions
\[
 P'_M(\mu)\to-p(\mu)e^{-\mu/2},\qquad
 {P''_M(\mu)\over M}\to-i p(\mu)e^{-\mu/2}. \tag{E}
\]
Hence a dimension-free second-derivative bound fails even on the
zero-sum hyperplane.  Local `wolframscript` evaluation independently gave

\[
 \begin{array}{c|cc}
 M&P'_M(\mu)&P''_M(\mu)/M\\ \hline
 64&-0.11487675-0.00001533i&0.00039099-0.11473642i\\
 256&-0.11487177-0.00002205i&0.00014819-0.11475143i
 \end{array}
\]
where \(p(\mu)e^{-\mu/2}=0.114938670136316\ldots\).  This agrees with
`second66.json` and its complex residuals, rather than only with magnitudes.

## What the valid first-derivative result does imply

Let \(H_b=\{x\in\mathbb R^b:\sum_i x_i=0\}\), and let \(T_c\) be
convolution on \(H_b\) by the projection of the iid \(Z_c\) vector.  Its
multiplier is \(P_c\) restricted to \(H_b\).  Plancherel and (B) give,
uniformly in \(b\) and \(c,d\in I\),
\[
 \|T_c-T_d\|_{2\to2}
 \le\sup_{\xi\in H_b}|P_c(\xi)-P_d(\xi)|
 \le C_I|c-d|. \tag{F}
\]
For a strongly differentiable path, Bochner FToC followed by Minkowski gives
the genuine path-maximal consequence
\[
 \left\|\sup_{c\in I}|T_cf-T_{c_0}f|\right\|_2
 \le\int_I\|\partial_cT_cf\|_2\,dc
 \le C_I|I|\,\|f\|_2. \tag{G}
\]
This is a useful unconditional *projected* \(L^2\) result.  It is not a
bound for \(g_c^s\), the conditional simplex kernel.  In particular, it
does not yield the required rare-deficit average \(o(L/N)\): projection

loses the conditioning normalizer, and passing to a conditional density
requires division by \(f_{b,c}(s)\).  Although \(u\mapsto u_+\) is

Lipschitz in \(L^2\), (G) provides neither a conditional/TV estimate nor the
needed rare-event rate after that normalization.

## Conditional Fourier algorithm

There is no algebraic score error in `conditional66.wl` under its stated
pointwise interpretation.  Hold \(a=1-\theta\) and the observed \(s\) fixed,
and define the **uncentered** tilted log-CF
\[
 L_c(x)=-\log(1-ix/a)+\log(1+e^{-ca+icx})-\log(1+e^{-ca}).
\]
Its fixed-(a) score is exactly the script's `sc`:
\[
 \partial_cL_c(x)=ap+
 {e^{-ca+icx}(-a+ix)\over1+e^{-ca+icx}}. \tag{H}
\]
The script displays/integrates the centered value
\(\lambda_c(x)=L_c(x)-im_cx\), where
\(m_c=a^{-1}+cp\).  At the selected evaluation point
\(s=bm_{c_0}\) and for the two-value frequencies with sum zero,
\[
 e^{-i\tau s}\exp\!\left\{\sum_iL_{c_0}(t_i+\tau)\right\}
 =\exp\!\left\{\sum_i\lambda_{c_0}(t_i+\tau)\right\}.
\]
When differentiating the original conditional inversion at fixed \(s,a\),
the score is \(\sum_i\partial_cL_c(t_i+\tau)\), precisely what the code
uses.  It would be incorrect to add \(-im'_c x\) to `sc` in this use: that
term differentiates the *centered representation* and is cancelled by the
derivative of the fixed-(s) Fourier phase.  Conversely, `sc` should not be
described as \(\partial_c\lambda_c\), and the calculation is not valid away
from the mean-matched point without restoring the phase.

The code remains a finite-pattern, finite-cutoff trapezoidal diagnostic, as
its JSON scope correctly says.  An imaginary normalizer residual alone is
not a quadrature/truncation error certificate; the coarse/refined comparison
is useful numerical evidence but cannot establish a uniform conditional
operator or maximal theorem.  `second66.wl` has no apparent algebraic error in its stated numerical scope. Earlier exploratory frequency scans remain outside the archive.
