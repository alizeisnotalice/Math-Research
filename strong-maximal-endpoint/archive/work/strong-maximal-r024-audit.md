# r024 hyperbolic bounded-amplitude mechanism audit — 2026-11-06

Audit based on the supplied r024 theorem statement and the existing model
lemmas; the formal report was not yet on disk at audit time. No repository/Git
change.

## Verdict

The proposed mechanism is valid. The earlier preliminary audit incorrectly
requested an upper-tail truncation; the formal r024 text and r002 §1 show it
is not needed.

For `r=exp(-mu)`, negative tilting decreases both the exponential variance
and the Bernoulli variance (its success parameter starts below one half).
Thus
\[
 \log E e^{-\lambda(Z-\mu)}\le \mu\lambda^2/2,\qquad\lambda\ge0.
\]
At
\[
q=n\mu-\sqrt{2\mu n\log(2A/\epsilon)},\quad A=1+\epsilon,
\]
Chernoff gives `P(S<q)<=epsilon/(2A)`, hence
`A P(S>=q)>=1+epsilon/2>1`.  The step terminal has exact modulus
`J=e^{-q} Phi_n(A)` and score at least `e^q/Phi_n(A)`.

For each requested scale `epsilon=1/n`, `1/sqrt(n log n)`, and
`sqrt(log n/n)`, both the deviation loss and
`log Phi_n(1+epsilon)=log(1+epsilon)+(n-1)log(1+log(1+epsilon))`
are `o(n)`.  This gives the restricted bounded-amplitude exponent at least
`B=e^mu`; the already proved unrestricted upper exponent gives equality.
In contrast, the constant full-cube/constant-rectangle family has score
`R(A)` and its capped near-one exponential rate is one.

Strict physical realization needs no upper cutoff of the scalar terminal.
The nondecreasing bounded step `h=A 1_[q,infinity)` is exactly allowed by
r002 §1: it has a finite lower zero region and a bounded final layer extending
to infinity. For each fixed n,
`f_L(x)=h(S(x)-L)1_[-1,1]^n` is already a bounded, compactly supported actual
input of amplitude at most A. Dominated convergence in r002 applies because
h and `Phi_n(h)` are bounded; the strict `epsilon/2` margin is retained for
sufficiently large physical depth. Thus a nonmonotone upper truncation is
neither required nor desirable.
