# Round 30 critical-angular scalar diagnostics

This directory is a minimal, reproducible numerical check of the scalar
envelope used in Round 30.  It is deliberately not a numerical proof of the
maximal-operator theorem or its all-input upper bound.

For \(n=2,3\), \(\gamma\in\{1/4,1/2,3/4,9/10\}\), and
\(u=2^{-k}\), \(k\in\{8,16,32,64,128,256,512,1024,4096\}\), the programs
check the unique nonnegative solution of

\[
 L-\log(1+cL)=\log(1/u),\qquad c=(n-1)/n,
\]

and evaluate \(F(u)=e^{-\beta L}(1+bL)\), with
\(\beta=1-\gamma\) and \(b=(n-1)\beta/(n-\gamma)\).  Values are normalized
by \(u^\beta\log(1/u)^\gamma\) and compared to the normalized leading
coefficient \(J=bc^{-\beta}\).  This is distinct from the raw cube mass
\(\kappa=2^n(n-\gamma)/(n\beta^2)\) in the slice formula.  The diagnostics
are the local envelope, a long-class core term, and a strict-hole mass:

\[
C=\log(1/u)^{1/(2n)},\quad u_+=uC^n,
\qquad \delta=1/\log(1/u),\quad u_-=u(1-2\delta)^n/\delta^n.
\]

The `upper_normalized` value is only \(C^{-\tau}F(uC^n)\) divided by the
normalizer: it omits the long-class denominator and short-class term, so it
is not a complete upper bound.  The lower expression is the current
strict-hole version; the Python table
also emits the old \((1-\delta)^n\) expression only as a non-evidentiary
asymptotic comparison.  A row is marked geometry-valid only when
\(C>2\), \(\delta<1/3\), and \(u<1\).  The proof may use the optimized
choice \(C=\log(1/u)^{\gamma/(1+n\gamma)}\); this fixed-choice diagnostic
does not attempt to certify that rate optimization.

Run from this directory:

```sh
python3 critical_angular_batch.py       # fresh or signature-guarded cache replay
wolframscript -file critical_angular_check.wl
```

Then run `critical_angular_independent.m` in MATLAB.  Python uses 120-digit
`Decimal` bisection and launches a synchronized first wave of 15 workers;
Wolfram uses 240-digit numeric `ProductLog[-1]`; MATLAB independently uses
log-domain double bisection, so no `exp(-H)` underflow occurs for \(k=4096\).
The generated `runtime_cache/` and full `round30_summary.json` are ignored
from archival on purpose; only the compact result summary is versioned.
