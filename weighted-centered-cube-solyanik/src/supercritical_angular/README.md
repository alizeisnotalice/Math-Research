# Round 31 supercritical-angular diagnostics

These programs check only the explicit scalar origin-hole formula.  They do
not numerically prove the corresponding maximal-operator statement for all
inputs.  For \(\beta=1-\gamma\), \(\tau=n\beta+d\), \(e=n-\tau\),
\(P=d/\gamma\), and the deliberately named tail exponent
\(r_0=d/e\), they check
\[
u(a)=a^{n\gamma/e}+\frac{n-1}{P}a(1-a^{r_0}),\qquad
F(a)=a^{\tau\gamma/e}+\frac{\tau(n-1)}{P(n-\gamma)}a^\beta(1-a^{r_0}),
\]
against \(F/u^\beta\to A\), where
\[
A=\frac{\tau}{n-\gamma}\left(\frac{(n-1)\gamma}{d}\right)^\gamma.
\]
`tail_power` is not the paper's radial exponent \(q\).

Run from this directory:

```sh
python3 supercritical_batch.py       # 120-digit fresh run or signature replay
python3 fraction_certificate.py      # exact rational certificates
wolframscript -file supercritical_symbolic.wl
```

Then run `supercritical_independent.m` in MATLAB.  The Python and MATLAB
solvers invert in \(L=-\log a\), using logarithmic sums; the earlier 700-step
direct bisection in \(a\) could not resolve \(2^{-16384}\) and is discarded.
Python runs 224 base rows plus 16 small-\(d\) refinement rows with a
synchronized first wave of 15 workers.  Cache and full JSON output are
reproducible local products and intentionally excluded from archival.

`fraction_certificate.py` proves both \(B>7/6\) through the stronger exact
comparison \(B^{10}>(7/6)^{10}\), and the finite strict-halo parameter
identity for \(a=s^{18}\), \(s=1/10\).  It never replaces a fractional power
by a floating approximation.
