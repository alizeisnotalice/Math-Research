# Round 32 flat concentration diagnostics

The programs check the scalar candidate for the translated cube concentration
of \(w=|x_1|^{-\gamma}\), \(\beta=1-\gamma\):
\[
M=a^\beta+(u-a)^\beta,\quad D=a^\beta+(1-a)^\beta,
\quad K=M'D-MD'.
\]
They numerically bracket the root in \(b=a/u\in(0,1/2)\).  A stored sign
bracket is **not** an interval-arithmetic certificate and proves neither the
global concentration formula nor uniqueness.

`flat_concentration_batch.py` uses 120 decimal digits on the \(u\to0\)
branch and \(120+\lceil k\log_{10}2\rceil+20\) digits (142--218 here) on
the \(u\to1\) branch; it records endpoint signs divided by \(1-u\).  It
preserves an earlier strict numerical bracket when a finite-precision midpoint
evaluates to zero.  It also labels an unresolvable asymptotic remainder rather
than calling it a pass.

Run from this directory:

```sh
python3 flat_concentration_batch.py
python3 fraction_half_certificate.py
wolframscript -file flat_concentration_symbolic.wl
```

Then run `flat_concentration_independent.m` in MATLAB.  MATLAB intentionally
skips near-one values that round to one and endpoint cases without a strict
double-precision bracket.  The earlier half-power check using an absolute
\(10^{-25}\) tolerance was too loose at tiny \(u\); the final Python code
checks the root in the dimensionless \(b\) scale and checks \(\Omega\) with
a relative \(\min(\Omega,1-\Omega)10^{-80}\) tolerance.  Full tables and
cache are reproducible local outputs and are not versioned.

`grid_lower_certificate.py` is separate from the scalar-root search.  It
checks exact parameters in the Cartesian-grid lower construction, including
\(\Lambda_{1/5}^5=4096/3125>1\), and the finite \(n=2\) instance.  This is
evidence for the rigorous lower construction supplied by the analytic proof;
it is not an exact-halo numerical search and did not arise from the earlier
half-power scalar diagnostic.
