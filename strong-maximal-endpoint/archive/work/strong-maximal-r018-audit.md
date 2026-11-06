# r018 independent theory audit — 2026-11-06

Scope: read-only audit of `rounds/r018/report.md`; no repository/Git changes.

## Verdict

Accepted. The new result is correctly limited to the ramp family, not promoted
to a general (C_n) equivalence.

## Constant and remainder audit

* `p=(mu-1)/mu` gives
  `kappa_3=2+mu^3*p*(1-p)*(1-2p)=mu*(3-mu)`. Thus the CDF correction has the
  stated sign and coefficient.
* At `x=-z+Delta/(sigma sqrt(n))`, Gaussian displacement is `Delta*f_G` and
  the Edgeworth correction is `a*(1-z^2)*f_G`, with `a=(3-mu)/6`. Thus
  `Delta=a*z^2+b` leaves `(a+b)*f_G`. All quoted errors are `o(f_G)`, since
  the absolute Fourier-CDF remainder is `O(1/n)` and `f_G~z^2/n`.
* The Fourier proof works on the full inner range: the cubic phase inequality
  is global, the quartic real part is absorbed by small epsilon, and the
  outer `/|t|` integral is exponentially small after extracting half of the
  continuous-exponential factor. Hence (4) is a valid absolute `O(1/n)` CDF
  expansion.
* The positive-mgf expansion retains the needed quadratic term:
  `U-1=beta+(mu/2-a)*f_G+o(f_G)`. The negative-mgf overshoot is
  `F-G<=mu/d*exp(-d^2/(2*n*mu))=O(z/n)=o(f_G)` at `d=d_n-O(log n)`.
  Hence `A-1=(mu/2-2*a-b)*f_G+o(f_G)` and
  `b0=mu/2-2*a=(5*mu-6)/6`.
* The dominated-convergence proof gives `I_n~sqrt(pi*n/2)`. The fixed-mu
  score has log coefficient `a-1/2=-mu/6`; adding r017's `(mu-1)/8` shift
  gain gives `-(mu+3)/24` exactly.

Actual Wolfram 15 algebra reproduced all four displayed constants. The
previously live 60/100-bit Wolfram output is now present for `n=256,1024`;
it is numerical corroboration only.

No blockers found in the Fourier/CDF estimate, constants, overshoot, or
the (I_n) asymptotic.
