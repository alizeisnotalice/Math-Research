# Round 030: bounded stationarity diagnostics

The existing `revision_v3` formula family and certificate are the input. No operator lower bound or new upper-bound mechanism is asserted. Python 3.15, SciPy 1.18.1 and mpmath 1.3.0 were used; no MATLAB session was touched. The computer reported 15 logical/physical CPUs and 24 GiB RAM; existing MATLAB processes were idle in the sampled `ps`. Six workers were used, each single-threaded in BLAS, with no persistent cache. The resource counts, PIDs, hashes and stops are in the JSON files. On macOS `ru_maxrss` is reported in bytes. Each dimension used 41–55 combined float objective/derivative evaluations plus one high-precision diagnostic evaluation; the limit is 128 float evaluations. The fixed 32-point bracket scan is new, not a replay of the earlier optimizer trace. It does not establish uniqueness or global minimization.

`stationarity.py` uses positive quadrature for B and J, exact analytic derivative expression, sign-change bracketing, and 60-digit positive quadrature at each float root. The t0 formula n/(log n - 3 log log n) is evaluated only when the denominator is positive and t0 > 1. At n=4,256,1024 respectively it yields approximately 9.8427, 629.9343, 911.6360; it is unavailable for 7,16,64. The computed local minimum roots have n/t approximately 1.1133, 1.0820, 1.2017, 1.5836, 2.1222, 2.7835 for n=4,7,16,64,256,1024. Their (F - n log n + n log log n)/n is approximately 2.3640, 2.2277, 2.0128, 1.7096, 1.4739, 1.2899. These are numerical diagnostics of F at located stationary points, not certified U values or asymptotic proof; no inference of global minimality is made. Large-n integrals use positive integrands to avoid cancellation in the signed polynomial expansion.

The n=7 root suggested t=808656059921/125000000000. `certify_n7.py` invokes the existing exact rational enclosure and finite-a routines. It checks the complete rational inequality F_new <= F_upper < 24.554960525964662 < F_old_lower and a finite-a bound with a=2^256 below the same number; see `certificate_n7.json` for endpoints and exact strict separation. The numerical improvement over the earlier point is approximately 3.87e-14 and is only a refined parameter in the same family, not a new mechanism. The separate original theorem is needed to translate this certificate to the original operator bound.

Reproduce from the clone root:

```
PYTHONDONTWRITEBYTECODE=1 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-dim-7-bounds/upper/round030/compute/stationarity.py
PYTHONDONTWRITEBYTECODE=1 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-dim-7-bounds/upper/round030/compute/certify_n7.py
```

The diagnostic JSON contains run-specific PID and timing fields and will change on reproduction. No output values alone prove the global infimum or a new asymptotic expansion.
