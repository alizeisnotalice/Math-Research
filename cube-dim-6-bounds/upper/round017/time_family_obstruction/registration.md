# C017 registered theory and scalar postprocessing only

Only this directory is writable. All earlier rounds and B015 inputs are
read-only. No subagents, MATLAB, stage, commit, push, full fan traversal,
new kernel call, new spatial point, new G value, or new T evaluation.

The exact spatial point is
(100000000000000,99999998999999,99999997999988,90268521812170,
53714088833376,157107163435)/100000000000000.
The existing B015 certificate supplies K and p at T0=6/125.

One scalar batch, at most 1200 seconds and 12 named rational checks:
read exact dyadic endpoints and verify their provenance; check positive
input intervals, the existing q(T0)<1/4096 and p(T0)<q(T0), the positive
coefficient A in the time-envelope formula, and its two rational endpoint
values at a=223/5000 and b=251/5000. Preferred common bound is 1/4100;
the fallback is 1/4096. Both choices are fixed now. No endpoint is adapted
after seeing results. Every named check and its exact result is retained.

The underlying time-envelope theorem, probability-mixture implication,
and an all-T rational root-certificate criterion are analytical work.
No all-T conclusion is permitted without its missing global certificate.
Any proposed exact G or rational p(T) root/interval computation requires
a separate future release; none runs in this batch.
