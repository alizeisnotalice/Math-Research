# A027: original potential-cost feasibility, not ruled out

**Conclusion:** existing evidence does not strictly exclude the original
body32 Green positive-part family from improving64. It excludes several
specific envelopes and one earlier g target. The actual C6 upper bound
remains64, with no new general-L1/all-height numerical certificate.

The accurate definitions and all-scale transfer are restored in
mechanism.md. In particular I is the integral of G(1,y)^(6/5) over the
FULL [0,1]^5 face, B=I^(5/6), and the ideal positive-part cost is
(6B/g)^(6/5). Competition is exactly B<(16/3)g. It is not the short-time
average cost T/(64 min K_T); failure of fixed short-time candidates does
not transfer to this nonlinear spatial truncation family.

## Strict new diagnostic

B023 and F023 independently agree on Gamma for
x0=(50,50,50,45,27,0)/50. This audit only postprocessed that saved rational:

    G(x0)=Gamma/(128 pi^5) ~ 1.75407245845747e-5,
    1/60000 < G(x0) < 1/55000 < 1/50000.

Off-S2 continuity makes g<=G(x0) an essential-infimum statement. Therefore
**g>=1/50000 is impossible**. The target g>=1/60000 survives this test but
has no global proof. The exact necessary I ceiling for any competitive
model is [Gamma/(24 pi^5)]^(6/5). No existing certified I lower bound has
been shown to reach this ceiling.

At g0=1/60000, the precise competition gate is

    B<1/11250,
    I<(1/11250)^(6/5) ~ 1.37599537155658e-5.

The existing short-time norm tail at epsilon=.001 is bounded by
 t0=(6/7)243^-1/6 pi^-11/6 (.001)^7/6.
If the two GLOBAL premises g>=1/60000 and ||H_.001||q<=1/14000 were proved,
then I<=(1/14000+t0)^(6/5), near1.29916954919234e-5, and the cost would be
near60.4266939170373. Both global premises remain unproved. This is a less
restrictive H-norm target than the older1/16000, whose conditional cost
was already52.8688; it is not a new C6 bound.

## Why the existing numerical record does not decide the family

The old conditional-coarea estimator with costs171–1434 used random/QMC
parameter integration and failed a variance stop; it was not an exact
fan evaluation. The later32-point generic G pilot and512-point H_.001
pilot DID evaluate all967680 fan/group terms at each input, but in
binary64; isolated fmpq/Arb reference points calibrate implementation,
not integration. The actual 512-point formal H norm is5.76520e-5. The
243-node Gauss functional has a rigorous interval and formal norm
6.95476e-5, but no closed remainder. Neither is an I bound. The malformed
early denominator257 generic-only pilot is explicitly unusable.

The true whole-face envelope B<7/5000 and the spectral global lower
G>=approximately6.35027e-8 are valid but noncompetitive. Local A006/A007
boxes do not upgrade that global g. B008's off-S2 and singular-tube
budgets are also complete analyses of particular upper allocations;
their failure is not a lower bound on actual B. The exact distinctions,
old point indices and source hashes are in pilot_inventory.md and
reuse_inventory.json. No old successful pilot or failed grid was rerun.

## Most informative next batch

Prefer a finite **G-specific** falsification pass before expanding I
sampling. The existing C013 generic Green/gradient engine already has
its old-point and independent derivative validation; reuse it. B015 z
is a suitable proposed seed because its nonzero-node preflight passed,
and B015 only evaluated short K, p and their gradients. Subject to the
main coordinator's later-work duplicate check, one new G+gradient pass
there and at most one preregistered rational descent point can challenge
g>=1/60000. Stop on failure or after the second point. Surviving points
would not establish a minimum or global g. The old A007 G-gradient point,
B010 exact-G points and B023 x0 must not be re-evaluated as new baselines.

A new H face batch should first meet a proof gate using the saved Gauss Q:
a one-sided remainder of at most roughly3.33571e-7 would prove the relaxed
H target; the existing mixed-moment theorem requires W<=roughly2.25120e-6.
These allowances were certified by scalar postprocessing, but no such
remainder is proved. B006–B008 generic derivative/tube envelopes already
fail; repeating their grids or recomputing their points is not useful.
next_batch.md provides the exact acceptance conditions and missing terms.
No next batch has been launched.

## Execution, failures and limits

Registration/hash preceded all new scalar arithmetic. There were zero
new kernels, native calls, positive-time values, grids, optimization,
MATLAB calls or subagents. The12-item scalar cap was exhausted conservatively:
first attempt computed4 items then stopped on a missing source path
component; second computed5 then stopped before item6 arithmetic because
fmpq cannot parse a decimal string. A final3-item script retained the
Gamma/G enclosure, exact g0/I gate and relaxed-target arithmetic. Total
actual scalar items:4+5+3=12. No further scalar computation was performed.
Unsaved failed-attempt values are not claimed as retained evidence.

results.json contains the three retained Arb512 records and the failure
history. final_postprocess.py is their sole working reproducer; the
obsolete failed script was removed with its hash retained. The original
sources stayed read-only, and no stage/commit/push occurred. The final
receipt pins every deliverable and the read-only inputs. A027 is complete
and frozen; the main agent retains final proof and integration responsibility.
