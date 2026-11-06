# `sm_receiver_alternating30.m` read-only audit — 2026-11-06

## Core transform and gradient: accepted

`logDensity(b,c,s)` is not the probability log-density despite the bottom
comment: its displayed formula omits the Gamma `exp(-s)` factor, so it is
exactly `log(e^s f_{b,c}(s))`.  Let `s_i=x_i-u_i`.  Therefore

`lr+sum(u)-sum(x) = log f_{b,c_1}(s_1)+log f_{b,c_2}(s_2)`.

Consequently `P` is the ordinary probability quadrature weight and its
mass test is the right diagnostic.  Also

`R=exp(lr+sum(u)-tauRef)`

is the correct relative likelihood for the `e^{-sum(x)}` modular measure
(with the common `tauRef` scaling), not a missing-shift error.

The outer gradient has the correct sign and shift.  The identity

`partial_c f_{b,c}(s)=b p_c[f_{b,c}(s)-f_{b-1,c}(s-c)]`

requires the lower-block term to use density argument `x-u-c` while the
probability conversion retains the original receiver geometry.  Passing
`ux=u+c` to `receiverAverageOrder` does exactly this: the extra relative
density factor cancels against its `+ux-sum(x)` conversion.  Hence
`b*p_c*(B-Bold)` is correct.  The explicit recomputation of `Btrial`
provides an additional discrete nondecrease check independent of the
gradient formula.

## Retention/history: accepted with two implementation qualifications

Each receiver ratio affects only its own average while `Hold` is fixed;
thus `oldUnderNew>=minOld` is a valid discrete invariant.  If the inner
candidate is rejected, the retained `H,w` remain associated with the last
accepted `profileSourceCs`; the final JSON records that source separately
from the current receiver `cs`, so it is not a provenance mismatch.

Two minor but real code-level qualifications:

1. On the *initial* `solveOrRetain([],...)` call, `useNew=isempty(Hold)`
   bypasses `scaleOK`.  If cap/scaling failed, the script could start with a
   zero/non-strict profile despite the prose saying that every accepted
   profile is strictly feasible.  The initial call should assert
   `inner.strict_scale_ok && min(P'*H)>=target` before continuing.
2. The inner retention test permits `J_new <= J_old*(1+2e-10)`, so it does
   not literally guarantee nonincreasing discrete cost.  This is harmless
   for a floating initializer but should be described as a tolerance, not
   as an exact invariant.

All conclusions remain only for the finite quadrature grid.  The code
already says that its averages, modular cost, union diagnostic, SQP result,
and 33/65-point comparisons are not continuum certificates.
