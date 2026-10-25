# A007: one-center time-moment batch

Continue body32 g certification; no new project and no global C6 claim.
Use only the owned directory `upper/round007/g_time_moments/` in the theory
clone. Other workers are active; no other files are edited, no MATLAB,
commits, pushes, or recursive agents.

The single center is the exact rational point from round003
`local_minimum/exact_local_reference_results.json`, denominator 1000000007.
That source supplies a G value, but no H_tau values at this center. The eight
time endpoints are tau=0.001*2^k, k=0,...,7. Eight new H nodes are permitted.
At most five of those same center/time nodes (tau<=0.016 and H_tau>1/60000)
also receive exact real fan gradients. There are no extra spatial nodes.
Thus even counting each gradient calculation as another evaluation, at most
13 evaluations are run, below the 16-node cap. Wall-time cap: 20 minutes.

Before kernel execution, calibrate the real divided-difference gradient on
analytic scalar anchors, including a confluent node-derivative identity.
The kernel gradient is real automatic differentiation of atan divided
differences; repeated nodes use the copied-node DD6 identity. No complex-step
operation is applied to a real-part formula.

For each epsilon, construct center time-band masses Delta H. Arb differences
are intersected with [0,infinity) in their rigorous endpoint bounds. Whole-box
J2 uses each band's own exp(r/tau), and the tail uses T^-2 exp(r/T) H_T.
The lower Hessian inequality D_h^2 H >= -||h||_2^2 J2/16 and the directly
computed center gradient give an actual continuous Taylor lower certificate.
The previously proved upper inequality H_ii <= (4225/4224) J2 is also recorded.
Actual clipping to [0,1]^5 is used in the Taylor linear term.

Stop after this finite batch, report accepted radii, exact arithmetic evidence,
and uncovered complement, then await primary scheduling. Do not rerun the
round006 eight boxes or any existing H baseline.
