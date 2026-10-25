# A015 finite positive-band Taylor calibration and cover interface

Only this directory may change. No new kernel, gradient, MC, grid, MATLAB,
agents, stage, commit or push. All old sources remain frozen.

Derive the positive band B_(a,T)=H_a-H_T <=K_T and its whole-cell moment
upper bound, negative Hessian and simplex Taylor lower bound. Body32 first;
the augmented scalar-Gram extension may be proved without numerical transfer.
Give a full global-cover acceptance contract with exact closed simplex
partitioning, same-plane tube leaves, source identity checks and no omitted
leaf. Do not assume chamber convexity or disjoint plane tubes.

Actual numeric budget: one old center only, T=4/125, a in
{1/1000,1/500,1/250,1/125}, and r in{1/250,1/125}, exactly eight cases.
Read A007 H/grad values and B010 H_(4/125)/grad at the identical center.
The calibration cell is the center's max-norm r box, clipped to x1=1 and
[0,1]^5; ordering is checked exactly. It is a convex cell with32 vertices.
Use only existing dyadic time endpoints between a and T. No infinity tail.

Primary matched target is m0=T/4000=1/125000, matching B009's old pure-Phi
radius991/500000=.001982. Also record the bare64 thresholdT/4096=1/128000.
Only test the two prescribed radii; do not bisect or optimize a new radius.
T=.032 has already been refuted elsewhere by B010: these local calibrations
cannot revive that T or produce a new C6 bound.

At most two extra exact simplex geometry tests: root vertices/strict common-
plane trap, and longest-edge midpoint splitting with exact union and volume
identity. Twenty-minute maximum. Stop after the eight cases and two tests,
including if every tested case fails; retain the entire global coverage gap.
