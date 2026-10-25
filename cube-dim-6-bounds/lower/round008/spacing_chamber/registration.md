# Round 008 D: second-spacing chamber and next actual raw endpoint event

Current comparison is the complete E007 ideal-union upper fraction and
36 field coefficients read from E's `audit/lower007/candidate_number_field/
candidate_results.json`. The displayed strict bound 2.06692413818 and the
stale main ledger are not the dynamic comparator. No replay of the already
certified 1e-6 step is performed.

Write s=a2, with all other R6-003 spacings and unit masses fixed. Start at
`s_current=(1000001/1000000)*(u^3+u^4*v^3)/4`. Rebuild all 4096 standard
length vectors in {1,2,3,4}^6, retaining all raw boxes that are nonempty at
s_current. No 96-box containment reduction is used for event analysis.

For the moving coordinate, an odd length l has interval
`[0,min(s/2,r-(l-1)*s/2)]`, and an even l has
`[max(0,l*s/2-r),s/2]`. For l>=2, positivity is equivalent to
`(l-1)*s<2r`; for l=1 it always holds at positive s. Therefore increasing
s cannot create a previously empty standard box. Boxes excluded by other
coordinates are empty for the entire parameter line. All necessary births
are thus settled analytically; there is no cap4-to-cap5/6 rerun.

Make one finite algebraic inventory of pairwise intersections of all
unclipped moving affine endpoints and cell faces for the retained raw
family. These endpoints have rational slopes in the absolute parameter s,
so every event lies in Q(u,v) with no general field inversion. Examine at
most the 12 nearest distinct events greater than s_current. At each,
compare actual clipped endpoint choices and their lexicographic left/right
orders. Skip events of invisible expressions that do not change the raw
arrangement. Stop at the first actual raw-arrangement event s_star.
This need not be the first kink of the final union volume; that distinction
will be explicit.

Normalize the chamber by `s(t)=s_current+(s_star-s_current)*t`, 0<=t<=1.
Export all raw boxes with affine field endpoints, exact ordered tie groups,
and adjacent gap checks at t=0,1. Persistent ties mean identical affine
functions; boundary contacts are allowed. Only coordinate 2 varies, so
the complete raw-union volume is affine throughout the chamber closure.

At most ONE new event parameter is evaluated with the full static exact
union. Compare its complete lower fraction with E007's complete upper.
If it passes, stop immediately and hand the candidate and chamber interface
to E. Do not scan subsequent events or optimize beyond the first promotion.
If the event inventory exceeds the 12-examination cap without an actual
event, report the obstruction rather than extend the list.

Hard limits: 12 examined event candidates, 1 full candidate union,
1,000,000 recursive union states, 1200 seconds globally. Stop on any
unseparated field comparison, invalid order proof, resource limit or
calibration disagreement. Arithmetic uses checked rational root isolation
at denominator 10^80. Hash E007, the old field code, round007 derivative
certificate, this registration and current source before execution.

Only `lower/round008/spacing_chamber/` is writable. No MATLAB, subagents,
Git commits/pushes or main-ledger edits. The final result records all input
hashes, counts, PID, execution handle and finished status.
