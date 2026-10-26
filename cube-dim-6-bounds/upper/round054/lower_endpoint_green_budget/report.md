# A054: full-time Green budget for all31 lower endpoint tubes

## Verified result

For the original body32 Cauchy frame, define G=integral_0^infinity p_t dt,
q=6/5 and rho=1/64. For each nonempty proper A containing coordinate1,
let B=A complement and

    E_A={x=(1,y), y in[0,1]^5:
          dist(x,span(1_A,1_B))<1/64,
          mean_B(x) in[0,1/4]}.

The registered complete positive envelope gives

    ||G||_(L^(6/5)(union_A E_A))
        <10293038485441/500000000000000000
        =0.000020586076970882.

This is a full positive-time bound on all31 specified low-mean tubes.
It is not only a short-time estimate and uses no newly evaluated kernel.
The right side is strictly less than(251/1000)*(4/48750). The computed
upper-budget ratio is approximately0.2508928130826033; this is a budget
comparison, not an estimate from below of the actual norm.

No high endpoint u>3/4, middle-domain long-time tail, outside-tube face
region, global Green floor or full I/g condition is paid by this result.
There is no new C6 bound; the current upper bound remains64.

## What changed mathematically

At a lower endpoint1_A, several signed pair planes can meet. The proof
classifies them completely: an incident plane must have support blocks
exactly A/B, with constant A signs and arbitrary B signs. There are
2^(m-1) such planes and2^m ordered pairs, all in class min(k,m). This
includes every mixed-B incident plane, which cannot be omitted by using
the two-term middle formula.

For nonnegative x, the distance to every incident signed plane is at
least dist(x,U_A). Every remaining plane has the uniform positive gap
3/(4sqrt(2))-sqrt(k)/64. The norm and every body-line distance also have
proved positive lower bounds. Both singleton cases k1 and m1 are
included, without full direction enumeration.

The already audited C037 global Green envelope is then split into all
incident ordered terms and all remaining terms. Every coefficient in
the background sum is a nonnegative integer. This avoids subtracting
correlated interval evaluations of Mtotal and keeps all992 ordered
pairs charged for each domain.

F034's genuine face Jacobian sqrt(km), maximal-coordinate cone fraction
1/k and four-dimensional normal radius give the volume and singular
integral upper bounds. At u0, positivity may clip the normal ball; the
proof only enlarges that clipping. It does not claim that the displayed
volume is the actual volume of E_A.

## Five types and union accounting

Here k is the size of A, m=6-k, n=binom(5,k-1). Decimal entries are only
readable displays of the full saved intervals.

| k | m | n | incident ordered pairs | singular norm budget | background norm budget | total U_k |
|---:|---:|---:|---:|---:|---:|---:|
|1|5|1|32|7.6982957655703e-6|4.3145846665898e-6|1.2012880432160e-5|
|2|4|5|16|6.6979797126855e-7|9.5970697773356e-7|1.6295049490021e-6|
|3|3|10|8|1.4197632185740e-7|4.3740575233891e-7|5.7938207419631e-7|
|4|2|10|4|4.8024996403263e-8|2.8364546673447e-7|3.3167046313773e-7|
|5|1|5|2|2.6685364985527e-8|2.5800538098165e-7|2.8469074596717e-7|

The background class counts are respectively
(160,480,320), (192,464,320), (192,480,312), (192,476,320),
and(190,480,320). Each row, including its incident count, sums to992.
The complete weights are(1,5,10,10,5), summing to31. No k versus6-k
identification removes weights. For the common function G,

    ||G||_union <= [sum_k n_k U_k^(6/5)]^(5/6).

The domains are permitted to overlap; no mutual disjointness is claimed.
The k1 type is the largest individual cost and its32 incident terms are
fully represented. The fixed union formula passes the specified
conditional whole-target comparison, so no radius change or scan was
performed.

## Source and evidence boundary

Theorem.md supplies the new endpoint classification, line/plane gap
proofs, exact face integration interface and positive union bound.
C037's three M values are consumed as their original complete dyadic
intervals. Their defining time integrals are never reevaluated here.

C037 manifest SHA:
ff19f756bbc40113cf28ad94bd5addb1d8e9cd9a5cc458298b14f688566fc911.
F037 final independent audit manifest SHA:
6a89dc230847d57b886adbf85e17189a925a310b27beafffdb5e4ed1bf250792.

C037 results preserve the historical status pending F035. Its frozen
audit_dependency.json discharges that status using the exact F035 source
binding, and F037 explicitly confirms it. This task checks that chain
without modifying or rerunning any prior source. The frozen F034 theorem
and independent geometry review are also pinned in registration.json.

The present computation is an exact interval evaluation of a proved
upper formula. It is neither a numerical approximation of the true
Green integral nor an independent reproof of the old M coefficients.
The parent retains final mathematical acceptance and integration.

## Actual finite execution

Start:2026-10-25T22:22:35Z. Registered deadline:22:42:35Z.
The fixed plan has eight configurations: one exact structural check,
five fixed geometric types, one union, and one pure saved-data audit.
There are zero kernel/native calls, zero old-M reevaluations, zero full
direction/fan enumeration, grids, optimizers or MATLAB calls.

An initial launch with /opt/homebrew/bin/python3 failed at importing
SymPy before run() and before creating any receipt or executing any
scientific configuration (tool chunkf7bdee). No code, parameter or source
hash was changed. Selecting the already installed mathematical runtime
resolved the environment issue:

    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B evaluate.py --run

That single actual producer run completed all seven configurations,
PID54112, exit0, elapsed0.006180208001751453 seconds, chunk502932.
It used Arb768 while preserving the complete Arb512 source-M uncertainty.
All output endpoints and sources are saved, not reconstructed from
printed decimal values.

The eighth configuration, check_saved.py, imports only Python standard
library exact arithmetic. It verifies all source bindings, the original
M intervals, five nonnegative ordered-pair ledgers, square-root distance
gates, and the true radial formula

    J_k^5 = face_factor^5*(5*pi^2/7)^5*rho^14.

The two norm contributions are certified by sixth powers; weighted union
powers are checked by fifth powers. Finally the stated rational upper
raised to6 strictly dominates the certified sum raised to5. The final
upper is also proved below(251/1000)*(4/48750). The audit passed on its
first invocation, chunk1bda52, elapsed0.004104666004423052 seconds.
No new Arb evaluation, M integral or spatial sample occurs in that audit.

All eight configurations are exhausted. This task stops; it does not
initiate another budget or change rho. The package contains only the
proof, fixed producer, fixed rational checker, exact results, receipts
and final hash manifest. No old source was edited and no commit/push
or delegation occurred.
