# R026 proposed finite root/candidate contract — await explicit main release

Theory clone baseb135609efca22428ac8f6a9c9485d114a4a03846, ownership only upper/round026/buffer_parameter_theory. No mathematical numerical calculation/MATLAB/commit/push by theory owner, no changes to prior rounds or lower branch. These files do not automatically authorize execution; main must release after source locking.

Functions and precision are frozen in parameters.json: n7, G(t)=t^7/(t+1)^5-1-7logt, F0(t)=[(1+1/t)^(7/2)+sqrt(1+7logt)]²,512-bit Arb with exact rational inputs. G may form its rational term exactly before converting to Arb. For F0use r³sqrtrwithr=1+1/t, positive sqrtL; the identity A+L+2sqrt(AL), A=r^7, is an internal algebraic cross-check, not an extra scientific row.

First hash/read R025 result.json, terminal.jsonl and manifest. Require statuscompleted, winner.t='11/2', and the terminal record id1matching that row's full interval/grid endpoint exactly. Reuse the saved F0(11/2)row verbatim with provenance; do not reevaluate that F0. Failure of this source gate stops rather than triggering a fallback recomputation. The new G(11/2)first-midpoint evaluation remains authorized and distinct.

Future fixed sequence:

1. Evaluate G(5),G(6), requiring strict negative/positive signs. These match the independent analytic bounds in proof(5). Unresolved sign stops.
2. Exactly80midpoint sign evaluations with exact Fraction endpoints/midpoints. Record each t, signed outward interval and updated bracket. An enclosure containing0stops immediately, with no automatic retry, changed bracket or precision expansion.
3. Take the exact midpoint t_cafter80successful steps, record final width2^-80and midpoint distance bound2^-81, and evaluate exactly one F0(t_c). Round its exact outward upper endpoint strictly upward on the10^-12grid: floor(U*10^12)+1over10^12.

At most83new scientific function evaluations are in scope, plus one old F0row accurately reused. The old row's same input appearing as a new Gsign must be counted by functional, not misreported as a second old-Fcalculation. Do not evaluate E, other t, a continuous optimizer, extra bracket-value bounds, kernel phase or MCsample. Log/sqrtintermediates and actual library call counts should be explicit.

For a claimed strict improvement compare the fresh candidate's actual upper endpoint with the reused F0(11/2)full lower endpoint, never a display or grid decimal. A nonseparated strict comparison is reported without new rows; the fresh candidate remains a legal upper expression if its own certificate is sound. The reused row remains a valid independent existing bound.

Archive source/program prehash, exact sign trace with83cap, raw candidate dyadic endpoints, exact rational ceiling and receipt. A checker may repeat the same fixed inputs as verification; it may not add scientific parameters. Report scalar CPU/wall/RSS without a CPU-saturation claim. Stop for main review after the finite contract, and distinguish parameter-family uniqueness from actualCnoptimality/lower bounds.
