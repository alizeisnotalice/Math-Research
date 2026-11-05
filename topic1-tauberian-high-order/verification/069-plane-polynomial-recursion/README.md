# R069 exact rational parametric plane recursion

This is the specified constant-capacity 2D min-cap model, not a proof of the nonlinear macro geometry or all-order Cn expansions. `plane_curve.py` uses only Python standard-library Fraction; reproduce e.g. `python3 plane_curve.py --delta 1/20 --depth 5` or append `--save-curves` to retain all final pieces' r,h,P polynomial coefficients (ascending powers of local t).

Each piece stores increasing r(t), h(t)=H(1-r(t)), and P(t)=∫_0^r h. The low map uses q=h and M=q²/2+q(1-r)+P. The high map uses q=1-r+h and M=2A-P+h²/2; each high piece is reversed in t and the list is reversed. New h=delta(q-M), and each new P is rebuilt by integrating h dq, with exact preceding-piece prefix. Every join's r/h and global endpoints are asserted. Mathematica `branch_identity.wl/log` checks the meeting, endpoint and direction identities under the stated derivative hypothesis.

`summary.json` records the final exact runs d=1/20,1/40,1/80 through j=5. All finite piece derivatives and the inequalities h>=0, h<=delta r pass a Bernstein-coefficient nonnegativity certificate; this is a sufficient finite input test, not an all-delta induction. At all three inputs I1 and I2 match R067 exact rational results. j=3 coefficients for all eight parameter pieces are retained in `plane-1-20-depth3-curves.json`; other depth outputs preserve exact rational moments and run resource metrics. Decimal gap ratios in summary are only convenience views of exact rational data, not a fitted proof.

Earlier probe output was generated before Bernstein audit and curve-export fields were added; those additions did not alter recursion or moments. Current depth5 outputs use the final auditable script. No artificial time limit is imposed.
