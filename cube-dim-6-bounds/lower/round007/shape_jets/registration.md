# Round 007 D: registered whole-union shape jets

This follows the completed round006 joint dimer screen. The only owned
directory is `lower/round007/shape_jets/`. No MATLAB, agent spawning,
commits or pushes. E separately audits the interface without duplicating
the direction scan.

Use exact arithmetic in the existing Q(u,v) representation, u^6=2,v^6=3,
for constant and first-order endpoint coefficients. Start from the exact
R6-003 a. Active axis is 3, with active lengths 1..8, inactive lengths
1..4 and the same central inactive starts as round006. Active starts are
all integers q=-4..6. Since a_min>1, r<=2 for ell<=4, a start outside
this range has a strict baseline gap from H=[0,2a_i] times the other
positive half cells and cannot enter for sufficiently small t. Active
lengths 5..8 are to be excluded only after exact checking
`((l-1)*a_min)^6>l*4^5`; those boxes remain strictly empty near t=0.

Perturb b_j=a_j*(1+t*h_j), p_i(k)=b_i*(k+t*E*(k mod 2)), weights
1+t*D*(-1)^k. The period source mass remains 2. One common witness radius
is the sixth root of selected mass divided by 2. The reflected union ratio
is always 16|U_H|; no density constraint is imposed on h here.

The 16 deterministic directions are + and - each of the six h_j, then
+ and - D (mass), then + and - E (gap). There is no optimizer or random
sampling. All baseline zero-width boxes are examined; a first-order
positive width is retained, including newborn thin boxes. No baseline
96-box containment representation is reused.

First run exactly four analytic self-checks: a separating interval touch,
a merging interval gap, a newborn one-thin-direction box, and two newborn
thin directions whose first-order area vanishes. Stop on any failure.

The separated-jet correctness theorem and a separately stated Lipschitz
extension for unresolved higher-order ties are in `theory.md`. Equal jets
do not assert equal true endpoints. Analytic germ keys detect and report
distinct endpoints sharing a jet. Such ties prevent claiming a complete
true endpoint order or an admissible finite step; they do not alter the
first-order volume once the Lipschitz lemma is accepted.

Limits: at most 16 directions after the four self-checks, globally
1,000,000 recursive union states and 1200 seconds. Stop at the first
strictly positive exact derivative, exporting coefficients/direction for
one subsequent small rational-step inward certification; do not perform a
step-size grid. Stop on arithmetic ambiguity, calibration mismatch, or
resource cap. Read/hash the current full ledger, frozen thresholds, exact
baseline certificate, arithmetic source, round006 report, and these sources
before any computation. Results record PID, time, node counts and hashes.
