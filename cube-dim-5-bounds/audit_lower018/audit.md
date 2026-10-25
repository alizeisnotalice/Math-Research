# Independent audit 018: the all-scale two-atom ceiling

## Verdict and scope

The proposed bound in `lower_round012/structural_obstruction.md` passes this
audit. For every dimension n >= 1, every alpha > 0, and every nonnegative
measure supported on at most two points with positive total mass,

    alpha |{M mu > alpha}| / mu(R^n) <= 3/2.

The original statement with distinct points and positive masses is therefore
valid. Coincident points and zero masses are harmless extensions. This is
an exclusion of a particular lower-bound input family. It does not give an
upper bound for arbitrary L1 functions, general bumps, or three-atom inputs.
No new C5 lower-bound candidate is certified here.

Target prose SHA-256:
`d3935bccb867e0652806fbcb847a05624f52ab1436de4571d560b7952c3e3c58`.

## Independent proof by homothety

This proof does not use the coordinate residual identities in the target.

First fix the exact level-set representation. For each positive-mass
nonempty subset I of the two atom labels, put

    s_I = (sum_{i in I} w_i / alpha)^(1/n),
    T_I = intersection_{i in I} {x : ||x-p_i||_infinity < s_I/2}.

Then `{M mu > alpha}` equals the union of these at most three open boxes.
Indeed, if x belongs to T_I, choose a radius strictly between the largest
distance to its selected atoms and s_I/2. The centered cube captures at
least their mass and has average strictly above alpha. Conversely, an
actual cube with average strictly above alpha captures a nonempty label
set J. Its radius r satisfies `2r < s_J`, and every captured point is at
distance at most r from x, so x belongs to T_J. This argument permits an
extra atom in a cube selected using I. It also works when both labels
occupy the same position. Use of closed receiver boxes changes only a
finite union of zero-measure boundaries; a box of zero width has zero
volume and its open version is empty.

Relabel the heavier atom as p0, translate it to zero, and independently
reflect coordinate axes so p1 = d with d_i >= 0. Let

    W = (w0+w1)/alpha, R = W^(1/n),
    L = (w0/alpha)^(1/n), lambda = L/R.

The singleton receiver cubes have total volume W. Their union S therefore
has volume at most W. The two-atom receiver box B has coordinate intervals
`[d_i-R/2, R/2]`, when these are nonempty, and volume

    V = product_i (R-d_i)_+ <= W.

If some d_i >= R/2, its corresponding side is at most R/2, and every other
side is at most R. Thus V <= W/2, including the empty-box case, and

    |S union B| <= W + V <= 3W/2.

Otherwise every d_i < R/2. Write Q = [-R/2,R/2]^n and let H be the heavier
atom's singleton receiver cube. Then B is convex, contains zero, and is
contained in Q. Since 0 <= lambda <= 1,

    lambda B is contained in B,
    lambda B is contained in lambda Q = H.

Consequently `|B intersect H| >= lambda^n V`. The heavier mass is at least
half the total mass, so `lambda^n >= 1/2`. The volume added by B to S is
therefore at most `(1-lambda^n)V <= W/2`, which proves the claimed bound.
Closed boxes were used for these volume statements only.

This proof treats alpha as arbitrary. It is consequently an all-scale
argument, with no discretization of alpha or of the displacement.

## Degenerate cases and separation consequence

- Unequal masses enter only through the exact ratio `lambda^n =
  w0/(w0+w1) >= 1/2`; no equality of singleton sizes was assumed.
- Any number of zero displacement coordinates is allowed. If all vanish,
  the two atoms coalesce and the level ratio is exactly one.
- If one mass is zero, discard it; the ratio is exactly one. Including its
  positive-mass pair box does not enlarge the singleton receiver cube.
- A coordinate equal to R/2 belongs to the first branch. A coordinate equal
  to R makes the pair box null; one larger than R makes it empty.
- At equality of singleton or pair-box endpoints, open and closed box
  descriptions differ only on zero-measure faces. The strict inequality in
  the maximal level set is preserved by the open-box representation.
- For a finite collection of one- or two-atom groups, the target's
  separation condition is sufficient: a witnessing cube has side strictly
  below `(Wtot/alpha)^(1/n)`, so it cannot capture atoms from two groups at
  that sup-norm separation. Subadditivity then gives the same 3/2 ceiling.
  The finite-total-mass argument does not settle large cross-group cubes
  in an infinite periodic measure.

## Actual finite exact checks

The registration fixed all fourteen inputs before the only execution.
They cover equal and unequal coincident atoms, a remote zero-mass atom,
the R/2 threshold from both sides and at equality, the R zero-width
threshold and its exterior, mixed zero and nonzero coordinates, extreme
mass imbalance, an exact algebraic singleton contact, and mass/level and
spatial rescaling. Full rational inputs are in `registration.json`.

The new program uses integer fifth-power comparisons at denominator 2^64
to bracket every required radical. It constructs rational boxes contained
in and containing the actual receiver boxes, including the algebraic
displacement in the singleton-contact case. It then partitions all axes
at all rational endpoints and sums disjoint cells using exact midpoint
membership. A separately implemented at-most-three-box inclusion-exclusion
calculation agrees for both the inner and outer unions in every input.
This code does not import or execute the target's programs.

All fourteen exact upper ratio bounds are at most 3/2. The largest displayed
upper endpoint is 1.0551002592053584 (equal masses with singleton contact).
Those decimals are for display only: every pass/fail comparison uses
`fractions.Fraction`. The coincident and zero-mass tests give the exact
interval [1,1]. Actual elapsed time was 0.659 seconds on Python 3.15.0 at
`/opt/homebrew/opt/python@3.15/bin/python3.15`. No MATLAB result is claimed
by this audit.

These finitely many checks corroborate implementation and degeneracy
handling. The all-scale result follows from the homothety proof, rather
than from the finite sample.

## Full standing C and independence

The full standing record was read from the main clone's
`cube-dim-5-bounds/lower_round003/a2_equals_a3_results.json`, SHA-256
`f385811e5460af63abc117804965461eed88d8ecb5b9d2cd30268387646f45d2`.
It contains the algebraic coefficient table, root-isolation data, and the
full rational ratio interval. Exact Fraction comparison verifies

    3/2 < C_lower_endpoint <= C <= C_upper_endpoint.

Thus this entire two-atom family cannot improve the standing C. The old
passed C computation was not rerun. Subsequent lower candidates must be
compared against that full record, rather than a rounded display value.

The auditor used a dedicated sparse clone at
`/Users/zhengzhihao/Desktop/T/Math-Research-cube5-lower-audit-cont`, based
on HEAD `3c73f4e7cd65d2f7671001be41de1f228f155b1c`. The main and candidate
clones were read only. All audit changes are under `audit_lower018/`.
Only the target prose, repository instructions, skill instructions, and
the full standing C record were read for mathematical inputs; the target
scripts/results were not read, imported, or executed. Independence means
an independently derived proof and independently written volume engine;
it is not a claim of blinded access to the theorem statement.

No Git commit or push was made, as requested by the assigning agent. The
audit clone is otherwise clean. Necessary files are this report, the
registration, the independent checker, its exact result, and the hash
receipt. There are no temporary search outputs to clean.

Reproduction:

    cd /Users/zhengzhihao/Desktop/T/Math-Research-cube5-lower-audit-cont/cube-dim-5-bounds/audit_lower018
    python3 independent_two_atom_check.py

The script validates the locked source and code hashes before executing.
No unrestricted search, general L1 upper-bound test, new lower-bound
candidate search, or post-breakthrough extra round was performed.
