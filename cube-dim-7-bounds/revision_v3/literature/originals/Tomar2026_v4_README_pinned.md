# Bounds for the weak type (1,1) constant of the centred maximal operator over squares

Certificates, checkers, search programs and Lean proofs for the paper *Bounds for the weak type (1,1) constant of the
centred maximal operator over squares* by Tejasvi Singh Tomar
(Zenodo, [doi:10.5281/zenodo.23554628](https://doi.org/10.5281/zenodo.23554628)). This repository is archived in the
supplement record [doi:10.5281/zenodo.23554629](https://doi.org/10.5281/zenodo.23554629).

## Theorem

Let M be the centred Hardy–Littlewood maximal operator over axis-parallel squares in the plane, and let c₂ be the
least constant such that λ |{Mf ≥ λ}| ≤ c₂ ‖f‖₁ for all f in L¹(R²) and λ > 0. Then

    1.7027237335 ≤ c₂ ≤ 3.59123.

The constant c₂ is entry 47a of the list of optimization constants collected by Tao and others
(https://teorth.github.io/optimizationproblems/). The covering argument gives c₂ ≤ 4. Yongxi Lin proved c₂ < 3.879,
with a proof checked in Lean (https://github.com/CoolRmal/centered-maximal-constant), and has also obtained
c₂ < 3.615749 with a kernel of the type used here (Section 1 of the paper describes the status of that bound). The
lower bound improves Lin's c₂ ≥ 1.6855099933, also proved in Lean in the same repository; Aldaz's earlier lower bound
is 1.6211. Both of our bounds are also proved in Lean (`lean/upper/` and `lean/lower/`).

## Proof outline of the upper bound

1. **Diamonds.** The map (u, v) ↦ (u + v, u - v) takes the diamond D = {|u| + |v| ≤ 1} onto the square [-1, 1]².
   So c₂ is also the weak type constant of the maximal operator over the dilates of D.
2. **Comparison principle (Lin).** For 0 < α < 2 let

       Af(x) = Σ_{j=1,2} ∫ [f(x + t e_j) + f(x - t e_j) - 2 f(x)] |t|^(-1-α) dt.

   If K is even and integrable, K ≥ 0, K ≥ 1 on D, and AK ≥ 0 away from the origin, then c₂ ≤ (1/2) ∫K.
   Section 2 of the paper gives a detailed proof.
3. **Closed forms.** The one-dimensional operator maps the cubic B-spline to a finite sum of powers |s - k|^(3-α),
   and A maps (|u| + |v|)^(-α) to C_α (|u| + |v|)^(-2α). The truncation of the power contributes one-dimensional
   integrals given by hypergeometric functions (Section 3).
4. **Certificate.** `cert/cert_h32_a4740_R3320_3p59123.json` gives the kernel

       K(u, v) = c Ψ(|u| + |v|) + Σ_{i,j} a_ij B(u/h - i) B(v/h - j),        Ψ(r) = max(r^(-α) - R^(-α), 0),

   with α = 47/40, R = 33/20, h = 1/32, c = 1977601376167/1099511627776 and 1181 nonzero coefficient classes, all
   dyadic. Then (1/2) ∫K = 3.5912243275064523... The format is described in `cert/README.md`.
5. **Verification.** Two programs check the claims (C1)-(C4) of `cert/README.md` in interval arithmetic. As in the
   paper, the first program is `checker/cert_run.py` and the second, written independently, is
   `checker/independent/hl_check.py`.
6. **Lean proof.** `lean/upper/` proves c₂ ≤ 3.59123 in Lean 4, with the numerical checks evaluated by the Lean kernel
   (see *Lean proof of the upper bound*).

## Checking the certificate

Both programs need Python 3.12 and python-flint 0.9.0 (`pip install python-flint==0.9.0`), which provides the Arb
library.

    python3 -I checker/cert_run.py cert/cert_h32_a4740_R3320_3p59123.json --Q 8 --maxdepth 7 --procs 16
    python3 -I checker/independent/hl_check.py cert/cert_h32_a4740_R3320_3p59123.json --claim 3.59123 --procs 16

The summary lines of the outputs are:

    (C3) (1/2) int K in [3.59122432750645234527147802 +/- 2.23e-27]  upper 3.591224327506452
    (C1/C2) K checks: 2628 knot cells, failures 0  t=0.3
    (C4) A K checks: delta = h/8, cells evaluated 219732, failures 0, min cell lower bound 1.346e-08  t=68.1
    PASS

    C3 mass: (1/2) int K in [3.59122432750645 +/- 2.35e-15] ; claim 3.59123 ; PASS
    C1 K >= 0: 1368 knot cells, 2025 leaves, min lower 0, failures 0
    C2 K >= 1 on D: 272 knot cells, 272 leaves, min lower [1.0000993 +/- 4.47e-8], failures 0
    C4 A K >= 0: 1368 knot cells, 443754 leaves, min lower 3.28818e-08, failing sub-cells 0
    OVERALL (c1,c2,c3,c4): PASS

Both programs end with a line that reports PASS or FAIL, and after a failure their exit status is 1. On a shared
384-thread x86-64 server with 16 processes the two checks took 68 s and 33 s. The full outputs are in `logs/`, together
with:

- another run of `cert_run.py` with finer cells (`--Q 16 --maxdepth 5 --only A`), which checks (C4) on 415,236 cells
  with no failures;
- a run of the second program with Taylor models on every sub-cell (`--tight`), whose least lower bound is 3.6e-5;
- the test logs of the second program;
- a comparison of the closed forms with quadrature of the definition of A, made at the smallest grid values of A K and
  at the local minimum 2.8035e-4 near (1.541724, 0.132603) found by a local search.

`logs/controls/` holds the outputs of the audit pipeline `checker/independent/audit_spec.py` on this certificate. These
are tests of the enclosures at its parameters, spot checks of the closed forms by quadrature, and four perturbed
certificates that the second program rejects: `H1` claims a bound below (1/2) ∫K, `H2` has AK < 0 at a known point,
`C1ctl` has K < 0 at a known point, and `C2ctl` has K < 1 at a known point of D. Each of `C1ctl` and `C2ctl` comes
with a run on the unperturbed certificate at the same depth (`C1ref`, `C2ref`), which passes. The first program also
rejects `H2`: with `--Q 8 --maxdepth 6 --Kdepth 4` it reports 1,069 failing cells (`first_H2.log`). In the same
directory, `full_scan.json` (9 × 9 points per cell) and `fine_scan.json` (33 × 33 points per cell) are scans of AK over
all knot cells; `band_scan.json` uses 65 × 65 points per cell near the curve r = R and comes from a short audit script
that is not included. To rerun:

    python3 -I checker/independent/audit_spec.py cert/cert_h32_a4740_R3320_3p59123.json 3.59123 OUT --procs 16
    python3 -I -c "import sys, json; sys.path.insert(0, 'checker/independent'); import audit_spec as A; \
        res, n = A.full_minima_scan('cert/cert_h32_a4740_R3320_3p59123.json', 16, n=32); print(n, res[:6])"

## The first program

`checker/hlcert.py` and `checker/cert_run.py` were written together with the search. `cert_run.py` checks (C1) and
(C2) on the knot cells of the grid h Z², computes an upper bound for (1/2) ∫K, and checks (C4) on sub-cells of the
knot cells in 0 ≤ v ≤ u. It prints the number of failing cells for each claim. On a sub-cell it bounds A K below as
follows:

- the spline part of AK is given by formula (3.3) of the paper; each power |u/h - m|^(3-α) is replaced by its Taylor
  polynomial of degree 2 at the centre with the Lagrange remainder, or, for the two knots bounding the knot cell, by
  its chord or its tangent at the centre, according to the sign of its coefficient;
- c AΨ is bounded below by termwise bounds at the corners of the sub-cell (C_α r^(-2α) and I are decreasing, and
  Lemma 3.6 gives monotonicity outside r = R), or by a Taylor model when the sub-cell lies entirely inside or entirely
  outside the curve r = R;
- the two polynomial parts are added and bounded below by their Bernstein coefficients.

The last line of the output reads `PASS` if no cell fails and `FAIL (n failing cells)` otherwise, in which case the exit
status is 1. The verification is the default run, `--only MKA`, in which M computes the bound for (1/2) ∫K, K checks
(C1) and (C2), and A checks (C4). Runs restricted with `--only` check only the listed parts. The default run prints the
certified constant: (1/2) ∫K lies in [3.59122432750645234527147802 +/- 2.23e-27], with upper bound 3.591224327506452.
The bound 3.59123 of the theorem is this value rounded up, so a default run that ends in `PASS` and prints a constant
of at most 3.59123 verifies the theorem. `cert_run.py` leaves this comparison to the reader; the second program makes
it (`--claim`).

## The second program

`checker/independent/` holds a second program, written separately, from the definitions in `cert/README.md` and the
proofs in the paper, without using the code of the first. Its documentation is in the docstring of `hl_check.py`.

- `hl_check.py SPEC --claim BOUND [--procs P] [--maxdepth D] [--tight] [--out REPORT.json]` checks (C1)-(C4). The
  spline part is handled with exact rational Bernstein nets, and AK with polynomial models of degree at most (3, 3).
  It rejects any input field that would change the kernel (it ignores the descriptive fields `lp_val`, `source`,
  `bound` and `note`), and it records the SHA-256 digests of the files it loads.
- `quad_A.py` evaluates A K at a point by rigorous quadrature of the definition of A, and `e2e_check.py` compares it
  with the closed forms of the paper at given points.
- `test_*.py` test the enclosures used by `hl_core.py` against direct evaluation; `run_tests.sh` runs them and writes
  its logs to `runs/tests/`. `test_c12.py` and `test_guards.py` regenerate the files in `testdata/synth/` and
  `testdata/G3_bad_key.json`, so the tests need a writable copy of the tree.
- `audit_spec.py SPEC CLAIM OUTDIR [--procs P]` runs the whole pipeline of `logs/controls/` on a certificate.

The docstrings of both programs use other names: `DEFS.md` means the format and the claims of `cert/README.md`,
`PROGRESS.md` means the proofs in Sections 2 and 3 of the paper, F1 means formula (3.3) together with Lemma 3.2,
and F2, F3 and F4 mean Lemmas 3.4, 3.5 and 3.6.

## Search

`search/` holds the programs that found the certificate and the logs and intermediate solutions of the runs that led to
it; `search/README.md` describes them. Only the certificate enters the proof.

## Lower bound (`lowerbound/`)

`lowerbound/c2lower_1p70272.json` describes a periodic measure: `lam` (here 1), the periods `Px` and `Py`, and the 50
atoms of the period cell `[0, Px) × [0, Py)` as triples `[x, y, mass]` of rationals; `reach` bounds the sides of the
boxes that the programs consider. The programs also use the rounding unit `D = 2^40`, which is not stored in the file.
By Lemma 6.2 of the paper, `lam |U ∩ C| / μ(C)` is a lower bound for c₂ when `U` is the union of the lattice translates
of finitely many rectangles, each contained in a rectangle `R_S(ρ_S)` of Lemma 6.1 with `0 < ρ_S < √(μ(S)/(4 lam))`.
The command

    python3 -I lowerbound/check_exact2.py lowerbound/c2lower_1p70272.json

enumerates the boxes `[x_a, x_b] × [y_a, y_b]` whose bounds are atom coordinates near the cell and whose sides are at
most `reach`. For each box it takes the set `S` of atoms in the box and `ρ_S = (⌊√⌊μ(S) D²/4⌋⌋ − 1)/D`, so that
`ρ_S < √(μ(S)/4)`, clips the rectangle `[x_b − ρ_S, x_a + ρ_S] × [y_b − ρ_S, y_a + ρ_S]`, which lies in `R_S(ρ_S)`, to
the cell, and computes the area of the union of the 5135 clipped rectangles of positive area exactly, in integer
arithmetic with a sweep line. It takes about two seconds, and its output ends with

    EXACT ratio lower bound = 482167717377842991882283179477574801392641890014291050323/283174368139841919615520431313000857449005056000000000000
    ratio >= 1.7027237335

and writes the same data to `c2lower_1p70272.check2.json` next to the certificate (our copy is in `lowerbound/logs/`).
The command

    python3 -I lowerbound/check_exact.py lowerbound/c2lower_1p70272.json

computes the same number with Python fractions, without the integer encoding, in about ten minutes; its output is
`lowerbound/logs/c2lower_1p70272.check_fraction.log`. These two programs use the same method and need only the Python
standard library.

`NOTE_lower_bound.md` states the two lemmas and outlines the check. It states Lemma 2 for the rectangles `R_S(ρ_S)`
themselves, while both programs use the smaller box rectangles above, which Lemma 6.2 of the paper allows.
`NOTE_addendum.md` documents the certificate format, the boxes and the window that the two programs use, the radii, and
one step of the proof of Lemma 2. In the certificate, the field `source` names the search output from
which it was rationalized (`lowerbound/search/c2lower_1p70272.search_result.json`), and `float_ratio` is the value that
the floating-point search reported; neither enters the check. `lowerbound/search/` holds the search programs and that
search result, together with `rationalize_seg.py`, which made the certificate from it; they are not needed for checking.

`lowerbound/frame_only.json` keeps only the two atoms of the certificate on the line y = 0, that is, Lin's lattice with
the new parameters. For it `check_exact2.py` gives 1.6447274428 (`lowerbound/logs/frame_only.log`), the value quoted in
Section 6 of the paper.

`lowerbound/independent/lb_check2.py` was written from `NOTE_lower_bound.md` without using the code of the two
programs above. It enumerates the sets `S` in another way, with pruning that loses no rectangle meeting the cell, takes
the rectangles `R_S(ρ_S)` themselves, and computes the area of their union exactly from coverage counts on the grid
formed by the sides of the rectangles. It needs NumPy (we used Python 3.12 and NumPy 2.5.4):

    python3 -I lowerbound/independent/lb_check2.py lowerbound/c2lower_1p70272.json 1099511627776 \
        482167717377842991882283179477574801392641890014291050323/283174368139841919615520431313000857449005056000000000000

The last two arguments are `D` and the claimed value. In about ten seconds it prints the same rational number and
`RESULT: PASS`. It reports 3407 rectangles rather than 5135 because it counts distinct rectangles `R_S(ρ_S)` meeting
the cell, while `check_exact2.py` counts one clipped rectangle for each box. The other files of
`lowerbound/independent/` are:

- `cert_lin.json`: Lin's lattice, with his parameters replaced by rational approximations. For it `lb_check2.py` gives
  1.6855099933518 with `D = 2^40` and 1.6855099933555 with `D = 2^60`, while Lin's bound is 1.6855099933555...;
- `lb_check.py`, a version without the pruning that reads `D` and a window size from its input file, and `test_lb.py`,
  which tests it on random sets of rectangles, on small measures with known values and on invalid inputs, and rewrites
  the files in `tests/`;
- `test_lb2.py`, which tests `lb_check2.py`: its union area against inclusion-exclusion, its exact agreement with
  `lb_check.py` on 12 random periodic measures, and its exit status for a claim above or equal to the value. It
  rewrites the files in `tests2/`;
- `logs/`: the outputs of these runs, and in `logs/controls/` runs of `lb_check2.py` that fail as they should: the claim
  raised by 10⁻³⁰, `D = 2^20`, the mass of the atom at (Px/2, 0) changed by 1 % either way, one small atom shifted by
  0.01, and certificates with an atom outside the cell, two atoms at the same point or a zero mass.

## Lean proof of the upper bound (`lean/upper/`)

`lean/upper/` is a Lean 4 project that proves

    theorem Hl47a.weakTypeConstant_two_le : Hl47a.weakTypeConstant 2 ≤ ENNReal.ofReal (359123 / 100000)

together with the strict inequality (`weakTypeConstant_two_lt`), the weak type bound itself (`isWeakTypeBound_two`) and
the form with the level sets {Mf ≥ λ} and real λ > 0 (`nonstrict_two`). The definitions, in
`lean/upper/Hl47a/Statement.lean`, import only Mathlib. They are the definitions of Lin's formalization of his bounds
(https://github.com/CoolRmal/centered-maximal-constant, commit 72c022b) and use the strict level sets {Mf > λ}, which
give the same constant c₂. The four theorems depend only on the axioms `propext`, `Classical.choice` and `Quot.sound`,
and the numerical data are checked by the Lean kernel (`decide +kernel`). In Lean, `#print axioms
Hl47a.weakTypeConstant_two_le` after `import Hl47a.Main` prints these axioms.

The project is a modified version of Lin's formalization of the partial balayage argument for squares
(https://github.com/CoolRmal/partial-balayage-lean, commit 0055cdf5b947f4e23b0b736a297ab38b48b1a270,
Apache License 2.0), in which he proved c₂ < 452/125 with his kernel. Lin's analytic chain is kept, with the constants
of our kernel (α = 47/40, R = 33/20, h = 1/32). The spline part of our kernel does not vanish outside the diamond r < R,
and the check of AK ≥ 0 there is new; Lin's checks of K ≥ 0 and of K ≥ 1 on D are replaced by checks for our kernel,
with K ≥ 0 checked on the whole square max(|u|, |v|) < 9/4. The numerical data, generated from the certificate (the data
inside r < R mostly in Lin's formats, the data outside it and for K in new formats), give lower bounds for AK on 32,027
rectangles that cover the part of r < R in 0 ≤ v ≤ u and on 1,751 rectangles that cover the rest of the 1,368 cells of
0 ≤ v ≤ u on which K is not identically zero, and lower bounds for K on 2,628 cells; the bound (1/2) ∫K < 3.59123 comes
from an exact formula for ∫K. `lean/upper/README.md` gives the credits, `FILES.txt` says for each file except
`README.md`, `FILES.txt` and `SHA256SUMS` whether it is unchanged from Lin's commit (232 files), modified or generated
in the format of one of his files (1,264) or new (808), `LICENSE-Lin` is the licence file of Lin's commit, and
`SHA256SUMS` lists the files.

```
cd lean/upper
lake exe cache get
lake build
```

builds the 2,299 modules. The kernel checks of the numerical data take about 64 CPU-hours, or about an hour on 88
hardware threads. The programs that generated the data files are not included; the Lean kernel checks the data,
so these programs need not be trusted.

## Lean proof of the lower bound (`lean/lower/`)

`lean/lower/` is a Lean 4 project, derived from Lin's repository (Apache License 2.0, commit 72c022b), that proves

    theorem CenteredMaximal.ofReal_le_weakTypeConstant_two :
        ENNReal.ofReal (17027237335 / 10000000000) ≤ weakTypeConstant 2

and the same statement with the exact rational number above (`ofReal_le_weakTypeConstant_two_exact`). The definitions
of the maximal function and of `weakTypeConstant` are Lin's, unchanged. They use the strict level sets {Mf > λ}, which
give the same constant c₂. The proof extends Lin's proof of his lower bound from his lattice to periodic measures with
finitely many atoms in each cell. Data generated from the certificate (`lean/lower/data/c2lower_1p70272_e2ca1616.json`,
the same file as `lowerbound/c2lower_1p70272.json`) are checked by the Lean kernel (`decide +kernel`), and both
theorems depend only on the axioms `propext`, `Classical.choice` and `Quot.sound`. `lean/lower/README.md` describes
the proof and the build (`lake exe cache get`, then `lake build`, about 25 s once the Mathlib build is in place), and
`lean/lower/PROVENANCE.md` lists the files taken from Lin's repository, unchanged or modified, with their checksums.

## Files that differ from the runs

- The outputs in `logs/controls/` were produced by the original audit files of the second program, whose SHA-256
  digests are recorded in `logs/controls/summary.json`. The test files, `e2e_check.py` and `audit_spec.py` in this
  repository differ from those originals only in their paths, which are now relative to the script and to `testdata/`,
  and in that `audit_spec.py` runs `sys.executable`. `hl_check.py` differs from its original only in that it exits with
  status 1 after a failure.
- The programs in `search/` are those that were run, except for the launch scripts, which were reduced to the two
  relevant commands.
- `lowerbound/NOTE_addendum.md` refers to the files by their archived names with hash suffixes, such as
  `check_exact2_abb311c7.py` for `check_exact2.py`.

## Licence

Apache License 2.0; see `LICENSE`.
