# Astra lower001: off-lattice satellite mass insertion

The registered batches are complete. They do not improve the standing exact
constant C from `../lower_round003/a2_equals_a3_results.json`. No record lock,
formal ledger edit, commit or push was made. The main clone was read only.

## Mechanism and finite scope

The new input adds a second periodic orbit of small mass ε to the baseline
rectangular unit-atom lattice. Its half-period displacement is selected among
the nine symmetry types `(first coordinate shifted or not, k shifted tail
coordinates)`, k=0,…,4, excluding the zero displacement. This includes mixed
captures of the two orbits at one common radius, not just the union of the
two separate halos. It differs from the previous equal-mass five-axis BCC
batch and from perturbing masses at existing lattice sites.

The first registered batch used ε=1/1024 and1/512 at fixed spacings, giving
18 floating diagnostics. All are below C. The closest is the orbit shifted
by half a period along one tail coordinate only: at ε=1/1024 the ratio is
approximately1.9940413966874526, a gap of0.0002750648810332.

Since the original spacings lie at capture kinks, inserting mass also changes
the relative scale even when two orbits merge. A separately registered
followup multiplies every spacing by(1+ε)^(1/5), preserving period density.
Its nine candidates at ε=1/1024 also all fall below C. The same one-tail
shift is closest, approximately1.9940671048365064, a gap of0.000249356732.
The merged-orbit scaling calibration was actually executed and agrees with
C's decimal display within2e-15; the exact invariance follows from dilation.
The finite differences saved in the diagnostic outputs are not certified
derivatives and do not exclude arbitrary displacements or other masses.

## Exact reconstruction of the nearest two cases

A third registration selects exactly the closest fixed-spacing and codilated
cases for rational reconstruction, to verify geometry and normalization.
Each spacing and every equality side length is rounded strictly inward to
denominator2^64 using integer fifth-root isolation. Thus the final inputs
and radii are rational and have strict level-one margins.

| Case | Exact-rational ratio, decimal display | Conclusion |
|---|---:|---|
| Fixed baseline spacing approximation |1.9940413966874522025005485898832207|Strictly below C's lower isolating endpoint|
| Density-preserving codilation approximation |1.9940671048365059233907734135104615|Strictly below C's lower isolating endpoint|

Both runs process4608 configurations, generate135 distinct boxes and retain
63 boxes. Recursive slicing and independent coordinate-mask integration
agree exactly. The respective slicing counts are3279 and3074; actual Python
times are approximately0.52s and0.48s. Every retained witness checks distinct
atoms, exact mass and the single strict radius inequality. The exact ratios,
every retained receiver box, every mass and radius, comparison intervals and
reference/source hashes are saved in the two exact result files.

`verify_saved.py` does not import either producer. It reconstructs the
selected coordinates by nearest-distance sorting, checks92,630 atom/receiver
endpoint inequalities in each case, recomputes the rational union by a
frozenset slab method, and reproduces both exact ratios. This is an
additional same-agent implementation check, not an independent agent audit.

`transfer.md` supplies the complete finite collar and smoothing argument:
period mass1025/1024; reflected receiver volume32V; collar mass
(1025/1024)(M+7)^5; received volume32VM^5; small smooth bumps at an enlarged
common radius that retains the exact strict margin. Consequently each
recorded A is a legitimate all-positive-scale L¹ lower bound, but is weaker
than the standing C. No new global exclusion or operator upper bound follows.

## What remains open

Half-period satellites have an unfavorable small-mass response at the
registered scales. Arbitrary off-lattice displacements, genuinely correlated
multi-orbit patterns, or a macroscopic multiscale construction remain open.
An exposed-boundary influence kernel could locate profitable generic
insertions before new search, but it was not implemented here. The initial
macro-density idea was rejected only heuristically because of overlapping
coverage; no theorem excluding multiscale constructions is asserted.

Reproduction, from the clone root:

```sh
python3 -B cube-dim-5-bounds/astra_lower001/satellite_screen.py
python3 -B cube-dim-5-bounds/astra_lower001/codilation_screen.py
python3 -B cube-dim-5-bounds/astra_lower001/certify_satellite.py fixed
python3 -B cube-dim-5-bounds/astra_lower001/certify_satellite.py codilation
python3 -B cube-dim-5-bounds/astra_lower001/verify_saved.py
```

NumPy is needed only for the floating screens; the exact paths use Python's
standard library. This round used local Python only and did not touch the
shared MATLAB instance. All retained files are necessary registrations,
final source, compact result evidence, the transfer proof or this report.
