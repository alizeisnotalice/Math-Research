# Finite tensor lattice: source-centred radial hardband pressure

2026-10-07. Scope: floating numerical pressure on a huge finite tensor source, not a theorem or a full lower-bound reproduction. Applied `math-l01-lower-bound-test-parameters`; read the source interface at `../../anisotropic_cube_flow/lower_models_interface_20261006.md` (relative to this directory).

## Fixed source and actual query

The complete input is the probability measure

    mu = L^(-n) sum_{z in {0,...,L-1}^n} delta_z,
    L=257, n=8,32,128, total atoms=257^n.

No atom set is expanded at the production dimensions. The exact log10 atom counts are stored as n log10(257) in each result. Every source sample draws all coordinates independently uniformly from the entire finite lattice. Boundary sources are retained. For n=128 the probability of at least one coordinate in {0,256} is about 63%; an interior-only or infinite-periodic experiment would be a different input.

Choose uniformly one of the 2n faces of C=[-1/2,1/2]^n, independently uniform free coordinates, then x=z+r omega with r=exp(t/n), 0<=t<=n log 2. The closed query cube has **side length R**, R in [1,2], and

    A(x,R)=prod_i #{k in {0,...,256}: |x_i-k|<=R/2}/257^n,
    u(x)=max_{1<=R<=2} A(x,R)/R^n.

The intended mathematical winner is the smallest maximizing R. The numerical implementation takes the smaller R only when **computed floating scores are exactly equal**; arbitrary mathematical ties and near score ties are not symbolically certified. Scores are evaluated at R=1,2 and every relevant coordinate integer arrival 2|x_i-k| in (1,2]. Between arrivals the count product is fixed and R^(-n) decreases; thus this is the continuous R-window maximal, not an R-grid approximation. Events with the same floating arrival value are all applied before testing a candidate. Per-coordinate zero counts at R=1 are retained until that coordinate's first arrival. Baseline, first-arrival transitions, ties and endpoint arrivals are explicit.

Production uses local dx=r omega and 2|dx_i-offset|, offset in {-2,-1,0,1,2}; z is used only to test support membership. Since |dx_i|<=1, this covers every possible captured integer. log(u)+n log L is computed directly; no L^n or microscopic lambda is exponentiated.

Fixed thresholds are log lambda=−n log L+c sqrt(n), c=−0.5,0,0.5,1. All c values are reported; none was selected after observing energy. The tested field is

    Y=1_{2 lambda<u<=4 lambda} (r/Rwinner)^n 1_{r<=Rwinner},
    phi_n(t)=E_{source,cone} Y,
    E_n=integral_0^{n log2} phi_n(t)^2 dt.

For normalized cone measure the Jacobian is dx=r^n dt dcone. There is no factor 2^(-n): volume(rC)=r^n and the radial change t=n log r supplies the factor n.

This is the hardband g1 relaxation. Neither original softFIRST, all prior scales, nor the rest of geom qualifications is certified. The finite atomic lattice resembles a basic tensor-lattice pressure source but is **not** A13/A14 with frozen jitter, positive finite spectral compensation, peak widths and synchronized winner packets. Nor is it a B60–B63 reproduction. No complete supplied lower family, positive-width L1 transfer, original weak quotient, or uniform energy order is claimed. Atomic-to-L1 smoothing would require separate fixed-dimension threshold/error analysis; no such transfer is performed here.

## Verification and repaired version

`validation_v2.json`: 360 random n<=4,L<=4 cases and 120 cases with ties, endpoints and lattice boundaries compare the event implementation with explicit all-atom sup-distance sorting. Both scores and winner radii agree within 2e−12. Additional local face-arrival checks use n=8,32,128, z_i=256, r in {1,sqrt(2),2}; the tagged source-face arrival equals r exactly in the local representation. These are implementation cross-checks, not interval certificates.

`internal_formula_check.py` / `internal_formula_validation.json`: 600 generic **interior** cases at n=8,32,128 compare against a separate formula. If U_i=|dx_i|, the second arrival is S_i=2 max(U_i,1−U_i), with forced face S=r; sort these and score k log2−n log S_(k). Winner differences are zero and maximum log-score difference is 1.14e−13. For nonface U uniform on [0,r/2], its CDF is 2(s−1)/r on 1<=s<=r and (s+r−2)/r on r<=s<=2. This check does not condition the production source on interior points.

The old `pressure_v1_superseded.py`, `round1.json`, `round2.json` and their profiles remain unchanged for audit. They are **superseded and invalid as energy evidence**: forming x=z+r omega and then subtracting integer k caused cancellation, converting the mathematically exact tagged source arrival r into r±ulp; `r<=Rwinner` then incorrectly removed some whole samples. Repaired production files are named `v2_round*`. The change materially increased energies.

V2 also handles structural c=0 hardband ties algebraically. If winner R=r, counts are 2^m, and t_j=j n log2/(nodes−1), then log(u L^n)=[m−j n/(nodes−1)]log2. Membership in (log2,log4] is decided by integer cross multiplication. For counts 2^m and winner 1 or 2, the corresponding integer score is used. The same thresholds remain strict below and closed above. Screen records retain how many floating classifications this changed and count any remaining near-threshold values. No numerical tolerance widens a gate. Arbitrary near score ties and continuous values still have floating arithmetic error; the experiment is not an interval proof.

## Frozen three-round protocol and uncertainty

| Round | Samples per independent side | Fine t nodes | Nested coarse nodes | Independent pairs |
|---|---:|---:|---:|---:|
| v2_round1 | 128 | 129 | 65 | 8 |
| v2_round2 | 512 | 257 | 129 | 8 |
| v2_round3 | 1024 | 513 | 257 | 8 |

After these three rounds, a single additional boundary control uses n=128, L=65537, 512 samples per side, 513/257 nodes and eight independent pairs, with the same four normalized c values. `control_preregistration.json` was saved before this run. No L scan is performed. The exact chance of a source having any boundary coordinate is 1−(1−2/L)^n; the saved screens give its empirical counterpart. The L change is a separate source and is not a fourth point in a fitted dimension law. The configurable L/dimension arguments added after round2 do not change the main estimator and passed rounds are not overwritten.

Each pair has two independent source/cone draws A and B, and the eight replicate pairs within each run use distinct fixed seeds. Within a side the same sampled paths are retained across t. The estimator integrates phi_A phi_B, so at each fixed node it is unbiased for phi^2. Same-sample phi_A^2 is saved only as a positive-bias diagnostic. Seeds and repeats are fixed in the script; `v2_round*_n*_rep*.npz` saves t, phi_A, phi_B. **Runs are not all mutually independent:** the tag suffix 1 makes the L65537 control share the seed namespace with v2_round1, so these different-L/different-sample runs may be coupled. The final main v2_round3 has suffix 3 and different seeds from the control; its comparison with the control uses different namespaces. Old/new repair comparisons intentionally use the same seeds. Final and partial JSON snapshots preserve all dimensions, c values, estimates, screens and timings. Old rounds are not overwritten.

Approximate marginal 95% intervals use the eight independent replicate energies with t_7=2.3646; eight replicates do not guarantee a normal/t approximation, and a zero interval from eight zero observations is **not** proof of zero energy. A separate bounded Hoeffding 95% interval is provided for each grid statistic Z in [0,n log2], using that range and eight independent pairs. These rigorous probabilistic grid bounds are generally very wide and are not simultaneous over all entries.

Coarse/fine trapezoid differences are computed from the same paths and saved per replicate. The hardband and winner gates are discontinuous, so these differences are resolution diagnostics, not deterministic quadrature error bounds. MC intervals cover the finite-node numerical target, not the continuum energy. Near-window-endpoint and winner-equals-r screens include intentionally exact cone-face/end-node events. Symbolic c=0 tie resolution does not constitute an interval certificate for arbitrary floating comparisons. No dimension-uniform envelope follows from these data.

Execution command (bundled Python with NumPy):

    /Users/zhengzhihao/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 pressure.py --validate --round 1

Run `--round 2` and `--round 3` separately. Existing round tags refuse overwrite. Final three-round numerical tables and completion receipts are in `summary.md` and `manifest.json` after all runs finish.
