# R004: a nonproportional compact kernel direction

2026-11-05. Status at freeze: the explicit candidate below has passed a new complete C1–C4 Arb run and a separate 192-bit C1/C2/direct-mass verification. Independent review is requested before upgrading the public bound. This is a feasible-kernel result, not a statement of optimality in this direction, in the spline family, or for the sharp weak constant.

## Explicit candidate and exact improvement

Let \(K_0\) denote the frozen scaled kernel from R1, with certificate SHA256
`d570ac3ffc1357e129183eda9579fde8a1bf3775e914efe6ef2d6be336348c99`.
Keep \(\alpha=47/40\), \(R=33/20\), \(h=1/32\), \(N=80\), the positive singular coefficient \(c\), and all other spline classes fixed. For the centered cubic B-spline \(B\), set

\[
P(u,v)=B(32u)B(32v),\qquad K_\varepsilon=K_0-\varepsilon P,
\qquad \varepsilon=1/200.
\]

Thus only `a["0,0"]` is decreased by exactly `1/200`. The candidate is
`computation/round004-kernel-direction/candidate-origin-eps1over200.json`, SHA256
`88af41536e1117d953bbc4594ee99f1319b79b75d98d740c5d3c82caaf1f306f`.
The certificate has 1181 coefficient classes, as before. In particular, this is not another proportional rescaling: the positive singular coefficient is unchanged, while a nonzero compact tensor spline is subtracted.

The standard centered cubic spline has \(\int B=1\). The actual Mathematica 15.0 calculation in `direction_mass.wl`/`.log` verifies this from the two polynomial pieces, and gives

\[
\frac12\int P=\frac{h^2}{2}=\frac1{2048},\qquad
\frac12\int K_0-\frac12\int K_{1/200}=\frac1{409600}.
\]

The separate verifier independently expands the full index set \(-N\le i,j\le N\), rather than reusing the checker's orbit-multiplicity sum. It evaluates

\[
\frac12\int K_{1/200}
=\frac{c\alpha R^{2-\alpha}}{2-\alpha}
 +\frac{h^2}{2}\sum_{i,j}a^{\rm candidate}_{ij}.
\]

Its exact Arb rational endpoints prove

\[
3.590865312932196537
<\frac12\int K_{1/200}
<3.590865312932196538
<\frac{3590865312933}{10^{12}}.
\]

No printed floating-point minimum or decimal ball string is used as an inequality proof. The endpoints are obtained by `lower().fmpq()` and `upper().fmpq()` and compared to exact rationals.

## Reuse audit and choice of direction

The upstream search `search/cp_lp3.py` is a floating-point cutting-plane LP, built on `cp_lp2.py`, `bernK.py`, `bspl_closed.py`, `sing2.py`, and `sing3.py`. It forms a coarse-to-fine spline map, maps both objective and constraints through it, lazily adds Bernstein rows for C1/C2, and adds point rows for the generator. Its pit scan is still finite sampling, not a proof. It requires NumPy/SciPy/HiGHS; local NumPy/SciPy are available but HiGHS is not installed. No large optimizer was installed or run in this round.

The search README records failed intermediate certificates (a narrow C4 dip, then a C1 support-edge failure) and tolerances/margins. The published control `logs/controls/summary.json` records a C4 violation from changing class `31,20` by `1/100000`, and C1/C2 negative controls. These show why coefficient mass removal cannot be accepted from a grid or LP status alone. They were read as interface/guard examples, not reused as proof of this candidate.

The byte-unchanged independent checker is reusable: it accepts rational strings through `Fraction`, computes C1/C2 by exact spline Bernstein nets plus rigorous singular tangents, C3 by direct mass, and C4 by joint polynomial models and remainder bounds. Source hashes are

- `hl_check.py`: `0a5ee85e9b6c250494aaf88abce55bb28a4ba5fc46075c531eeefaa0cc2888f8`;
- `hl_core.py`: `c3d864b1368ad8150fd3366bf77f687a9a3cdd49b924a20d576aa678aed42ebf`.

The original source clone and the R1 portable package were not edited. All new files are confined to the two assigned round004 directories. Imports use the already frozen R1 vendor copy, and Python is run with `-I -B`.

The origin direction has simple support geometry. \(P\) is even in each coordinate, symmetric under exchange, nonnegative, and supported in \([-1/16,1/16]^2\). C1/C2 can change only in the three fundamental knot cells with \(p\le1\). Also

\[
AP(u,v)=h^{-\alpha}\bigl[(L_\alpha B)(u/h)B(v/h)
 +B(u/h)(L_\alpha B)(v/h)\bigr].
\]

Consequently \(AP=0\) when both \(|u|,|v|\ge1/16\); in the fundamental region, only C4 cells with \(q=0,1\) can change. Nevertheless this round numerically recomputed **all** new C4 cells; it does not inherit unchanged cells from the R1 report.

## Actual computation record

`explore.py` evaluated 1294 rational points with Arb at 192 bits, using exact rational coordinates. This was only a screening calculation. The smallest observed point cap for \(\varepsilon\) was about 0.01488994; no inference about untested points or a continuum minimum was made. The single endpoint \(1/200\) was then fixed and retained unchanged during all complete checks.

`make_candidate_and_c123.py` performed the first full C1/C2/C3 check at 192 bits. The new C2 lower endpoint is strictly larger than 1; its diagnostic decimal is about 1.0000000446281452801. C1 covers 1368 knot cells and 2025 leaves; C2 covers 272 cells and 272 leaves.

The new **complete** run used the unmodified upstream checker/core, python-flint 0.9.0, Python 3.12.13, their default 128-bit precision, 10 pool workers, `maxdepth=8`, `kdepth=10`, no `--cells` restriction, and `--only c1,c2,c3,c4`. The checker PID was 5132. Its full report and stdout are `full-c1234-report.json` and `full-c1234.stdout.log`.

| Check | Actual result |
|---|---|
| C1 | PASS; 1368 knot cells, 2025 leaves, zero failures |
| C2 | PASS; 272 cells, 272 leaves, zero failures |
| C3 | PASS; halfmass strictly below `3590865312933/1000000000000` |
| C4 | PASS; 1368 knot cells, 443757 leaves, zero failures |
| Provenance | All C4 workers match the candidate SHA and the main process's recorded source/library provenance |

The wrapper measured 24.94 seconds wall time, 183.51 seconds child user CPU and 1.86 seconds child system CPU. Sampled peak tree RSS was 1384800 KiB and peak summed CPU was 998.3%. The 10-worker setting is not a claim of permanent full utilization. Its 11 observed direct-child PIDs include multiprocessing administrative processes. A misleading original metadata field was renamed accordingly; `runtime-metadata-clarification.json` binds the executed wrapper/raw record and the final normalized metadata, with the one-field reversible change. No numerical record, candidate, checker, or core was changed by that clarification.

`verify_result.py` performs a second cheap verification: exact candidate identity, new report SHA/provenance and coverage, full expected C4 cell set plus result length (excluding omission/duplicate substitution), all zero failures, total leaves, fresh 192-bit C1/C2, independent direct mass, and the point obstruction below. It does not replay C4 leaves. The full run is an actual complete run, but its report is per-knot-cell aggregate output, not a separately checkable per-leaf polynomial certificate. Readers can rerun the byte-unchanged checker if they wish to replay the numerical computation.

## Continuum, distribution, tail, and weak-type transfer

The candidate remains a finite cubic spline perturbation of the same truncated singular power. Hence it is even and integrable: \(\alpha=47/40<2\) gives integrability of the singularity in dimension two, and the spline part is compactly supported. The perturbation is globally \(C^2\) and compactly supported, so its coordinate-jump integral is well-defined and continuous; it introduces no origin atom or interface distribution. The singular coefficient is unchanged.

The complete C1 check supplies global nonnegativity and C2 supplies the diamond majorization. The complete C4 cell set includes all fundamental knot cells where the singular or spline part can be nonzero, including spline regions beyond \(r=R\). The candidate's support is contained in the square \(\max(|u|,|v|)\le9/4\). Wherever the full kernel vanishes, C1 gives the nonnegative jump integral directly, proving C4 in the remaining exterior, including the unbounded coordinate strips of the generator.

The checker proves the pointwise inequality off the axes and the interface \(r=R\). Lemma 3.1 of the supplied Lin–Tomar paper applies to this same finite-spline/truncated-power class: on every compact set away from the origin the difference integral has an integrable dominating function and represents the distribution \(AK\). The null axes and interface therefore do not leave an unverified distributional boundary term. Thus \(AK_{1/200}\ge0\) as a distribution on \(\mathbb R^2\setminus\{0\}\). We do **not** require \(AK\) to be integrable near the origin, or a finite signed measure minus an origin atom; those stronger assumptions are inappropriate at this exponent.

Theorem 2.2 of that paper then gives

\[
C_2\le\frac12\int K_{1/200}<3.590865312933.
\]

This is the weak constant for general \(L^1\) inputs and every positive scale, via the paper's Appendix A transfer and monotone truncation. The linear diamond-to-square transformation has determinant of absolute value 2; it multiplies both input mass and receiving area by 2, so the constant is unchanged. Strict/non-strict level conventions give the same global weak constant by threshold limits. No local Lean build, new Lean theorem, kernel-family optimum, or equality with the sharp \(C_2\) is claimed.

## A precise scope for the one-dimensional direction

Because each admissibility condition is affine in \(K\), and both \(K_0\) and \(K_{1/200}\) are feasible, every \(K_\varepsilon\) with \(0\le\varepsilon\le1/200\) is feasible, with halfmass decreasing exactly by \(\varepsilon/2048\).

The 192-bit control at the exact point \((u,v)=(35/64,5/256)\) proves

\[
AP(u,v)>0,\qquad
A K_{3/200}(u,v)<-0.0000102147.
\]

Both coordinates are positive and \(u+v<R\), so the generator is continuous near this point. Hence the negative value gives a negative open neighborhood, not merely an isolated null-set observation. Since the value decreases with \(\varepsilon\), **every** \(\varepsilon\ge3/200\) on this one affine ray is infeasible. The gap \((1/200,3/200)\) is unresolved; the finite screening does not establish its optimal endpoint. The exact negative upper endpoint and positive direction lower endpoint are retained in `point-control.json` and rechecked in `verify_result.py`.

## Reproduction

From the project root, the already installed small environment is sufficient:

```sh
computation/kernel-baseline/venv/bin/python -I -B computation/round004-kernel-direction/verify_result.py
```

To replace the recorded full C4 run by a fresh one, copy the round004 directory first if preserving its frozen outputs, then run the unchanged checker directly with a new report output path:

```sh
computation/kernel-baseline/venv/bin/python -I -B scratch/kernel-public-package/vendor/hl_check.py \
  computation/round004-kernel-direction/candidate-origin-eps1over200.json \
  --claim 3590865312933/1000000000000 --procs 10 --only c1,c2,c3,c4 \
  --maxdepth 8 --kdepth 10 --out /tmp/round004-kernel-fresh.json
```

For another machine, use Python with `python-flint==0.9.0`; the mathematics does not require matching the old absolute paths or binary hashes. New environment provenance should be recorded separately, while the recorded historical worker digests are checked against their own recorded main provenance.
