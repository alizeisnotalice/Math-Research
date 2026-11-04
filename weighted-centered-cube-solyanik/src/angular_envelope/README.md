# 第24轮：精确角向稳定包络的计算接口

完整全p>1证明在paper/angular_envelope.tex。Python只需标准库；WolframScript15和MATLAB R2026a分别运行其脚本。精确有限证书与非整数数值诊断分开报告。

## 在本目录复现

- `python3 decimal_solver.py`：独立90位候选求解器。
- `python3 fraction_interval_verify.py --workers 15`：有理根夹逼与H区间。
- `python3 verify_certificates.py`：独立重放最小四点证书；附参数fraction_summary.json可重放103项全量证书。
- `python3 direct_exact_families.py`和`python3 integer_p_fraction_diag.py`：有限精确构造。
- `wolframscript -file cubic_envelope_guards.wl`与`wolframscript -file general_convexity_guards.wl`：10+4项精确符号守卫。
- `wolframscript -file general_p_highprecision.wl`：45项260位诊断；检查内支残差、外支反演与行数，失败非零退出。
- MATLAB执行matlab_envelope_independent.m与matlab_general_p_diagnostic.m。

# Round 24 fixed-d cubic-envelope certificates

This scratch package checks the explicit envelope in `round24.md`.  It does
not replace its one-variable proof.  The fixed convention is
`0<d<2`, `D=E f^3-1`, and `H(d)` is the stated inner-root/outer-boundary
function.

## Fraction interval verifier

`fraction_interval_verify.py` is the final 15-process verifier.  It covered
103 rational values of `d`: 37 interval inner-branch values, one exact-root
inner value, and 65 outer-branch values, including dyadic approach to zero,
the switch, and two.  The archive fresh run used 15 actual workers: 72670, 72671, 72672, 72673, 72674, 72675, 72676, 72677, 72678, 72679, 72680, 72681, 72682, 72683, 72684. All 103 certificates passed.

For every non-exact inner `d`, it uses only Fraction arithmetic.  It brackets the root
in `(1/3,1/2)` by `h_d(left)<0<h_d(right)`, where
`h_d=r(q)*(d-mathfrak_d(q))`, and then uses positive rational interval
arithmetic to enclose `H(d)=G_d(q_d)`.  The largest final root interval is
below `2^-98`; the largest final H interval is below `2^-95`.  The added
`d=mathfrak_d(3/8)` hits its rational root exactly and is stored as an exact
`H` value rather than a degenerate strict-sign interval.  For `d>=2/3`,
it records the exact rational outer expression `d*(4-d)/(2-d)^2`; in
particular the switch has exact value `H(2/3)=5/4`.

The first orchestration attempt incorrectly put second-wave tasks behind a
15-worker barrier and ended with 11 `BrokenBarrierError` records.  No
mathematical certificate failed.  The final runner synchronizes only its first
15 tasks.  The later exact-root revision reran all final certificates under a
new signature and writes `fraction_summary.json` with zero failures.  Both
logs are retained; `fraction_interval_v3.log` is the final p=3 run.

## Exact direct families and numerical candidates

`direct_exact_families.py` passed 29 exact rows:

* 7 rational `q in (1/3,1/2)` parameter rows with `d=mathfrak_d(q)`, exact
  stationary identity, strictly positive binary values, and `D=H_param`.
* 15 strict-positive outer rows, `q=d/2+2^-k`, checking fixed-`d` approach
  from above to the outer infimum.
* 5 fixed one-dimensional equal-atom checks, including `d=0`, verifying
  `D=3d^2` for `0<d<1`; one literal unit-atom row; and one explicitly labelled
  nonnegative zero-value boundary audit excluded from the strict-positive
  class.

`decimal_solver.py` is a separate 90-decimal numerical candidate solver on
nine representative points.  It is not used as a certificate.

## Mathematica and MATLAB

`cubic_envelope_guards.wl` was actually run in WolframScript 15.0.0.  Ten
`TrueQ` guards cover the binary formula, strict convexity formula, boundary
derivatives, `partial_d G`, stationarity/parameter identities, outer formula,
and switch value.  Any failed guard exits nonzero; all passed.  See
`cubic_envelope_guards.log`.

`matlab_envelope_independent.m` was statically checked and actually run in
MATLAB R2026a Update 4.  It independently bisects the root and inverts the
parametric/outer map on 98 diagnostic rows.  It reports zero failures, maximum
root residual `4.996003610813204e-16`, and maximum inverse error of the same
size.  This is floating diagnostic evidence only.

## General-p numerical diagnostic

The general-p extension is separate from the p=3 certificates.
`integer_p_fraction_diag.py` supplies 20 exact rational derivative and
strict-convexity rows for `p=2,3,4,8,16`; it is not a noninteger root
certificate.  `general_p_highprecision.wl` was run after exact parameter
construction with 260-digit arithmetic and 1,000 bisection steps.  It covers
`p=1.1,1.5,1.9,2,2.1,3,4,8,16`, five points each including `d=2/p` and a
`10^-12` inner-side point.  All 45 diagnostic rows passed; the largest inner
derivative residual was about `1.91e-259`.

`matlab_general_p_diagnostic.m` independently ran 45 double-precision rows.
It attempts every inner solve first, resolves 25 ordinary inner roots and all
outer inverse checks with zero failures (largest resolved residual
`1.78e-15`, outer inverse error `2.22e-16`).  Only two subquadratic
near-switch rows are labelled `indeterminate_near_switch` after their actual
computed `q-d/2` falls below `128` machine epsilons; they are not presented as
resolved strict roots.  In particular the `p=1.9` near-switch row is actually
resolved.  The high-precision Wolfram run covers the two double-degenerate
rows.

## Limit

The scripts distinguish direct exact constructions, rational interval
certificates, and floating candidate checks.  They do not use a numerical
grid to prove global minimality, uniqueness, or sharpness; those claims remain
those of the analytic argument in `round24.md`.

## 主代理归档裁决

原scratch的102个有理d未显式包含d=2/3；归档版本v3已增补该点并实际新算103项。构成37区间内支、1个精确有理根、65外支。原始102输出不被当作已覆盖此分界的证据。独立重放器不导入生成器：以临界三次式的有理符号夹根，再检验四次多项式的Bernstein系数全部非负，从而证明整个根区间上的H上下界。全103项和最小四点均可单独重放。

20项整数p有限诊断只核查所列点的公式正性。29行构造与其他参数集可能重合，不相加作新输入。MATLAB一般p的45行最终为25已解内支、2未定近切换点、18外支；未定行不是严格根通过。Wolfram在exact参数生成后用260位覆盖对应退化。早期80位精度不足、barrier失配和保守跳过的记录保留本机，未算入最终通过。

完整表、失败清单及检查点可再生，不入Git；最小证书与run_summary.json保留核验接口。
