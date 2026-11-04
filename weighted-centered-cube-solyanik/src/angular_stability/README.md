# Round 22: finite-probability \(L^p\) stability archive

This package concerns only positive finite probability vectors \((r_i,f_i)\)
with \(\sum r_i=1\), \(\sum r_if_i=1\).  It contains no outer-halo,
maximal-average, or weighted-geometric calculation.

## Python 15-process diagnostic

`lp_stability_batch.py --workers 15` was run under Python 3.15.0 on macOS
27.2.  BLAS-related environment variables are set to one before imports. The
final version is `lp-angular-v2-sharp-spike-scale10-20261105`; each checkpoint
requires both this version and a SHA-256 signature of the deterministic input
family before it can be reused.  `lp_batch_v2.log`, `lp_summary.json`, and 15
checkpoint files preserve the run.  The actual computed worker PIDs were

`60224, 60225, 60226, 60227, 60228, 60229, 60230, 60231, 60232, 60233,
60234, 60235, 60236, 60237, 60238`.

All 15 groups were freshly computed, producing 1,500 rows: the ten requested
values \(p=1.1,1.25,1.5,1.9,2,2.1,3,4,8,16\), three balanced vectors, three
hole vectors, three spikes, and one deterministic positive random angular
vector per \((\text{seed},p)\).  Every input is rational and is normalized
to exact weighted mean one before floating evaluation.  `failed_inputs.json`
is the empty failed-input list; the program exits nonzero if any validation
fails.

For float calculations, \(E f^p-1\) is evaluated as the centered Bregman
sum
\[
 E\{f^p-1-p(f-1)\},
\]
and the \(1<p\le2\) right side uses `log1p`/`expm1`.  This avoids direct
subtraction of nearby values when \(S\) is close to one, but it remains
floating numerical evidence.  All requested diagnostics passed:

* \(1<p\le2\): \((E|f-1|^p)^{2/p}\le(S^2-1)/(p-1)\).
* \(p\ge2\): \(E|f-1|^p\le Ef^p-1\) and
  \(E|f-1|^2\le(Ef^p-1)/(p-1)\).

At the balanced \(\delta=2^{-16}\) rows, the normalized low-\(p\) sharp
ratio ranged from `0.9999999998959525` to `0.999999999992724`.  At the final
spike scale \(M=2^{10}\), the \(D_p^p/(Ef^p-1)\) ratios for
\(p=3,4,8,16\) were `0.9970760`, `0.9960995`, `0.9922142`, and `0.9844889`.
For the \(q=2^{-9}\) hole vectors, the normalized variance ratios were
`0.9980507`, `0.9967516`, `0.9924932`, and `0.9846029`.  These values are
consistent with, but do not prove, the stated sharp-sequence limits.

## Independent exact integer audit

`exact_integer_audit.py` independently reimplements the inputs and verifies
Fraction arithmetic for \(p=2,3,4,8,16\).  It checked 120 rows, including
scalar Bregman inequalities and their averaged forms, with zero failures;
`exact_integer_audit.json` records the exact ratios of spike and hole
sequences.  The main batch separately contains 750 exact-integer rows.  Thus
the exact audit is corroboration rather than reliance on a floating tolerance.

## Independent MATLAB diagnostic

`matlab_lp_independent.m` was statically analyzed and actually executed in
MATLAB R2026a Update 4 (`26.1.0.3312084`).  It uses an independently written
centered-Bregman float evaluation on 200 vector cases (20 seeds times all ten
\(p\)'s, rotating balanced/hole/spike/random families).  It reported zero
failures and maximum violation `1.39e-17`; data are in
`matlab_lp_summary.json`.  The code does not clear the user workspace or reset
MATLAB paths.  MATLAB results are finite floating diagnostics, not interval
certificates.

## Scope

Finite vectors cannot prove the theorem for all probability weights.  The
rigorous angular-stability proof, including any limiting and sharpness claims,
must be checked independently.  This archive only gives reproducible rational
test inputs, exact integer cases, stable floating checks for noninteger \(p\),
and explicit numerical sharp-sequence probes.

## 复现与证明接口

主批`python3 lp_stability_batch.py --workers 15`；独立整数审计`python3 exact_integer_audit.py`；符号复算`wolframscript -file scalar_checks.wl`；MATLAB运行`matlab_lp_independent.m`。原始完整表及失败清单可再生，不入Git；`run_summary.json`保留最小统计。

正文paper/angular_stability.tex的全域不等式和sharp局部L1常数均完整证明。有限概率向量可在高维无原子锥测度中实现为简单角向函数，但不能当作固定一维两原子角向测度的完整穷举。浮点容差通过不等于严格证书；严格Fraction结论仅限已列整数p与有限数据。经典Lp工具不作本项目新发现。

去重口径：1500行中含跨15种子重复的固定控制；按(p,有理概率律)合并后共240种不同测试，其中整数p为120种（主表750行）。不把重复控制当作独立证据。独立120项审计与主表相应唯一整数数据范围一致，代码另写，未共享其计算函数。
