# Round 21 computation archive

All work here uses the outer Lebesgue halo only. No script averages with
respect to the weight when constructing a halo.

## Exact certificates

`exact_main_certificate.py` uses a side-effect-free copy of the final
round-20 Fraction endpoint/cell routine and verifies
\(H_{1/4}(0,1)=(-1,2)\), modulo endpoints.  For
\(b=1+64\mathbf1_{(1,2)}\), \(v=|x|^{-1/2}\), it records
\(w(H)/w(E)=(69+2\sqrt2)/3\), the strict square gap
\(45^2-2\cdot22^2=1057\), and the endpoint residual
\(4128575/8388608>0\) at \(t=2^{-16}\).  The mixed expression is called a
lower bound and the pure expression an upper bound, not an equality.  It also
checks the witness-only condition \(B>19+22\sqrt2\): integer \(B=51\)
passes and \(B=50\) fails because \(31^2<968<32^2\).  This does not claim a
sharp global threshold.

`scaled_bump_certificate.py` is the exact BG11 family, not the earlier
discarded fixed-height scaling: \(y=2^{-k}\),
\(w_\epsilon=1+|x|^{-1/2}+(64/y)\mathbf1_{(y^2,2y^2)}\),
\(E=(0,y^2)\), \(0\le k\le60\).  It verifies bump \(L^1\)-mass \(64y\)
and the witness formula \((3y+66+2\sqrt2)/(y+2)\), which is at least
\((69+2\sqrt2)/3\).  The difference has numerator
\((60+2\sqrt2)(1-y)\); strict base comparison uses the square gap above.
No optimality claim is made.

`family_windows.py` makes 72 exact dyadic checks (\(n=1,2,3\), six rational
\(\tau\)'s, and \(B/a=4,16,64,130\)) of the genuine sufficient comparison
\[
 \tfrac12t^{\tau/n}>2^n(B/a)\frac{t}{1-t},\qquad t^{1/n}\le1/8.
\]
It gives a certified \(t_0=2^{-Lk_0}\), with \(1\le k_0\le10\), and the
reason it remains valid for every \(0<t\le t_0\).  The output separately
records optional second-order predicates without using them to exclude a
first-order window.  It is exact algebra conditional on the analytic
geometric lower and background upper estimates, not a proof of those
estimates.

## Python batch

`outer_batch.py --workers 15` generated `outer_batch_v3.log`,
`outer_summary.json`, and 15 checkpoint records. Python was 3.15.0 on macOS
27.2. BLAS thread environment variables are set before imports. The exact
Fraction geometry uses 15 actual worker PIDs:

`57162, 57163, 57164, 57165, 57166, 57167, 57168, 57169, 57170, 57171,
57172, 57173, 57174, 57175, 57176`.

The final run was freshly computed under
`outer-homogeneous-v3-bg5-cacheguard-20261105`: each checkpoint includes a
SHA-256 encoding of its exact Fraction configuration, so reuse requires both
the signature and version. The summary labels computed/reused state and both
current and checkpoint PID.

The 15 configurations have global backgrounds
\(b=1+\sum h_i\mathbf1_{I_i}\), hence \(1\le b\le B\), and two separated
singularities. All 75 fixed \(E,\alpha\) rows passed common-halo convexity;
the largest `mixture - max(component ratios)` was `0.0`. All 75 also passed
the direct BG5 check \(b(H)/b(E)-1\le2B(1-\alpha)/\alpha\); largest residual
was `-0.08080808080808081`. There are 268 local-hole rows at both sites,
\(\epsilon=2^{-k}\), \(4\le k\le12\). Largest absolute local difference
was `0.76884` at the first tested scale and `0.04467` at the last. These are
numerical consistency data only, not a proof of convergence or eventual
equality of global constants. The program exits nonzero on validation failure.

## Independent MATLAB check

`matlab_outer_independent.m` was statically analyzed and actually run with
MATLAB R2026a Update 4 (`26.1.0.3312084`). It checked 60 finite local-hole
examples using a separately written floating endpoint/cell halo routine,
primitive singular integrals, and base-plus-step background integrals.
`matlab_outer_summary.json` reports zero failures, maximum convex slack
`-0.003325050403336416`, and final-scale local difference
`0.027897731080297028`. The first execution exposed a one-sided-cell row
concatenation error, which was corrected, statically rechecked, and then
rerun successfully. MATLAB evidence is numerical, not an interval proof.

## Limits

No calculation here proves eventual equality of global Solyanik constants.
The fixed \(\alpha=1/4\) contrast specifically rules out an all-\(\alpha\)
extrapolation. The scripts are reproducible diagnostics; analytic local
domination and finite-separation arguments remain proof obligations.

## 归档运行命令

`python3 exact_main_certificate.py`、`python3 scaled_bump_certificate.py`、`python3 family_windows.py`只需标准库；`python3 outer_batch.py --workers 15`重建全表与检查点。Mathematica15运行`wolframscript -file symbolic_check.wl`；MATLAB运行`matlab_outer_independent.m`，不清空用户会话。`run_summary.json`为原运行最小摘要；完整生成表与日志不入Git。

归档版本另把BG11比值的两个有理系数由实际成本/光晕质量逐项重算，再比对其符号公式，防止只打印目标公式；父代理已实际重新运行。原始试探中未接回H9的泛化下界、错误的固定高度缩放族、MATLAB拼接失败均已纠正，不列为有效证据。

归档目录实际新算15任务、记录15个不同worker，随后续跑复用15任务且零新算；缓存续跑不计为新独立证据。二阶审计36个(n,τ,q)组合中24个满足当前严格范围，未通过二阶条件不影响首项窗口。
