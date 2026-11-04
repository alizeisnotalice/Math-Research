# 第23轮计算接口：三次全局L1界

严格定理和任意概率空间证明在paper/global_l1_cubic.tex。这里只复核有限输入和代数；一般p不由网格候选推断为定理。

## 运行

需要Python3标准库、WolframScript15以及MATLAB R2026a。两个子目录的检查点分开，避免同名job文件互相覆盖。在各自目录执行：

- binary/：`python3 binary_l1_batch.py --workers 15`；`python3 p3_exact_audit.py`。
- finite/：`python3 multi_search.py --workers 15`；`python3 multi_exact_verify.py`；`wolframscript -file p3_multi_symbolic.wl`；MATLAB执行`matlab_multi_coarsegrain.m`。

两套批程序BLAS线程数设1，实际15进程计算；版本和输入签名匹配才复用缓存。不同运行和缓存重读的进程不额外算新证据。目录内生成的完整表、失败清单、检查点及日志可再生，Git仅保存run_summary.json最小统计。

## 二点候选搜索

最终273863个唯一(p,q,a)网格点：p3为114258，p=1.1,1.5,2,4,8各31921；实际原始PID61674–61688。中心化Bregman浮点评估，零失败仅表示程序规则通过；p2最小值0.9999999999999921属数值误差。p3最小候选2.7725424859377363接近已自足证明的C3。其他p只有候选，不作sharp主张。独立Fraction二点程序另核18项恒等式和根夹逼。

## 一般有限向量的最终审计


The final 15-process search ran 155 nonconstant vector rows across dimensions
3, 4, 8, 16, and 64: 60 random rows, 45 embedded critical rows, 45
within-group-spread critical rows, and 5 literal-unit-atom rows.  Its actual
PIDs were `66171` through `66185`, it had zero failures, and its lowest float
candidate was `2.7725424860315533` in an embedded `n=3` row.

Rows were deduplicated by merging equal `f` values and sorting the resulting
`(value,total_probability_mass)` pairs.  By this law-level definition the 155
search rows are **92 unique nonconstant probability laws**.  The final output
is `search_summary.json`; logs, checkpoints, and the empty
`search_failed_inputs.json` are also retained.

The independent Fraction verifier regenerates the same 155 nonconstant rows;
they are an audit, not another 155 inputs to add to the count.  It additionally
runs five all-one constant rows solely through the zero-distance `T=d=0`
branch.  Thus it has 160 rows, **93 unique probability laws** total, and the
same **92 unique nonconstant laws**.  There were zero failures.

For every nonconstant Fraction row it verifies the rational certificate
`A=8T-11d^2 >= 0` and `A^2-125d^4 >= 0`.  Since `d>0`, this is equivalent to
`A >= 5 sqrt(5) d^2`.  It also coarse-grains with the low set `{f<1}` and its
full complement `{f>=1}`; it checks exactly that `T>=T_bin` and that `d` is
preserved.  The complement intentionally includes literal unit atoms.

## Independent local-tool checks

`p3_multi_symbolic.wl` was executed with WolframScript 15.0.0.  It guards the
binary formula, endpoint derivative, critical coefficient, the `d=2q(1-a)`
identity, and the cubic Jensen remainder.  All six guards were true, and any
false guard exits nonzero.  Output: `p3_multi_symbolic_v2.log`.

`matlab_multi_coarsegrain.m` was statically checked and actually run in MATLAB
R2026a Update 4.  It has 150 diagnostic rows: 120 nonzero-distance rows and
30 constant rows; zero failures; minimum positive-distance ratio
`2.772542486031552`.  These MATLAB rows contain repeated fixed controls and
are reported only as rows, not as added unique laws.  Output:
`matlab_multi_summary.json`.

## Limit

Critical samples keep strict positive `a=2^-8,2^-16,2^-32`; the excluded
`a=0` configuration is only a limit.  The package supplies reproducible
finite evidence and exact finite certificates, not a replacement for the
general proof.


## 归档入口再次验证

主代理在归档目录重新执行两组新算、独立有理审计及符号脚本，并检验缓存续跑。归档Wolfram脚本在原6项守卫之外增补一维恒等式，最终7项TrueQ守卫；失败以非零退出。向量长度不是原问题空间维数。旧v1/v2/v3中间日志保留在本机scratch，不混为最终v4统计。
