# 第30轮计算记录：负角幂临界对数端点

## 范围

本记录只复核正文所用的标量局部包络、长类核心项和严格孔质量的
渐近口径，**不**以有限算例证明全输入的极大算子上界或下界。完整
证明、适用参数及速率见 `030_critical_logarithmic_endpoint.md`。

取 \(n=2,3\)，\(\gamma=1/4,1/2,3/4,9/10\)，
\(\beta=1-\gamma\)，\(c=(n-1)/n\)，
\(b=(n-1)\beta/(n-\gamma)\)。对
\(u=2^{-k}\)、\(k=8,16,32,64,128,256,512,1024,4096\)，以
120 位十进制二分独立求解
\[
 L-\log(1+cL)=\log(1/u),\qquad
 F(u)=e^{-\beta L}(1+bL).
\]
归一量为 \(u^\beta\log(1/u)^\gamma\)，比较常数为
\(J=bc^{-\beta}\)；它不同于切片质量式中的原始 cube 质量
\(\kappa=2^n(n-\gamma)/(n\beta^2)\)。长类核心项使用
\(C=\log(1/u)^{1/(2n)}\) 和 \(u_+=uC^n\)；当前严格下界孔使用
\(\delta=1/\log(1/u)\) 和
\(u_-=u(1-2\delta)^n/\delta^n\)。仅当 \(C>2\)、\(\delta<1/3\)、
\(u<1\) 时把行标为几何有效。72 行中 48 行满足此条件；其余行仍
是标量反演诊断，不能作为见证几何的证据。

## 实际运行

从 `src/critical_angular` 实际运行：

```sh
python3 -m py_compile critical_angular_batch.py
python3 critical_angular_batch.py       # fresh: 72 tasks, 15 synchronized workers, 0 failures
python3 critical_angular_batch.py       # replay: 72 signature-matched cache hits, 0 failures
wolframscript -file critical_angular_check.wl
```

MATLAB R2026a 还实际运行了 `critical_angular_independent.m`；代码分析
为零问题，72 行、48 行几何有效，反演残差最大为 \(4.55\times10^{-13}\)。
Wolfram 使用 240 位数值 `ProductLog[-1]`，所有 72 个残差经
`TrueQ[residual < 10^-120]` 守卫通过，并独立符号核对切片体积、质量
及 \(\kappa b=2^n(n-1)/(n\beta)\) 的系数恒等式。Python 120 位结果的最大残差为
\(3\times10^{-117}\)。完整表、运行缓存与日志仅保留在本机 scratch，
未纳入版本库。

在最大尺度 \(k=4096\) 的两例为：

| \((n,\gamma)\) | \(J\) | 局部 | 长类核心 | 严格孔类 |
|---|---:|---:|---:|---:|
| \((2,1/2)\) | .4714045208 | .4723385335 | .4720089023 | .4706869356 |
| \((3,9/10)\) | .09917902324 | .09977776981 | .09965280877 | .09900687444 |

数值显示这些显式量从两侧接近 \(J\)，但不是渐近等式的证明。尤其
表中的长类核心仅为 \(C^{-\tau}F(tC^n)\) 的归一值，尚未包含
\(A_C(1-\Omega_w)\) 分母与短类项，不能被读作完整上界。

## 独立性、修正与局限

Python 使用 15 进程 `Decimal` 单调二分；Wolfram 以数值反支
`ProductLog[-1]` 求解；MATLAB 使用双精度但完全在对数域二分，避免
\(2^{-1024}\) 和更小尺度的下溢。三者复核的是同一明示公式，因而
一致性不能替代公式本身的解析证明。

早期 MATLAB 直接构造 \(u=e^{-H}\) 在 \(k\ge1024\) 下溢，已废弃；
早期 Wolfram 对符号恒等零再数值化会触发额外精度警告，已改成有限
精度反支及 `TrueQ` 残差守卫。早期下界诊断的 \((1-\delta)^n\) 只保留为
Python 输出中的标签化比较；正式计算和正文严格孔均用
\((1-2\delta)^n\)。

可复现源码、参数、摘要及 SHA-256 在
`src/critical_angular/README.md` 与 `round30_minimal_summary.json`。该
程序的固定 \(C\) 只检查主项口径；正文用
\(C=\log^{\gamma/(1+n\gamma)}\) 平衡相对余项，未由这里的有限表声称
其最优性。
