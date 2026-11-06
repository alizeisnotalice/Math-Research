# Fan--Grama--Liu：条件 Bernstein 鞅的相对正态尾核验

日期：2026-11-06（Asia/Shanghai）。本记录直接核验 Fan--Grama--Liu 的 arXiv v3 全文及定理的证明段落；不把普通 Berry--Esseen 的绝对误差或相邻的指数不等式替作相对尾展开。

## 已核验的来源与定理

**已核验版本。** Xiequan Fan, Ion Grama, Quansheng Liu, *Cramér large deviation expansions for martingales under Bernstein's condition*, arXiv:1210.2198v3（14 Sep 2014，26页）；期刊版 *Stochastic Processes and their Applications* **123** (2013), 3919--3942，DOI [10.1016/j.spa.2013.06.010](https://doi.org/10.1016/j.spa.2013.06.010)。

原文第 2 节设
\[
X_k=\sum_{i=1}^k\xi_i,\qquad
\langle X\rangle_k=\sum_{i=1}^k\mathbb E(\xi_i^2\mid\mathcal F_{i-1}),
\]
其中 \((\xi_i,\mathcal F_i)\) 是平方可积鞅差阵列。条件为：

\[
\tag{A1}
\left|\mathbb E(\xi_i^k\mid\mathcal F_{i-1})\right|
\le {k!\over2}\epsilon^{k-2}
\mathbb E(\xi_i^2\mid\mathcal F_{i-1}),
\quad k\ge3,
\]
和
\[
\tag{A2}
\left|\langle X\rangle_n-1\right|\le\delta^2\quad\hbox{a.s.},
\qquad 0<\epsilon\le\tfrac12,\quad0\le\delta\le\tfrac12.
\]

Theorem 2.1（原文 (8)）对每个固定 \(\alpha\in(0,1)\)、所有 \(0\le x\le\alpha\epsilon^{-1}\)，给出右尾上界
\[
{\mathbb P(X_n>x)\over\overline\Phi(x)}
\le e^{c_\alpha(x^3\epsilon+x^2\delta^2)}
\bigl[1+c_\alpha(1+x)(\epsilon|\log\epsilon|+\delta)\bigr].
\]
Theorem 2.2（原文 (10)）给出存在绝对 \(\alpha_0>0\)，当 \(0\le x\le\alpha_0\epsilon^{-1}\)、\(\delta\le\alpha_0\) 时，
\[
{\mathbb P(X_n>x)\over\overline\Phi(x)}
\ge \exp\{-c_{\alpha_0}[x^3\epsilon+x^2\delta^2
+(1+x)(\epsilon|\log\epsilon|+\delta)]\}.
\]
于是原文紧接的式 (14) 是对数相对误差的双边界。定理证明的第 5--6 节使用指数倾斜；例如式 (54)--(56) 先得到倾斜尾下界，再将 \(\overline\Phi(\underline\lambda)\) 与 \(\overline\Phi(x)\) 比较。这确认定理确实是相对 Gaussian tail，而非绝对 Kolmogorov 界。

## 对所给尺度的直接代入

设条件 (A1) 以 \(\epsilon_n=C_\epsilon n^{-1/2}\) 成立，并且给定的平方变差条件精确地是
\[
|\langle X\rangle_n-1|
\le \delta_n^2=C_V n^{-1/4}(\log n)^{1/4}\quad\hbox{a.s.}
\]
（这里 \(C_V\) 是固定常数）。则定理中使用的参数为
\[
\delta_n=\sqrt{C_V}\,n^{-1/8}(\log n)^{1/8},
\]
而不是右端 \(C_Vn^{-1/4}(\log n)^{1/4}\) 本身。

对全部 \(0\le x\le2\sqrt{\log n}\)，Theorems 2.1--2.2 的共同误差项一致满足
\[
\begin{aligned}
x^3\epsilon_n
 &=O\!\left(C_\epsilon n^{-1/2}(\log n)^{3/2}\right),\\
x^2\delta_n^2
 &=O\!\left(C_Vn^{-1/4}(\log n)^{5/4}\right),\\
(1+x)(\epsilon_n|\log\epsilon_n|+\delta_n)
 &=O\!\left(C_\epsilon n^{-1/2}(\log n)^{3/2}
 +\sqrt{C_V}\,n^{-1/8}(\log n)^{5/8}\right)=o(1).
\end{aligned}
\]
同时 \(2\sqrt{\log n}\le\alpha_0\epsilon_n^{-1}\) 且 \(\delta_n\le\alpha_0\) 对所有充分大的 \(n\) 成立。因此，**如果 A1、A2 和鞅差条件已逐项严格证明，则**
\[
\sup_{0\le x\le2\sqrt{\log n}}
\left|{\mathbb P(X_n>x)\over\overline\Phi(x)}-1\right|\longrightarrow0.
\]
这是可用于该窗口的统一相对尾结论。原文没有给出本项目模型的 A1/A2；那些必须独立建立。

## 不可省略的适用条件

1. **可预测二次变差。** A2 约束的是 \(\langle X\rangle_n=\sum_i\mathbb E(\xi_i^2\mid\mathcal F_{i-1})\) 的几乎处处偏差。仅有 \(\sum_i\mathbb E\xi_i^2=1\) 或 \(\operatorname{Var}(X_n)=1\) 不足以套用定理。
2. **条件矩的形式。** 原文 A1 是 \(\left|\mathbb E(\xi_i^k\mid\mathcal F_{i-1})\right|\) 的界，即条件**有符号幂矩的绝对值**；它不是 \(\mathbb E(|\xi_i|^k\mid\mathcal F_{i-1})\) 的陈述。后一个以同一右端成立时是更强的充分条件，因为 Jensen/三角不等式给出 A1；无条件矩界则不是。
3. **归一化。** 定理中的总和已经按单位可预测方差归一化。若实际变量为 \(\eta_i\) 且总方差约 \(n\)，必须先取 \(\xi_i=\eta_i/\sqrt n\)，并相应把条件 Bernstein 常数转成 \(\epsilon_n=H/\sqrt n\)。
4. **不由标题推断。** 本结论来自 v3 的 Theorems 2.1--2.2 和式 (14)，不是 2015 年同作者的条件 MGF/最大不等式论文。

## 结论边界

在所给尺度下，原文定理的渐近条件足以推出 \(x\le2\sqrt{\log n}\) 的相对误差 \(o(1)\)，但不能单凭本条文献完成任何“全体输入”的第三项：还须在所构造的停时/低均值控制中证明鞅差、A1 的同一 \(C_\epsilon\)，以及 a.s. 的 A2。并且 \(\alpha_0,c_{\alpha_0}\) 是存在性常数，论文没有给出可用于有限 \(n\) 的数值阈值。MATLAB R2026a 对 \(C_\epsilon=C_V=1\) 计算的理论误差包络在 \(n=10^4,10^6,10^8\) 分别约为 11.93、5.92、2.98，说明其趋零很慢，不能把渐近判定冒充有限维认证。

## 来源定位

- [arXiv 摘要和版本/期刊信息](https://arxiv.org/abs/1210.2198)：v3、作者、26页、期刊卷页和 DOI。
- [arXiv HTML 原文](https://arxiv.org/html/1210.2198v3)：第 2 节的 (A1)、(A2)、Theorems 2.1--2.2、Corollary 2.1 和式 (14)；第 5--6 节的 (54)--(56) 显示相对尾比较的证明链。

```bibtex
@article{FanGramaLiu2013,
  author  = {Fan, Xiequan and Grama, Ion and Liu, Quansheng},
  title   = {Cram{\'e}r Large Deviation Expansions for Martingales under Bernstein's Condition},
  journal = {Stochastic Processes and their Applications},
  volume  = {123},
  pages   = {3919--3942},
  year    = {2013},
  doi     = {10.1016/j.spa.2013.06.010},
  eprint  = {1210.2198},
  archivePrefix = {arXiv},
  primaryClass = {math.PR}
}
```
