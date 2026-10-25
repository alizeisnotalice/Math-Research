# 第二批：审计 (2^5) 的 Besicovitch 覆盖证明

检索/阅读日期：2026-10-25。问题是：已登记的公开上界
\(C_5\leq 2^5=32\) 的**原始证明**在 Lebesgue 五维方体情形能否留下一个可量化的正余量 \(\delta\)，从而直接得到 \(C_5\leq32-\delta\)。本批没有重新筛选带隐常数的 \(O(n\log n)\) 工作；新增定向网络检索数为 0（限制为至多 4），新增原始来源数为 0。它重新取得并阅读了已登记的 Aldaz 原始预印本 v2 的关键证明段。

## 原始来源、版本及准确常数

J. M. Aldaz, [*Kissing numbers and the centered maximal operator*, arXiv:1811.10478v2 (10 Mar 2020)](https://arxiv.org/pdf/1811.10478v2) 的实际下载 PDF（20 页）SHA-256 是

```
d42a1d6518298227e2af47876761d47225e38b28e2f7726f91c0cc2498ab12e7
```

已阅读位置与所给内容如下。

1. Definition 1.1（p.2）定义 Besicovitch family：不同球互不含对方的中心，并定义重叠常数 \(L(X,d)\)。
2. Theorem 2.5、Lemma 2.6（pp.5--6）：若 \(L(X,d)\le L\)，贪心选择的 balls 覆盖 level set，且 \(\sum_{B\in\mathcal C'}1_B\le L\)；其式 (2) 的最后一步为
   \[
     \mu\{M_\mu f>t\}\le {1\over t}\int\!\Big(\sum_{B\in\mathcal C'}1_B\Big)f\,d\mu
     \le {L\over t}\int f\,d\mu.
   \]
   Theorem 1.2 及 Theorem 2.10 说明该 uniform constant 是对所有允许背景测度的上确界。
3. Theorem 3.2（pp.10--11）给出 \(L(\mathbb R^d,\|\cdot\|)=H^*(d,\|\cdot\|)\)；Corollary 3.5（p.11）对 \(\ell_\infty^d\) 明确给出 \(H^*=2^d\)，因此 general-measure uniform 常数恰为 \(2^d\)。代入 \(d=5\) 正是 32。

这解释了现有 \(C_5\le32\) 的来源：Lebesgue 只是该 uniform 定理允许的一个背景。Theorem 3.3 的 sharpness 构造则使用离散背景测度；它**不能**反推 Lebesgue 的等号。

## 五维的准确“无余量”几何核验

不能通过“32 个 cube 只会在边界处相交，Lebesgue 测度为零”来从上述证明扣掉任何 \(\delta\)。令 \(\varepsilon\in\{-1,1\}^5\)，取 32 个闭 \(\ell_\infty\) 球（共同半径 1）

\[
 Q_\varepsilon=B_\infty\!\left({3\over4}\varepsilon,1\right).
\]

它们全部包含 \(0\)；若 \(\varepsilon\ne\eta\)，则在某坐标符号不同，故

\[
 \left\|{3\over4}\varepsilon-{3\over4}\eta\right\|_\infty={3\over2}>1.
\]

所以任何 \(Q_\varepsilon\) 的中心都不在另一个中，它们构成 Definition 1.1 的 Besicovitch family。更强地，

\[
 \bigcap_{\varepsilon\in\{-1,1\}^5}Q_\varepsilon
   =[-1/4,1/4]^5,
 \qquad
 \lambda^5\!\left((-1/4,1/4)^5\right)=(1/2)^5={1\over32}>0.
\]

因此重叠函数在一个正 Lebesgue 测度集上恰取 32，而非仅在一个边界点上取 32。最小复核程序 [verify_strict_upper_overlap_cube5.py](verify_strict_upper_overlap_cube5.py) 使用 `Fraction` 运行，生成 [strict_upper_overlap_cube5_certificate.json](strict_upper_overlap_cube5_certificate.json)。它实际验证 32 个 cubes、最小中心 \(\ell_\infty\) 距离 \(3/2\)、以及共同开交体积 \(1/32\)。

## 证明余量审计结论

在 Lemma 2.6 的覆盖链中，平均值条件可写成严格不等式，但其严格差没有统一下界，不能把逐个 cube 的 `>` 变成与 \(f,t\) 无关的 \(\delta>0\)。同时，上一节的构型表明重叠估计 \(\sum1_Q\le32\) 本身在 Lebesgue 正测度区域严格取等；故无法利用“边界零测”把它替换为 \(31\)，或替换为任何 \(32-\delta\)。

这只否定了**从该 Besicovitch multiplicity 证明直接榨取 Lebesgue 严格余量**的路线；它不排除有另一种使用 Lebesgue 平移不变性、输入质量分配或多尺度结构的证明得到 \(C_5<32\)。要闭合严格改进，缺少的步骤正是一个对所有 \(f\in L^1(\mathbb R^5)\)、所有阈值和所有尺度的、绕开此 32 重正测度配置的 global covering/charging inequality。

## 实际运行摘要

运行：`python3 literature/verify_strict_upper_overlap_cube5.py`。

输出证书中记录：`number_of_cubes=32`、`minimum_pairwise_linf_center_distance=3/2`、`common_overlap_volume=1/32`。这是几何诊断/反例，不是新的 \(C_5\) 上界或下界证书。
