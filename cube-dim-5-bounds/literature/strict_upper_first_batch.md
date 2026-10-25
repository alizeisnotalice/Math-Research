# Lebesgue 严格上界：首批核验（登记）

## 问题与验收标准

固定 \\(d=5\\)，令 \\(M_5\\) 为 Lebesgue 测度上、以每个点为中心、坐标轴平行且**所有坐标共用同一半径**的立方体极大算子。首批只接受能推出

\\[
  \\lambda^5\\{M_5f>\\alpha\\}
  \\le U_5\\alpha^{-1}\\lVert f\\rVert_1,
  \\qquad U_5<32,
\\]

且明示适用于一般 \\(L^1(\\mathbb R^5)\\) 输入、全部正半径的公开原始文献结论。隐含常数、渐近 \\(O(d\\log d)\\) 而无可追踪数值系数、欧氏球、非中心或一般度量空间结论均不计作本题严格上界。

## 已执行查询（2026-10-25，8 条上限内）

1. `Stein Strömberg Behavior of maximal functions in R^n for large n cubes weak type explicit constant pdf`
2. `centered Hardy Littlewood maximal cubes Lebesgue weak type 2^n strict upper bound Stein Stromberg`
3. `Melas centered Hardy Littlewood maximal operator cubes higher dimension weak type constant`
4. `Aldaz centered maximal operator cube weak (1,1) upper bound 2^n Lebesgue`
5. `Naor Tao random martingales localization maximal inequalities cube maximal weak type n log n pdf`
6. `Besicovitch covering centered cubes overlap constant 2^n Lebesgue maximal weak type`
7. `"centered" "cube" maximal "weak type" "n log n"`
8. `"centered maximal operator" cubes "5" weak type`

检索来源：公开网页/出版商索引、arXiv 及作者公开 PDF；检索返回只用于定位候选，定理适用性须以下列原文阅读为准。

## 首批原件候选（最多 6 件；待逐件阅读）

1. E. M. Stein and J.-O. Strömberg, *Behavior of maximal functions in \\(\\mathbb R^n\\) for large \\(n\\)*, Arkiv för Matematik 21 (1983), 259--269, DOI 10.1007/BF02384321.
2. A. Naor and T. Tao, *Random martingales and localization of maximal inequalities*, JFA 259 (2010), 731--779, DOI 10.1016/j.jfa.2009.12.009; author PDF / arXiv:0912.1140.
3. A. D. Melas, *The best constant for centered Hardy--Littlewood maximal inequality*, Ann. Math. 157 (2003), 647--688, DOI 10.4007/annals.2003.157.647.
4. A. S. Besicovitch, *A general form of the covering principle and relative differentiation of additive functions, I*, Proc. Cambridge Philos. Soc. 41 (1945), 103--110, DOI 10.1017/S0305004100022453.
5. J. M. Aldaz, *The weak type (1,1) bounds for the maximal function associated to cubes grow to infinity with the dimension*, Ann. Math. 173 (2011), 1013--1023, DOI 10.4007/annals.2011.173.2.10 (此前已读；此轮不重复作为 ``2^d`` 来源)。

## 停止条件

达到 6 件原始文献或以上 8 条查询的全部候选已处置即停止首批；若没有满足验收标准的显式 \\(U_5<32\\)，只报告“本批未找到”，绝不推出不存在。

## 实际原文/题录核验结果

### 1. Naor--Tao（全文，相关但不产生数值门槛）

作者公开 PDF（41 页，SHA-256 `0d67149610ae88ceccc02822be7e1dc921a692b334eb09e4c5e86ca55857fe8c`）的导言 p.2 和 Corollary 1.2（印刷 p.5--6）确实把 \\(n\\) 维赋范空间的 centered ball maximal operator 的弱 \\((1,1)\\) 上界描述为 \\(O(n\\log n)\\)。对强 \\(n\\)-microdoubling 空间，原文的式 (15) 是

\\[
\\lVert M\\rVert_{L^1\\to L^{1,\\infty}}\\lesssim K n\\log n.
\\]

它包含 \\(\\ell_\\infty^5\\) 的几何背景：半径放大 \\(1+1/5\\) 时 Lebesgue 球体积比为实际计算得到的

\\[
 (1+1/5)^5=2.48832.
\\]

但是符号 \\(\\lesssim\\) 的绝对常数在该表述中未显式给出；所以这不是可审计的 \\(U_5\\)，更不能推出 \\(U_5<32\\)。不能以 \\(5\\log 5\\approx8.04719\\) 代替被隐藏的乘子。

### 2. Stein--Strömberg（题录已核验；本批未获可读原文）

Crossref 的 Arkiv 1983 卷目给出正确 DOI 为 [10.1007/bf02384314](https://doi.org/10.1007/bf02384314)（此前容易混淆成错误的末位 `321`）。题名在其元数据中有印刷拼写 `functios`，卷 21(1--2), 259--269。Naor--Tao 的原文导言将该工作归为一般 \\(n\\) 维赋范空间 \\(O(n\\log n)\\) 结果；该二手定位不替代对 Stein--Strömberg 定理常数的原文认证。本批出版商 PDF 端点没有返回文章 PDF，故**没有**从此条导出数值 \\(U_5\\)。

### 3. Melas（全文，明确排除）

Annals 原文 PDF（42 页，SHA-256 `f82530d32bab47bfb8a16e767484d5768c307af6cf1987c2a1c450ff6a938879`）摘要和 p.647--648 都明确研究 \\(\\mathbb R\\) 中一维 centered interval operator，并给出其精确常数为 \\(12C^2-22C+5=0\\) 的大根 \\(C\\approx1.5675208\\)。它不是五维 cube 的弱 \\((1,1)\\) 上界。逐坐标迭代一维算子只产生复合弱型/Orlicz 控制，不能在这里合法地推出同为 \\(L^1\\to L^{1,\\infty}\\) 的五维常数，故排除。

### 4. Besicovitch 覆盖（题录已核验；不作为 Lebesgue 改进）

1945 原文的 Crossref 题录已核验（哈希见下），但本批未获全文。它是覆盖方法的历史来源，不是本题所需的、带显式 \\(U_5<32\\) 的 Lebesgue cube 定理。既有台账中已经单列了 Aldaz 的公开一般测度 \\(\\ell_\\infty^d\\) sharp \\(2^d\\) 结果；该“对所有测度”结论不能改写成 Lebesgue 下的等号或严格不等式，本轮亦未重复把它登记为新来源。

## 首批结论

在登记的 8 条检索和上述 4 项原始候选中，**未找到**可由原文精确代入、同时覆盖一般 \\(L^1(\\mathbb R^5)\\) 与全部尺度的严格公开界 \\(U_5<32\\)。这只是受限首批的否定检索结果，不是此类结果不存在的证明，也不改动现有 \\(C_5\\le32\\) 门槛。

## 题录与可复核源哈希

| 项 | 核验入口/版本 | 实际阅读 | SHA-256 |
|---|---|---|---|
| Naor--Tao | [作者 PDF](https://web.math.princeton.edu/~naor/homepage%20files/HLM.pdf); JFA 259 (2010), 731--779, DOI [10.1016/j.jfa.2009.12.009](https://doi.org/10.1016/j.jfa.2009.12.009) | 全文（定义、导言、(15)、参考文献） | `0d67149610ae88ceccc02822be7e1dc921a692b334eb09e4c5e86ca55857fe8c` |
| Melas | [Annals PDF](https://annals.math.princeton.edu/wp-content/uploads/annals-v157-n2-p08.pdf); Ann. Math. 157 (2003), 647--688, DOI [10.4007/annals.2003.157.647](https://doi.org/10.4007/annals.2003.157.647) | 全文（摘要、p.647--648、参考文献） | `f82530d32bab47bfb8a16e767484d5768c307af6cf1987c2a1c450ff6a938879` |
| Stein--Strömberg | Arkiv för Matematik 21 (1983), 259--269, DOI [10.1007/bf02384314](https://doi.org/10.1007/bf02384314) | 题录；原文 PDF 端点本批未成功 | Crossref 卷目 `1191e2d947739a1050d4dfde5d754969a3bf1c16cdb6b45ba1af19fdcfa2eacb` |
| Besicovitch | Proc. Cambridge Philos. Soc. 41 (1945), 103--110, DOI [10.1017/S0305004100022453](https://doi.org/10.1017/S0305004100022453) | 题录；未以其导出数值 | Crossref 记录 `b812f72806db4c613e8a5f999cad2e2d977342575038d2f0bbc67a9878332244` |

Python `Decimal` 实际计算（80 位精度）输出：`(1+1/5)^5 = 2.48832`、`5*ln(5) = 8.047189562170501873003796666130938...`；后者只说明量级表达式的裸因子，绝非上界证书。
