# 强极大算子端点常数：原始文献核查（2026-11-06）

## 核验范围和结论

本记录只核验三类公开、可追溯的结果：一维及强极大算子的锐 \(L^p\) 常数，固定维数的\(L\log L\)端点弱型/小阈值结果，以及 Solyanik（Tauberian）常数。它们都不能单独证明本项目关于 \(n\to\infty\) 的固定阈值 Orlicz 常数 \(C_n\) 的次指数渐近。

| ID | 状态 | 与本项目可安全使用的关系 |
|---|---|---|
| L01 | 已核验全文 | \(\|M_S\|_{p\to p}=c_p^n\)，其中 \(c_p\) 是显式方程的唯一正根；这是一般输入上界路线的锐 \(L^p\) 输入。 |
| L02 | 已核验版本记录、预印本全文 | 固定 \(n\) 的 \(L_{\widetilde\Phi_n}\to L_{\widetilde\Phi_n}^{1,\infty}\) 弱型和小阈值极限；模数和弱范数均不同于本项目。 |
| L03 | 已核验预印本全文 | 示性函数的 \(\alpha\uparrow1\) Tauberian 常数；可作近阈值测试族的背景，不能直接代替一般函数的 Orlicz 常数。 |

## L01：Grafakos--Montgomery-Smith 的锐强型结果

原稿的首页定理给出：对 \(1<p<\infty\)，一维**非中心** Hardy--Littlewood 最大算子 \(M_1\) 的 \(L^p\) 算子范数 \(c_p\) 是
\[
(p-1)c_p^p-pc_p^{p-1}-1=0
\]
的唯一正根（论文第 1 页，方程 (1)、主定理）。论文随后在第 3 页定义轴平行、各边独立的非中心强极大算子 \(N_n\)，证明 \(N_n\le M_1^{(1)}\circ\cdots\circ M_1^{(n)}\)，并以张量化近极值函数证明反向不等式。因此推论为
\[
\|N_n\|_{L^p(\mathbb R^n)\to L^p(\mathbb R^n)}=c_p^n.
\]

**证明位置。** [作者公开 PDF](https://stephenmontgomerysmith.github.io/preprints/uncentered.pdf)，第 1 页主定理和第 3 页 "Corollary"；在线文本的方程 (1)、主定理与第 3 页推论分别可直接定位于该 PDF 的第 1、3 页。论文的张量下界说明这里确为强矩形算子，不只是坐标复合的上界。

**本项目的合法用途。** 令 \(b_p=p\log c_p\) 后，题设中的方程
\[
p(1-e^{-b_p/p})=1+e^{-b_p}
\]
只是上述根方程的代换，故可作为待做渐近展开的正确起点。这个定理本身是强型 \(L^p\) 定理；它不表述固定阈值 \(1\) 的 \(\int \Phi_n(f)\) 最优常数，更不表述其 \(n\)-渐近的第二项。

## L02：Qin--Wu--Xue 的端点量与本项目不同

Qin--Wu--Xue 在定义 1.1 使用
\[
\widetilde\Phi_n(t)=t\,[1+(\log^+t)^{n-1}],
\]
并定义新的弱空间 \(L^{1,\infty}_{\widetilde\Phi_n}\)。其定理 1.1 研究的是该弱空间的限制弱型行为；第 2 节的 (2.5)、(2.9) 给出小阈值 \(\lambda\downarrow0\) 的上下极限，主系数为 \(2^n/(n-1)!\)。期刊版本信息为 *Banach Journal of Mathematical Analysis* 15 (2021), Article 69, DOI 10.1007/s43037-021-00154-6。

这与本项目必须分开：

1. 本项目使用 \(\Phi_n(t)=t(1+\log^+t)^{n-1}\)，而该文使用 \(\widetilde\Phi_n(t)=t[1+(\log^+t)^{n-1}]\)。当 \(n>2\) 且 \(t>1\) 时两者不相同。
2. 该文的常数属于 \(L_{\widetilde\Phi_n}\to L^{1,\infty}_{\widetilde\Phi_n}\) 的弱准范数，且首先取固定 \(n\)、阈值趋于零；本项目的 \(C_n\) 是固定阈值 1 后再考察 \(n\to\infty\) 的单一分布比值。
3. 因而 \(2^n/(n-1)!\) 可作为小阈值/示性函数对照，不是本项目 \(C_n\) 的已证等式、也不应与 \(B^n\) 的指数率比较成同一个端点常数。

**证明位置。** [arXiv 全文](https://arxiv.org/abs/2104.03497)：定义 1.1、定理 1.1、Corollary 1.1，以及第 2 节 (2.5)、(2.9)。[期刊版本记录](https://link.springer.com/article/10.1007/s43037-021-00154-6) 核对了作者、卷、文章号、日期与 DOI。

## L03：Solyanik/Tauberian 量也须单列

Hagelstein--Parissis 定义的锐 Tauberian 常数是
\[
C_{\mathcal B}(\alpha)=\sup_{0<|E|<\infty}
\frac{|\{M_{\mathcal B}\mathbf1_E>\alpha\}|}{|E|},\qquad0<\alpha<1.
\]
文章明确以 \(C_S\) 指代轴平行矩形（强极大算子）的该常数，并在其 Theorem 1 回顾 Solyanik 型结果，在 Theorem 2 处理坐标一维最大算子的复合。此量只允许输入为示性函数，研究 \(\alpha\uparrow1\)；它不能给出一般 \(f\) 的固定阈值 Orlicz 常数，也不能允许把维数相关隐常数忽略后转用到 \(n\to\infty\)。

**证明位置。** [arXiv 全文](https://arxiv.org/abs/1310.3771)，第 1 节方程 (1) 与第 2 节 Theorem 2；文中第 1 节还明确定义 \(C_S\) 为轴平行矩形对应的 Tauberian 常数。

## 可直接导入的 BibTeX

```bibtex
@article{GrafakosMontgomerySmith1997,
  author  = {Grafakos, Loukas and Montgomery-Smith, Stephen},
  title   = {Best Constants for Uncentred Maximal Functions},
  journal = {Bulletin of the London Mathematical Society},
  volume  = {29},
  number  = {1},
  pages   = {60--64},
  year    = {1997},
  doi     = {10.1112/S0024609396002081}
}

@article{QinWuXue2021,
  author  = {Qin, Moyan and Wu, Huoxiong and Xue, Qingying},
  title   = {The Limiting Weak-Type Behaviors of the Strong Maximal Operators},
  journal = {Banach Journal of Mathematical Analysis},
  volume  = {15},
  pages   = {69},
  year    = {2021},
  doi     = {10.1007/s43037-021-00154-6}
}

@misc{HagelsteinParissis2013,
  author       = {Hagelstein, Paul A. and Parissis, Ioannis},
  title        = {Solyanik Estimates in Harmonic Analysis},
  year         = {2013},
  eprint       = {1310.3771},
  archivePrefix= {arXiv},
  primaryClass = {math.CA},
  url          = {https://arxiv.org/abs/1310.3771}
}
```

## 检索日志与限制

- 2026-11-06（Asia/Shanghai）：直接打开 Grafakos--Montgomery-Smith 作者公开 PDF，读取全文的主定理、证明及强算子推论。
- 直接打开 Qin--Wu--Xue arXiv 全文与 Springer version-of-record 页面；前者用于定义、定理和公式定位，后者用于题录字段。
- 直接打开 Hagelstein--Parissis arXiv 全文，读取定义和 Theorems 1--2 的范围。
- 本记录没有对公开文献作穷尽性优先权检索，也没有对项目内部的 \(B\) 或三项公式作网页反向检索；因此不作 "文献首创" 的结论。


## 第5轮：独立非同分布的 Berry–Esseen 绝对误差

直接核验 [Vershynin, A friendly proof of the Berry-Esseen theorem, arXiv:2602.06234v1](https://arxiv.org/html/2602.06234v1)，Theorem 0.1 及第5节证明。
该版本日期为2026-02-05，作者 Roman Vershynin；所用版本固定v1，未混用后续版本的标题。
定理允许独立、非同分布、均值零、总方差归一为1的变量；误差由三阶绝对矩之和乘绝对常数控制。
因此对本项目总方差至少n、逐坐标三阶中心矩一致有界的三角阵，给出一致 O(n^-1/2) 的绝对分布误差。
第5轮只在高斯尾严格大于该绝对误差的窗口使用，不据此宣称临界尾的相对误差。
阅读范围：直接读取作者论文HTML的定理及证明；查阅日2026-11-06。

历史来源 Berry 1941 的 DOI 10.1090/S0002-9947-1941-0003498-3 可追溯；AMS原PDF读取受阻，未声称读过其原证明。
本轮可直接核验的证明来源使用上述作者预印本，而非受阻原文的二手摘录。
实际检索：非同分布 Berry–Esseen 三阶绝对矩；定向查找 Vershynin 题名与 arXiv v1。


## 第8轮：Hagelstein--Parissis Theorem 2 的显式维数依赖（2026-11-06）

这里仅考察该文的坐标复合算子 \(M_1\cdots M_n\)，并用 \(M_S\leq M_1\cdots M_n\) 把所得示性函数上界转给强矩形算子。原文 Lemma 1 的一维常数是完全显式的：对 \(0\leq\gamma<\alpha<1\)，
\[
 |\{M_{\rm HL}(\mathbf1_E+\gamma\mathbf1_{E^c})>\alpha\}|
 \leq \left(1+4\frac{1-\alpha}{\alpha-\gamma}\right)|E|.
\]
令 \(\alpha_j=1-(1-\alpha_1)^j\)。其 Lemma 2 逐步迭代而不丢弃常数，得到
\[
 |\{M_1\cdots M_n\mathbf1_E>\alpha_n\}|
 \leq\left(1+4\frac{1-\alpha_1}{\alpha_1}\right)^n|E|.
\]
故 Theorem 2 对**每个** \(0<\alpha<1\)（没有 \(n\)-一致的小邻域阈值假设）给出
\[
 C_{1\cdots n}(\alpha)\le U_n(\alpha)^n,\qquad
 U_n(\alpha)=1+4\frac{(1-\alpha)^{1/n}}{1-(1-\alpha)^{1/n}}.
\]
同样有 \(C_S(\alpha)\le C_{1\cdots n}(\alpha)\)。原文在这之后只写
\(C_{1\cdots n}(\alpha)-1\sim_n(\alpha^{-1}-1)^{1/n}\)；符号 \(\sim_n\) 明言隐常数可以依赖于 \(n\)，原文未给出可跨维使用的阈值或常数。

将 bounded-core reduction 的阈值写为
\[
q_n=\alpha_n^{-1}-1=
\left(\frac{e^{A\sqrt{\log n/n}}}{1-2/\sqrt n}-1\right)^{1/2},
\qquad \delta_n=1-\alpha_n=\frac{q_n}{1+q_n}.
\]
则显式定理实际给出的 \(n\) 次方根是
\[
U_n(\alpha_n)=1+\frac{4}{e^{L_n/n}-1},\qquad L_n=-\log\delta_n.
\]
对固定 \(A>0\)，
\[
q_n\sim\sqrt A\,(\log n/n)^{1/4},\quad
L_n=\tfrac14\log n-\tfrac14\log\log n-\tfrac12\log A+o(1),
\]
进而
\[
U_n(\alpha_n)=1+\frac{16n}{\log n}(1+o(1)).
\]
这个根趋于无穷，远大于 \(B=3.591121476\ldots\)，所以 Theorem 2 的显式上界在该阈值完全不能提供 \(B^n\) 级控制。

更直接地，若只要求该显式根小于 \(B\)，必要且充分的代数条件是
\[
(1-\alpha)^{1/n}<\frac{B-1}{B+3}=0.3931230043\ldots,
\quad\text{即}\quad
1-\alpha<\left(\frac{B-1}{B+3}\right)^n.
\]
这要求阈值余量指数小于 \(e^{-0.9336327280\ldots n}\)；当前 \(\delta_n\asymp(\log n/n)^{1/4}\) 相差指数级。故不能以文中固定 \(n\) 的 \(\sim_n\) 断言替代这一显式计算。

**本机数值核对。** MATLAB R2026a 以 \(A=3\) 直接代入上述精确公式：

| \(n\) | \(\alpha_n\) | \(U_n(\alpha_n)\)（即该界的 \(n\) 次方根） |
|---:|---:|---:|
| 16 | 0.290370384712 | 185.589477818 |
| 64 | 0.422742294889 | 464.909351057 |
| 256 | 0.531459562798 | 1349.68770965 |
| 1024 | 0.623291718450 | 4194.50003364 |
| 4096 | 0.701143511562 | 13564.2520734 |

这只是对原文显式上界的数值代入，不是对真实 \(C_S\) 或 bounded-core 问题的下界/反例。

**原始位置。** [Hagelstein--Parissis arXiv:1310.3771](https://arxiv.org/abs/1310.3771)，第 2 节 Lemma 1、Lemma 2、Theorem 2 及其紧随的证明。为避免网页数学渲染遗漏公式，本轮直接下载并读取作者上传的 arXiv 源码 `HagelsteinParissis2013Revised.tex.tex`；上述三个公式分别是原文的 (2)、Lemma 2 和 Theorem 2。


## 第13轮：条件Bernstein鞅的相对中偏差

直接核验Fan, Grama, Liu, *Cramér large deviation expansions for martingales under Bernstein’s condition*，arXiv:1210.2198v3（2014-09-14，26页），期刊SPA123(2013),3919–3942，DOI10.1016/j.spa.2013.06.010。
[原文第2节](https://arxiv.org/html/1210.2198v3) Theorems2.1—2.2、式(14)及第5—6节证明(54)—(56)已读取。
A1要求每个k≥3的条件有符号幂矩绝对值由(k!/2)epsilon^(k−2)乘条件方差控制；条件绝对矩的同界是更强的充分条件。A2要求单位归一化鞅的可预测二次变差偏差≤delta²几乎处处，而非仅无条件方差正确。
定理给正态尾比率的双边控制，对数误差至多常数乘x³epsilon+x²delta²+(1+x)(epsilon|log epsilon|+delta)。固定绝对alpha0的范围x≤alpha0/epsilon、delta≤alpha0内合法。
第13轮严格建立epsilon=32/(sigma sqrt n)、delta²=O(n^-1/4(log n)^1/4)，故x≤2sqrt(log n)内相对误差O(n^-1/8(log n)^5/8)。须把delta取为A2右端的平方根。原文未数值化绝对常数，不能生成有限维进入渐近区间的证书。
版本题录另见[arXiv版本记录](https://arxiv.org/abs/1210.2198)。查阅日期2026-11-06；不以该概率论来源主张本项目端点渐近的公开优先权。


## 第27轮：对数凹的边缘积分与卷积封闭性

András Prékopa, *On logarithmic concave measures and functions*, Acta Scientiarum Mathematicarum 34 (1973), 335–343。已直接读取[作者机构保存的原文](https://rutcor.rutgers.edu/Prekopa/pdf/SCIENT2.pdf)第4节Theorems6–7（印刷页342–343）：非负对数凹函数的凸域边缘积分仍对数凹，两个对数凹函数的卷积仍对数凹。第27轮仅使用这两个一般封闭性；特定Gamma–Bernoulli核的四阶起点由本轮有理多项式证书另证。查阅日期2026-11-06，不据此主张特定密度结论的公开优先权。


## 第54轮计算与内部来源

有限输入的几何转移接续r002/r009/r028；完整KKT、Ward与响应接续r050/r052/r053。本轮证明和证书见rounds/r054/，不新增公开优先权归属。定向算术使用 [Python Decimal官方文档](https://docs.python.org/3/library/decimal.html) 的exp/ln正确舍入和next_minus/next_plus契约，整数幂采用重复区间乘法。Wolfram 15的120/240位重放仅作数值独立交叉验证。
