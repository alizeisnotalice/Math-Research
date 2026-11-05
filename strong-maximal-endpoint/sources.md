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
