# n=8 中心立方体弱型常数：检索记录

日期：2026-10-27（Asia/Shanghai）。实际检索工具：Codex `web.run`；未调用未连接的学术数据库。目标主题限定为 Lebesgue 测度上、轴平行中心立方体、弱型 (1,1)、固定维 n=8。检索范围有限，不声称穷尽。

## 实际检索与阅读

1. 查询 `Aldaz 2000 proposition 1.4 centered Hardy-Littlewood maximal cube A_n weak type constant pdf`；打开 Czech Digital Mathematics Library 的原始 PDF：<https://dml.cz/bitstream/handle/10338.dmlcz/127554/CzechMathJ_50-2000-1_14.pdf>。核对题录首页、定义段、Prop. 1.4 与证明（PDF 纸本页 110–111 / PDF 页 8–9）。结论是每 n≥2 的显式下界公式。用 Python `Fraction` 的区间算术，将 2^(1/8)、(3/2)^(1/8) 各自夹在 10^-50 网格相邻整数根之间，代入原文公式并向外传播，得到可复核粗区间 1.78390524162816 < A8 < 1.78390524162817。根界来自整数八次幂比较，不由十进制近似回造。
2. 查询并打开 arXiv v2 版本页 <https://arxiv.org/abs/1811.10478v2>（确认作者、题名、v2 日期 2020-03-10），再打开同版 HTML <https://arxiv.org/html/1811.10478v2>，阅读定义和 Corollary 3.5（HTML 行 296–301）。其一致弱型常数对中心 ℓ∞ balls 正好为 2^d，故对 d=8 给出上界 256，适用于 Lebesgue；该锐性陈述是对所有背景测度统一常数的锐性，不是 Lebesgue 最优常数的刻画。
3. 查询 `Aldaz 2011 Annals centered cubes weak type dimension monotonicity`；打开 Annals 论文 PDF/页面 <https://annals.math.princeton.edu/wp-content/uploads/annals-v173-n2-p10-p.pdf> 与 <https://annals.math.princeton.edu/2011/173-2/p10>。检索结果摘要明确原文定义为全部 r>0 的中心 ℓ∞ balls，且记载 C_(d+1)≥C_d；用户提供的资料包也记录原文第 1013 页引用 AV07 Th. 2。此轮另外给出独立的乘积截断证明，见 `cube-dim-8-bounds/baseline/round001/inheritance.md`。
4. 检索 Aubrun 2009 的原始期刊页/作者 PDF：<https://cml.centre-mersenne.org/articles/10.1142/S1793744209000067/>，<https://math.univ-lyon1.fr/~aubrun/recherche/cubicmaximalfunction/maximalfunction.pdf>。阅读期刊元数据与摘要：最优常数随 n 的下界增长快于 (log n)^(1-o(1))。这是渐近结论；此次没有提取可用于 n=8 的显式门槛/数值。
5. 查询 Iakovlev–Strömberg 2013；检索到出版方 PDF <https://intlpress.com/api/bgcloud-front/resource/pdf/volume/1806603729623216130-1806603729623216130-88dbe922fef07a892fbe2fa55c340767.pdf>。搜索摘要给出论文报告 C_n≥c n^(1/4)，其中 c 与维数无关；web.open 该 PDF 返回 HTTP 403。另查到 Deleaval–Guédon–Maurey 后续论文 PDF <https://www.numdam.org/item/10.5802/afst.1567.pdf> 对该路线的复述及其自身 n>256 条件；这不足以直接给 n=8 可用数值。本记录因此写“未提取”，不写“无固定维结果”。
6. 查询 Tomar DOI `10.5281/zenodo.23554628`；web.open <https://zenodo.org/records/23554628> 报 URL 不可访问，搜索引擎也没有可靠命中该 DOI。本地用户资料包中 `cube-finite-dim-bounds/literature/tomar-v4-source.json` 记录官方 API 在 2026-10-25 取回的版本 4 元数据及 PDF SHA-256 `aa32aab64b2d973d3bbcbdd12e309125c28247284c739758c24abfc2bb9e2ae6`，且记录仅阅读二维 squares 论文 p.1、p.2 定理 1.3、p.10 Lemma 6.1–6.2、pp.10–11 证明。本轮把它作为继承资料包的可追溯二维候选，不声称重新访问或重新阅读全文。

## 纳入/排除与局限

纳入：Aldaz 2000 固定维下界，Aldaz 2011 单调性，Aldaz 2020 一致上界，Aubrun / Iakovlev–Strömberg 渐近定量来源。Tomar DOI 项纳入候选账目但限定为二维背景。

排除：欧氏球新文、Poisson/heat kernel、强型或其他背景测度的结果，不直接视为中心立方体 Lebesgue 弱型常数的新 n=8 界。Aldaz Corollary 3.5 虽给最优的一致背景常数 256，其锐性不能转化为 Lebesgue 上界的最优性。当前固定维下界主要来自内部 n=7 证书沿维度单调性传递，并非上述渐近文献的 n=8 取值。
