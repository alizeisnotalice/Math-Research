# 六维公开文献账本：首轮定位

检索日期：2026-10-25（Asia/Shanghai）。对象固定为 Lebesgue 测度上的中心、轴平行、共同半径方体极大算子；所有尺度及一般 \(L^1(\mathbb R^6)\) 输入。

## 适用的公开数值门槛

* **下界。** Aldaz, *A remark on the centered \(n\)-dimensional Hardy--Littlewood maximal function*, Czechoslovak Math. J. 50 (2000), 103--112, Proposition 1.4，给出每个 \(n\ge2\) 的构造下界
  \[
  A_n=\frac32a^{n-1}-\frac12\left[(3a-2(3/2)^{1/n})^{n-1}-((3-2^{1/n})/2)^{n-1}\right],\quad
  a=(1+2^{1/n})/2.
  \]
  在 \(n=6\)，`check_aldaz_n6.py` 用纯整数/有理区间证明
  \[
  C_6\ge A_6,\qquad
  A_6\in(1.763857472539509489679047135,\;1.763857472539509489679047136).
  \]
  因而可安全写成严格有理门槛 \(C_6>1.763857472539509489679047135\)，但不得将命题本身误写成 \(C_6>A_6\)。
  原论文证明先给原子测度见证，再以常规平滑传递到 Lebesgue 的一般 \(L^1\) 弱型常数；因此范围相符。
* **上界。** Aldaz, *Kissing numbers and the centered maximal operator*, arXiv:1811.10478v2（2020）。Theorem 1.2 把全体 \(\tau\)-additive 局部有限 Borel 测度上的统一中心弱型常数识别为 Besicovitch 常数；Theorem 3.2 识别该常数与严格 Hadwiger 数；**Corollary 3.5** 专门给出 \(L(\mathbb R^d,\|\cdot\|_\infty)=2^d\) 且该统一界达到。Lebesgue 测度适用，因此 \(C_6\le64\)。这是公开可核实的有效起点，但本轮没有证明它是六维最强公开上界。

在同一 Aldaz 公式中，严格有理区间算术也证明 \(A_5<1.749<1.763<A_6\)。所以公开的 Aldaz 五维公式经维数单调性传递，不能改进其直接六维公式。

## 不可当作六维准确门槛的高维结果

* Aubrun (2009) 给出随维数渐近的 \((\log n)^{1-o(1)}\) 型下界；Iakovlev--Strömberg, *Lower bounds for the weak type (1,1) estimate for the maximal function associated to cubes in high dimensions*, Math. Res. Lett. 20 (2013), 907--918, DOI 10.4310/MRL.2013.v20.n5.a7，给出 \(c n^{1/4}\) 型下界。已核的摘要/综述没有给出可在 \(n=6\) 代入的明示绝对常数，因此不能将它们转写为六维精确下界，也不能据此与 Aldaz 的 \(A_6\) 比较。
* 球、独立矩形尺度、强型 \(L^p\) 或其他背景测度的结果均不进入本账本。
* Lin 的 `Two partial balayage principles for weak-type estimates` 本轮实际读取目录、总表及 §4.3。其方体条目仅是二维轴平行正方形；一般维数的公式均属于欧氏球、Poisson 或热核，未出现一般 \(n\) 方体，更没有 \(n=6\) 方体上界。因此该网页提供方法背景，不进入 \(C_6\) 数值上界账本。

## 本批实际来源与复核

* 本轮直接下载并实际阅读原始 PDF，而非将资料包操作记录当作新证据；字节 hash、URL、版本、页码和阅读结论见 `primary_source_receipt.json`。Aldaz 2000 PDF 的 Proposition 1.4 在印刷页 110（PDF 第9页），其紧前 Lemma 1.1 在印刷页104--105（PDF 第3--4页）给离散见证到 \(L^1\) 的传递。Kissing v2 的 Theorem 1.2 在 PDF 第2页，Theorem 3.2、3.3 与 Corollary 3.5 在 PDF 第10--11页。
* 本轮用出版方/档案页与 arXiv 检索复核 Aldaz 2000、Aubrun 2009、Aldaz 2011、Iakovlev--Strömberg 2013 的题录和对象范围。该检索是定向批次，不能证明不存在其他更强的固定六维结果。
* 第二批恰用六条定向查询：Lin 一般方体范围、`C_6`、六维固定维数、Aldaz \(A_5\)、五维固定维数、六维固定维数。前四条没有定位可用于本算子的额外固定六维常数；后两条仅回到 Lin 的二维方体与无关命中。该负结果只限定本批检索范围。
* 运行：`/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python literature/check_aldaz_n6.py`。停止条件：得出 Aldaz 六维严格有理包围；该条件已满足。该脚本不建立 Iakovlev--Strömberg 的未显式常数。

## 六维研究的比较口径

公开冻结起点暂记为 \(A_6\) 和 \(64\)，但在将其作为“当前最好公开界”前，需要由主账本将本轮的定向文献结果与更广的最新原始来源继续交叉核验。内部 R4 证书经维数单调性可以另行形成更强的**内部**下界种子；它不是公开发表结果，也不应混入本表。
