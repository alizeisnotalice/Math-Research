# 第六轮 E：两轴 parity 质量的独立严格接口

本件只处理六维 Lebesgue、中心轴平行、单一共同半径、一般 L¹、全部正尺度的下界接口。D 的本轮首批 12 个数值候选没有越过冻结完整 R6 上端，因此 E 未对这些候选运行完整并集体积认证，不登记新 C6 界，也不重跑旧 baseline 或 length-cap 扩展。

冻结文件为主 clone 的 `coordination/phase_20261025/frozen_thresholds.json`，实际读取 SHA-256 为 `e0528774b908a5135e24337dacfd6c15c5cc3e156da4578cd694fdece5eecd6d`；比较对象是 `lower.exact_interval[1]` 的完整有理数，不是简单显示下界。该分数完整记录在 `interface_selfcheck_results.json`。

## 1. 模型与周期反射

令所有 a_h>0，格点位置 A k=(a₁k₁,…,a₆k₆)，取不同的 i,j，

\[
w(k)=1+\delta(-1)^{k_i+k_j},\qquad -1\le\delta\le1.
\]

定义 p_i=p_j=2，其余 p_h=1。按物理向量 p_h a_h 的平移保持 w。一个粗周期含四个剩余类，质量为 4，因为四个 parity 权重的和为 4。每一坐标的反射 k_h↦−k_h 也保持 w。

记

\[
H=\prod_h[0,p_ha_h/2],\qquad
C=\prod_h[-p_ha_h/2,p_ha_h/2].
\]

H 的 64 个坐标反射像填满 C，内部互不相交。若 U⊂int H 的各点分别有一个共同半径 cube 包含质量至少 α(2r)⁶，则无限周期模型的形式比值是

\[
\alpha\frac{64|U|}{4}=16\alpha|U|.
\]

这里无需再除以格的协体积；a_h 已体现在物理体积 |U| 中。α=1 时就是 D 的 `16*union`。无限测度只用于说明几何；第 5 节给出一个实际有限 L¹ 输入。

## 2. 块质量与共同半径

给定整数起点 q 和正整数长度 ℓ，取指标块

\[
B(q,\ell)=\prod_h\{q_h,\ldots,q_h+\ell_h-1\}.
\]

令 S(q,L)=∑_{k=q}^{q+L−1}(−1)^k=(−1)^q 1_{L odd}。展开乘积可得准确质量

\[
m(q,\ell)=\prod_h\ell_h+
\delta S(q_i,\ell_i)S(q_j,\ell_j)
\prod_{h\ne i,j}\ell_h.
\]

该式包含起点 parity；只提供 ℓ 不能确定 m。零质量块直接丢弃。对非零块取 r=(m/α)^(1/6)/2，则块质量除以 cube 体积为 α；每一坐标都必须使用同一个 r。

选定块全部原子被 cube 包含的中心域为

\[
\prod_h[(q_h+\ell_h-1)a_h-r,\ q_ha_h+r].
\]

与 H 相交后是一个候选盒。允许 cube 包含其他原子，因为 w≥0；不需要证明所选块是 cube 内全部原子。有限盒族的并集只需作为最大函数超水平集的子集；本件不声称枚举了全部块或全部尺度。

## 3. Fraction 六次根、内外包与严格余量

`correlated_fraction_checker.py` 独立使用 Python 标准库整数与 Fraction，没有调用 D 的 volume 代码或任何作者模块。默认 spacing 直接按正运算包围

\[
u^6=2,\quad v^6=3,\quad
a=(1+u^4/2,(u^3+u^4v^3)/4,(1+u)/2,\ldots,(1+u)/2).
\]

对每个非负有理数 z，通过整数二分找到 n，逐项核验

\[
(n/D)^6\le z<((n+1)/D)^6.
\]

真实候选默认 D=10⁵⁰；未使用浮点根。下端 r₋ 给出准确有理共同半径，并准确验证 α(2r₋)⁶≤m。负的格点指标通过有符号区间乘法处理。

设 a_h∈[a_h⁻,a_h⁺]。候选中心区间的下端向上、上端向下包围，再在 G=10³⁵ 网格上各额外收缩一格。对每个保留盒，checker 重新核验两个中心端点与源指标块两个端点的所有差严格位于 (−r₋,r₋)，并核验盒严格位于 H 内部。坐标凸性保证中间的所有源格点也包含。记录全部盒的最小认证余量 μ>0。

同时用根上端 r₊、spacing 外包和向外网格取整构造外盒。由此得到所列理想盒族并集体积的有理下端与上端；即使上端低于门槛，也只排除这一个所列候选，不是 C6 上界或整个机制的排除。

整数并集采用端点切片：在首坐标全部整数端点间分条带，对贯穿条带的盒投影递归求并集。条带内部不交；末维按区间合并准确求长。先去重与有理包含只减少计算量。最终除以 G⁶ 即为准确有理体积。内外两次求并集共用 10⁶ 递归节点、20 分钟预算。若预算耗尽，输出 `budget_stop`，不声称通过。

## 4. 粗周期 collar 与准确边界损失

对所用块及其反射，令

\[
S_h=\max_{\rm used\ blocks}\max(|q_h|,|q_h+\ell_h-1|),\qquad
C_h=\lceil S_h/p_h\rceil.
\]

只接收粗周期指标 t∈[−K,K]⁶ 内的 64 个反射盒。其所需原子指标可写为 p_ht_h+s_h，|s_h|≤S_h。保留的有限源集合按完整粗周期取

\[
t'_h\in[-K-C_h,K+C_h],\quad b_h\in\{0,\ldots,p_h-1\},\quad
k_h=p_ht'_h+b_h.
\]

该集合在坐标 h 的指标范围是

\[
[-p_h(K+C_h),\ p_h(K+C_h+1)-1].
\]

因为 p_hC_h≥S_h，所有所需原子都在此范围。完整粗周期使有限源总质量准确为

\[
4\prod_h(2K+2C_h+1).
\]

接收盒总体积准确为 (2K+1)⁶64|U|。因此边界损失为

\[
L_K=\prod_h\frac{2K+1}{2K+2C_h+1}\longrightarrow1.
\]

这里不存在被忽略的 parity 边缘质量；它由完整周期的选择消掉。collar 是源指标的保留条件，不是对所有尺度的截断或 C6 上界。

## 5. 显式有限 L¹ 平滑与严格门槛

令 μ 为第 3 节余量，η=μ/2，

\[
b_\eta=(2\eta)^{-6}{\bf1}_{[-\eta,\eta]^6},\qquad
f=\sum_{k\in\text{finite sources}}w(k)b_\eta(\cdot-Ak).
\]

这是有限非负和，属于 Lebesgue L¹，其范数就是上节有限源质量。每个指定 cube 包含所需原子的整团平滑质量，因而所有接收点的 Mf≥α。

若有理内缩比值 R₋=16α|U| 严格大于冻结完整上端 T，则取

\[
\beta=\alpha\frac{R_-+T}{2R_-}\in(0,\alpha).
\]

接收点均有 Mf>β。checker 从 K=1 开始倍增，以 Fraction 验证

\[
\frac{\beta}{\alpha}R_-L_K>T.
\]

左端正是该显式有限 L¹ 输入在水平 β 的弱型比值下界。故 C6≥左端>T。有限源无需实际枚举；有限集合、权重、K、η、β 均有明确公式。完整全部正尺度最大算子包含这些有限个共同半径，传递到全尺度下界没有额外假设。

## 6. 实际运行与复现

本轮仅运行 `--selfcheck`：6 个有理根例；120 个块质量公式与直接逐原子求和比较；10 个 periodmass=4 例；K=0,C=1 的 2916 原子有限粗周期质量；准确已知两六维盒并集 127；单位 spacing、δ=1/3 的两个单原子盒测试整个内外包与 finite-L¹ 流程。最后一个仅是接口例，目标 1/10，并非本研究的新 C6 候选；实际输出 K=1、源质量 8100、collar=(1,1,0,0,0,0)、损失 9/25、正 μ、正 finite margin。

```bash
PYTHONDONTWRITEBYTECODE=1 python3 cube-dim-6-bounds/audit/lower006/correlated_fraction_checker.py \
  --selfcheck \
  --frozen /path/to/cube-dim-6-bounds/coordination/phase_20261025/frozen_thresholds.json
```

未来只对 D 晋级的最佳一个候选运行：

```bash
PYTHONDONTWRITEBYTECODE=1 python3 cube-dim-6-bounds/audit/lower006/correlated_fraction_checker.py \
  --candidate /path/to/candidate.json \
  --frozen /path/to/frozen_thresholds.json
```

候选 JSON 必需项：`axes` 为两个 0-based 整数；`delta`、`level` 为精确有理数字符串（默认 level=1）；`spacing` 为 `frozen_R6_003`；`blocks` 中每项含六个整数 `start` 与六个 1..4 整数 `lengths`。所有输出限定在本件目录。真实新界还需主代理整合模型、准确比较和独立复核。
