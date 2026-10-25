# CI4：四维正 face octant 的连续积分诊断（未闭合）

2026-10-25。目标仍未证明：

\[
 I_4=\int_{[0,1]^3}G(1,y,z,w)^{4/3}\,dy\,dz\,dw\le19/100000.
\]

固定八条未单位化 body 方向、Lebesgue 背景、共同半径、一般 L1/全尺度；没有 tau、g 或成本搜索。本文不更新 U4=16。数值候选 .000183520694694 仅来自主输入；本工作没有将有限点、Gauss 阶差或浮点精度当作连续认证。

登记 `ci4_registration.md` SHA256 `db6154cf407f97687f928d253a5ab6d0cc4c41f605297abcabdcaf417cbc6bdb`。本批只做一个初始有限几何诊断及一次有明确失败原因的相关性修订，均为同一组 351 区域，没有扩网格。完成了本机 Python Fraction / python-flint Arb 的实际运行；未占用共享 MATLAB 会话。

## 结果与未支付接口

- 精确几何：351 个凸多面体，最多 8 顶点，全部谱节点符号在各区域内固定。切分、交点、facet hull、体积均使用有理数；没有以浮点凸包作认证。
- 边角遗漏：epsilon=10^-10 时，正 face octant 的真实方向线帽预算小于 `2.673612529406e-12`。这是真实遗漏积分的**上界**，其 Arb 包围的下端不是遗漏积分下界。
- 唯一相关性修订：33,840 对边际区间看似为零的 cross denominator 中，24,952 对由联合 affine 顶点 hull 解除了假零；8,888 对仍在同一顶点真为零。
- 共享 DD 评价在 103 个区域给出有限上下表达式，体积约 .027402960601490015；248 个区域仍有个别 cone 的真极点，体积约 .13926370606517666。这里的“有限/未解”是**逐 cone 表示**状态，不是将 Green 真实奇性扩大为这些区域。
- 103 个区域的非负下界全部退化为 0；其范围积分上界表达式约 `9.11779508732045e23`。这不是实际积分值，也不是反例。靠近被切去 corner 的宽多面体中，个别分母的极小值被独立相乘，再与其他 cone 相消，导致巨大过估计。
- 可由已审径向包络给全积分的有限但无用包围 `0 <= I4 < 310004`（更紧的上界表达式约 310003.20233705469，准确上端在 JSON）。因此本批没有达到 .00019，也没有发现 I4>.00019 的反例。

要闭合目标，仍需在 248 区域上的**全 cone 联合取消**，以及对其余区域中强烈相消、宽范围依赖的连续上包围。把这些 node-zero 边界删成零测集合不解决邻域积分预算；把它们算成真实方向线帽则错误。另一条主负责的正 time cutoff / 全 face Lq 预算接口已收到，但没有在本批改路线、重复其常数枚举或冒称完成其 H-epsilon 连续积分。

## 精确覆盖及体积账

G 对 signed permutations 不变，故正 octant 被 6 个排序区域覆盖，交叠仅为零测 ties；取

\[
 T=\{0\le w\le z\le y\le1\},\qquad |T|=1/6.
\]

240 条三角化 rays 中，恰有 30 个规范符号的 node hyperplanes 穿过 T 内部。初始切掉

\[
 C_\epsilon=\{(y,z,w)\in T:(1-y)+(1-z)+(1-w)<\epsilon\}.
\]

这只是 corner 三维 simplex 的一个排序片，准确体积 epsilon^3/36；它包含于以 (1,1,1) 为中心、半径 epsilon 的帽，不把 simplex 当球。其余域按 30 个超平面逐个准确二分，区域数为

`2,3,4,6,8,9,10,11,19,21,27,29,34,41,48,55,64,70,87,105,119,135,144,165,180,198,229,255,309,351`。

最终体积准确为

\[
 \frac{5999999999999999999999999999999}{36000000000000000000000000000000}
 =\frac16-\frac{\epsilon^3}{36}.
\]

每次 convex polytope 半空间切分保留该侧旧顶点及穿越边交点。边由共享的两个独立 supporting plane 识别；交点为 Fraction。体积 checker 对每个 incident supporting plane 收集 facet 顶点，经准确二维 Andrew hull 排序、三角化，并累加内部质心 pyramid 的准确行列式。它逐 facet 验证全体顶点位于同侧，并核总覆盖体积等于上式。分割的两侧内部不交，故这里同时有覆盖和体积账，非采样命中检验。

## 单 octant 帽与全部 faces 的换算

在 face x1=1 上，八条 body 线交于 q=(1,±1,±1,±1)。在 [0,1]^3 内最近交点总可取 q+=(1,1,1,1)，记 r=|(y,z,w)-(1,1,1)|。对任意交点 q，沿 face 的 z=x-q 满足

\[
 \operatorname{dist}(x,\operatorname{span}q)^2
 =|z|^2-(q\cdot z)^2/4\ge |z|^2/4.
\]

因此 d(x,S)>=r/2；纯 body 没有平行此 face 的方向线，所以此结论不需要先限 r<=1。030 已审热核估计可取

\[
 D=1/(48\pi^4),\quad B=180224/\pi^3,
 \qquad G(x)\le D+B/r.
\]

由 (u+v)^(4/3)<=2^(1/3)(u^(4/3)+v^(4/3))，整个三维球半径 a 的预算为

\[
 E_4(a)=2^{1/3}4\pi\{D^{4/3}a^3/3+(3/5)B^{4/3}a^{5/3}\}.
\]

corner 球在正 face octant 中只占一个球 octant，故 I4 的帽上界 E4(epsilon)/8；排序片再占其 1/6，故在 T 内上界 E4/48。全部 cube boundary 是 8 faces ×8 signed octants，所以总 face 积分为 64 I4，精细帽预算是 8 E4。030 的通用 64 E4 完整球预算依然有效，只是更粗；不能把“64 个 corner”与每 corner 的一个球 octant 再漏算或重算。

256bit Arb 实际结果：

- E4(10^-10)/8 = [2.6736125294051134082674553039069986238606867683181986603865e-12 +/- 1.28e-87]；
- E4(10^-10)/48 = [4.4560208823418556804457588398449977064344779471969977673109e-13 +/- 2.13e-88]；
- 全 faces 精细帽 8E4 = [1.7111120188192725812911713945004791192708395317236471426474e-10 +/- 8.16e-86]。

准确有理端点在 `ci4_caps_results.json`。同一纯 body 径向包络在整个正 octant 成立，而它包含于 q+ 半径 sqrt(3) 球的一个 octant，故 I4<=E4(sqrt(3))/8。这给前述无用但真实的全域上界，不是将局部帽预算外推为数值猜测。

## 表示失败的准确实例与依赖

`ci4_correlation.py` 实际调用全节点 residue，重现

\[
 x=(1,1,1,0),\quad d(x,S)^2=3/4,\quad G(x)=17/(640\pi^3).
\]

此点不在真实方向线上。cone `(60,14,61,20)` 的 ray14=(-1,1,0,0)、ray61=(-1,1,0,1) 却同时给 node0。全部 cones 合并后的 zero poles 准确为空。这是逐 cone 不可直接包围、但全 G 有限的明确例子；没有按零节点猜 sign(0)。

数学表示依赖主 `upper/u11_green.py` SHA fdfa422f2105710c84e7e43fdff1a4dcbef20dcc4e8a69d5a42116db3010b99d 和 `upper/u11_representation.md` SHA 63f759e850f3864390cf2af6bf39c57ead78130421389a2cfc5e853335063ea5；完整 U11 表示独审仍 pending。帽引理依赖主 `proofs/round011_face_caps.md` 与 030 已审版本。

连续 DD 代数未在本目录重写，调用 upper_theory 的 `u11_minimum_dd.py` SHA 2ed3070118111b93b7f5091e8d66145d42a49ac6297ba4f4974e148cf6743797。它以单个正/负节点的 1/3 公式和 2/2 公式分别给准确单侧界，并对每个 p+n 分母及 2/2 numerator 取同一 vertices 的 affine hull。Python Fraction 内核共享；最后 pi^-3 及正 cube roots 用 Arb 包围。共享来源与算法依赖已披露，本批不声称独立重证该 evaluator。

## 复现与停止

在本 clone 运行（都只写本目录已登记结果，caps/integral 需显式 --save）：

```sh
PY=/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python
D=cube-finite-dim-bounds/upper/continuous_integral_n4
$PY $D/ci4_partition.py --save $D/ci4_partition_results.json
$PY $D/ci4_correlation.py
$PY $D/ci4_caps.py --save $D/ci4_caps_results.json
$PY $D/ci4_integral.py --save $D/ci4_integral_results.json
```

原始 Green、030帽来源在主 clone 只读；共享 DD 在 upper clone 只读，代码内固定路径和 SHA。移入主目录后应将该依赖路径改为主整合的同 SHA 文件后另记 checker hash。运行中修正过一次初诊 pair 朝向统计 predicate（16920 为一半计数，最终33840）及 Arb.root(0) 返回 nan 的实现问题（明确将 0^(4/3) 取0）；这两处未改变分区或增加研究修订。

本批没有触达任何 STATE/总账、其他任务源码、Git staging/commit/push。未成功，故没有向 independent_audit 宣称可认证目标或锁定新 U4。停止机械加密；后续若继续，应先有跨 cone 连续取消或另一完整连续积分机制，再另行登记。
