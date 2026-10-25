# 第六轮 E：D Astra 二聚化的独立有限 L¹ 接口审计

结论：所述二聚化源、共同半径及 `16*union` 归一化成立。active 长度扩到 8 后，可统一使用 **active 粗周期 collar 3、其余五轴 collar 2**。若连未通过交集筛选的起点枚举块都保留，active collar 4 更保守。下面给出严格范围证明，不继承旧的 ell≤4 collar。

D 的 `dimer_results.json` 状态为 `completed_no_promotion`，24 个候选均无晋级；E 只读取状态、范围与源码，没有重新计算其数值并集，没有准确候选体积认证，也没有新 C6 界。本件不把有限批次无晋级解释为整个机制的排除。

## 1. 已核对的模型与 spacing scaling

以冻结 R6-003 的正 spacing a 为基准，令 0≤ε<1、|δ|<1，s⁶=1+ε 且 s>0。选 active 轴 i，

\[
b_i=a_i s/(1+\varepsilon)=a_i(1+\varepsilon)^{-5/6},\qquad
b_h=a_hs\quad(h\ne i).
\]

因此

\[
\prod_{h=1}^6b_h=\prod_{h=1}^6a_h\,s^6/(1+\varepsilon)=\prod_h a_h.
\]

active 位置与质量是

\[
P_i(k)=b_i(k+\varepsilon(k\bmod2)),\qquad
w(k)=1+\delta(-1)^k.
\]

其余轴为 b_h Z。active 平移 k↦k+2 对应物理平移 2b_i，权重不变；其余轴周期为 b_h。一个完整粗周期含两个指标原子，质量 `(1+δ)+(1−δ)=2`，物理体积为 2∏b=2∏a。因此这条路径的密度固定。即使不固定密度，下述 `16*union` 公式对任意正 b 仍成立；spacing scaling 影响的是物理接收盒体积。

由于 active 两个相邻间隙分别是 b_i(1+ε)、b_i(1−ε)，0≤ε≤1 时 active 位置单调；ε=1 允许端点重合，仍可按指标原子的两个非负质量相加。此次筛选只用 ε<1。

## 2. 只反射五个 inactive 轴

采用

\[
H_i=[0,2b_i],\qquad H_h=[0,b_h/2]\ (h\ne i).
\]

只对 inactive 坐标反射 k_h↦−k_h。它们保持 active 位置与权重不变，产生 32 个内部互不相交的 H 反射像，覆盖完整粗周期

\[
[0,2b_i]\times\prod_{h\ne i}[-b_h/2,b_h/2].
\]

对于 U⊂int H，若指定 cube 平均至少 α，则周期形式比值准确为

\[
\alpha\frac{32|U|}{2}=16\alpha|U|.
\]

α=1 即 D 的 `16*union`。同时 ε≠0、δ≠0 时，不能从旧格点模型继承 active 轴反射；本件完全没有使用该反射。无限周期源只说明几何，第 6 节直接构造有限 L¹ 输入。

## 3. 原子个数、质量、极值与共同半径

给定六维笛卡尔指标块 B(q,ℓ)，active 长度 ℓ_i∈1..8，inactive 长度 ℓ_h∈1..4。它有 ∏ℓ_h 个指标原子，至多 8·4⁵=8192 个。active 连续段质量是

\[
m_i=\begin{cases}
\ell_i+\delta(-1)^{q_i},&\ell_i\text{ odd},\\
\ell_i,&\ell_i\text{ even}.
\end{cases}
\]

总质量为 M=m_i∏_{h≠i}ℓ_h。active span 为

\[
b_i\begin{cases}
\ell_i-1,&\ell_i\text{ odd},\\
\ell_i-1+\varepsilon(-1)^{q_i},&\ell_i\text{ even}.
\end{cases}
\]

由单调性，active 最小/最大位置就是 P_i(q_i)、P_i(q_i+ℓ_i−1)。inactive 极值为 q_hb_h、(q_h+ℓ_h−1)b_h。若记极值为 L_h≤U_h，则中心盒条件为

\[
x_h\in[U_h-r,L_h+r]\quad\text{对全部六轴},\qquad r=M^{1/6}/2.
\]

这是一个所有坐标共用的半径，cube 体积准确为 M。cube 包含整个选定块，其他原子的非负质量只会增加平均值。未来准确认证时可选根的有理下端 r₋≤r，再严格内缩接收盒并记录到各 cube 面的最小余量 μ>0。这里没有以六个不同半径替代共同半径。

## 4. active 长度 8 的统一保留范围

下面在更宽的闭参数范围 0≤ε≤1、|δ|≤1 上给保守界；零质量块丢弃。

冻结 spacing 全部 a_h>1，E 用纯 Fraction 六次根包围实际核验了这一点。因为 (1+ε)⁵≤2⁵<2⁶，

\[
b_i=a_i/(1+\varepsilon)^{5/6}>a_i/2>1/2.
\]

active 偶数长度至多 8；奇数长度至多 7，因此 m_i≤8。于是 M≤8192，并且

\[
r\le8192^{1/6}/2<9/4,
\quad\text{因为 }(9/2)^6>8192.
\]

考虑任何与 H 有非空接收交集的 active 块，末指标记 v=q_i+ℓ_i−1。某个中心 x_i∈[0,2b_i] 的 cube 包含两端点，必有

\[
P_i(q_i)\ge-r,\qquad P_i(v)\le2b_i+r.
\]

又有 0≤ε(k mod2)≤1、r/b_i<9/2，故

\[
q_i\ge-r/b_i-\varepsilon(q_i\bmod2)>-11/2,
\qquad v\le2+r/b_i<13/2.
\]

两者都是整数，得到 **q_i≥−5，v≤6**。整个连续指标块均在 [−5,6]。注意此处没有使用 active 反射；若错误继承旧 ell≤4 的指标 [−2,2]，不能得到本结论。

inactive 采用 D 源码的中央起点：奇数 ℓ=2k+1 时 q=−k；偶数 ℓ=2k 时 q=1−k。对 ℓ∈1..4，它们及 inactive 反射后的全部指标落在 [−2,2]。

## 5. 起点枚举 supersets 与一般 collar 公式

D 的 active 起点循环使用 r_*=[P(ℓ_i+δ)]^(1/6)/2，其中 P=∏inactive ℓ_h，并加了一格缓冲：

\[
q_{\min}=\lfloor-r_*/b_i-\varepsilon\rfloor-1,\quad
q_{\max}=\lceil2+r_*/b_i-(\ell_i-1)\rceil+1.
\]

在实际 δ≥0、δ≤1 范围，r_*<(5/2)，因为其六次幂的质量至多 9·4⁵<(5)⁶。因此 r_*/b_i<5，得到 q_min≥−7、q_max+ℓ_i−1≤8。该循环所有块的指标都位于 [−7,8]；数学上起点条件必要的块都位于该循环 superset。浮点实现仍只作发现工具，本件不把浮点筛选误差变成准确排除。

对任意有限块族，更直接的准确公式是不依赖长度 cap 的

\[
p_i=2,\ p_h=1\ (h\ne i),\qquad
S_h=\max_B\max(|q_h|,|q_h+\ell_h-1|),\quad
C_h=\lceil S_h/p_h\rceil.
\]

因此，本轮**筛选后非空接收盒**统一可取 C_i=3、C_h=2 (h≠i)；若全部起点循环 supersets 都纳入，取 C_i=4、C_h=2。已有具体块数据时，可按一般公式获得更小且仍充分的 collar。

## 6. 有限源、准确边界损失与一般 L¹

接收粗周期指标 t∈[−K,K]⁶。所需源指标在 active 轴是 2t_i+s_i，在 inactive 轴是 t_h+s_h（含允许的 inactive 反射）。保留完整源指标粗周期

\[
t'_h\in[-K-C_h,K+C_h],\quad
b_h'\in\{0,\ldots,p_h-1\},\quad k_h=p_ht'_h+b_h'.
\]

这里 b_h' 是剩余类指标，勿与 spacing b_h 混淆。所保留的指标区间为

\[
[-p_h(K+C_h),\ p_h(K+C_h+1)-1].
\]

p_hC_h≥S_h 保证所有所需原子保留。active 位置仍使用原定义 P_i(k_i)，不是把二聚化原子重新放回均匀格点。每完整粗周期权重和为 2，故有限源准确总质量为

\[
2\prod_h(2K+2C_h+1).
\]

取指定接收盒到其共同 cube 面的最小严格余量 μ>0，η=μ/2，定义单位质量 bump b_η=(2η)^{-6}1_{[-η,η]^6}，在每个保留源位置按其权重平移相加得到 f。此 f 是有限非负和，属于 Lebesgue L¹；其范数准确等于上式，即使某些 bump 重叠也是如此。每个指定 cube 包含所需原子的整个 bump，故所有接收点 Mf≥α。

若未来独立内缩并集给 R₋=16α|U|>T（T 是冻结完整上端），则取 β=α(R₋+T)/(2R₋)<α。接收点均 Mf>β，接收体积为 32|U|(2K+1)⁶，其有限弱型比值下界是

\[
\frac{\beta}{\alpha}R_-\prod_h\frac{2K+1}{2K+2C_h+1}.
\]

有限 K 足够大时该式严格大于 T；可用 Fraction 倍增 K 准确检验。全部正尺度最大算子包含这些指定共同半径，故传递直接适用于原六维全尺度一般 L¹ 问题。这只是未来晋级候选的严格接口；本轮没有 R₋>T 的新候选。

## 7. 实际独立运行

E 只运行 `audit_interface.py`。它没有导入 D 的任何程序，没有调用 union 函数，没有重跑 baseline 或 coalesced endpoint 的体积。实际检查 **96 个 active 连续原子区间**（每区间 1..8 原子），覆盖三种 active 轴、两种 ε、两种 δ，包含负起点、全部 active 长度 1..8、全部 inactive 长度 1..4。逐区间以 Fraction 核对质量、端点、span、周期、r⁶ 包围、spacing scaling 根包围及必要的 [−5,6] 范围。另做 6336 个端点 collar 比较和三个粗周期准确质量比较。墙钟上限 600 秒，实际运行约 0.007 秒，`pass`。

有限损失函数另做一个**条件算术例**：假设 R=1、T=1/2、μ=1/100、C=(3,2,2,2,2,2)，准确找到 K=32、β=3/4，finite ratio−T=1388073439/148061634372>0。此例没有提供一个新的几何 witness，不是新 C6 候选。

完整输入哈希、源码哈希、范围及实际结果见 `interface_audit_results.json`。复现：

```bash
PYTHONDONTWRITEBYTECODE=1 python3 \
  cube-dim-6-bounds/audit/lower006/dimerization/audit_interface.py
```

源码仅复用 E 已有 Fraction 根/spacing 与有限损失辅助函数；输入路径指向本轮 D 独立 clone 和主 clone 的只读冻结文件。所有输出仅在 `audit/lower006/dimerization/`。未使用 MATLAB，未暂存、提交或推送；完成后停止本项计算。
