# Near-flat cube 几何：对称赢家切换、后验刚性与嵌套碰撞约束

2026-10-07，nearflat_cube_geometry，父任务指定 6.1-sol high。只新增本稿及同前缀后处理文件，不改主 TXT 或他人稿。使用 M03 的亏损/量词审计和 M02 的支撑证书边界；已读两 skill 的 SKILL、method、provenance、cube-interface。以下证明自含，不引用未核外部稳定性定理。

**取得的一般引理：双向扰动不仅约束 source column，还约束 receiver 上不同候选尺度的完整来源后验差。精确乘性局部极大点的真赢家在超阈集上几乎处处唯一；近全局极大点满足显式 hinge-regret 约束。嵌套 cube 将后验差变成真实捕获后验的碰撞量。尚未由此控制 K 到 sqrt(n) polylog(n)。**

已读 log_max_variation、其 independent_review、source_occupation_tail_geometry 和 winner_envelope_entropy_bridge。本轮全文查重本目录及 cube_general_20261003 下的 winner uniqueness / posterior tie rigidity / hinge 关键词；旧稿有 source 平坦、双接收点交集身份、posterior 非产品及一般 hinge 工具，没有检出本稿这个双向 log-max 切换与嵌套 posterior-collision 公式。不将来源一阶条件、K 与 weak 常数的等价或 entropy 身份重新计为进展。父任务的全空间 additive obstacle 和 nearmax 弱障碍界由其它工作者负责，本文不重复证明。

## 1. 原完整输入及 receiver 后验

固定完整有限正 Borel 来源 μ，W>0，全部原边长 a≤R≤b≤2a，或固定共同有限原尺度目录。所有 receiver 积分仍为 Lebesgue。令

\[
U_R(x)=h_R*\mu(x),\quad M(x)=\sup_RU_R(x),\quad
E=\{M>\tau\},\quad I=\tau|E|,\quad Z=1+n\log(b/a).
\]

保留原可测真赢家 R(x)。连续闭 cube 窗口的最大值/最小赢家可测性采用已审通过的 log-max 原稿；finite family 不添加任何新尺度。对 x∈E，另选任意可测原合法候选 L(x)，要求 U_L(x)>0，定义

\[
q(x)=U_L(x)/M(x)\in(0,1],\quad \delta(x)=-\log q(x),
\]
\[
\pi_{R,x}(dy)=\frac{\mathbf1_{Q(x,R)}(y)\mu(dy)}{\mu(Q(x,R))},\quad
\pi_{L,x}(dy)=\frac{\mathbf1_{Q(x,L)}(y)\mu(dy)}{\mu(Q(x,L))}.
\tag{1}
\]

这些是由同一完整 μ 产生的真实捕获后验，不是均匀 cube 来源；没有坐标独立。若某候选无正响应，就在该点取 L=R，不能定义零分母后验。原 Φτ=τ∫log_+(M/τ)，全输入上确界 K=sup Φτ/W≤Z。

## 2. 对称切换下界：一阶 kink 减二阶 log 曲率

任取空间上的 μ-可测 |v|≤1。令 μ±=(1±tv)μ，0<t<1，p_R=π_Rv，p_L=π_Lv；保持完整来源的所有物理位置与原合法尺度。原/候选响应对新输入分别为

\[
M(1\pm tp_R),\qquad qM(1\pm tp_L).
\]

新 true max 至少取这两个合法 query 的较大值。对 E 上原 log 和新 log_+ 相减，利用 log_+≥log，且 E 外原项为零而新项非负。Taylor 下界
log(1±tp)≥±tp−t²/[2(1−t)²] 对 |p|≤1 有效。于是两个方向相加给

\[
\boxed{
\Phi_\tau(\mu_+)+\Phi_\tau(\mu_-)-2\Phi_\tau(\mu)
\ge \tau\int_E
\big(t|\pi_{R,x}v-\pi_{L,x}v|-\delta(x)\big)_+dx
-\frac{t^2I}{(1-t)^2}.}
\tag{2}
\]

准确的点态代数为
max(a,b−δ)+max(−a,−b−δ)=(|a−b|−δ)_+。
不对上确界求导，不假设 winner 唯一或 threshold plateau 零测，不需要选出的尺度对新 μ 保持获胜。可积性由 I≤ZW、|p_R|,|p_L|≤1 及原窗口包络保证；(2) 中不要求 ∫δ 有限。

若 μ 是全局 ε-near-max，即 Φτ(μ)≥(K−ε)W，则新质量 W±=W±t∫v dμ，故 Φτ(μ±)≤KW±。质量增减在相加后准确抵消，得到

\[
\boxed{
\tau\int_E\big(t|\pi_Rv-\pi_Lv|-\delta\big)_+dx
\le 2\epsilon W+\frac{t^2I}{(1-t)^2}.}
\tag{3}
\]

该量词是对每个同一空间 source test v、每个原合法可测候选 L 都成立。不是允许 v=v(x,y) 随 receiver 另选；那会绕过同一个来源扰动。

若 L 也是真赢家，δ=0，记 m=I/W≤Z，则

\[
\tau\int_E|\pi_Rv-\pi_Lv|dx\le B_\epsilon(t)W,
\quad B_\epsilon(t)=2\epsilon/t+t m/(1-t)^2.
\tag{4}
\]

例如 0<ε≤Z/4，取 t=√(ε/(4Z))≤1/4，粗放松给 Bε≤6√(εZ)。更直接保留 m 的 (4) 可减小常数。近赢家用 (3)，其 δ 是真实全输入 response gap，不是免费宣布某个网格覆盖全部 nearwinner。

## 3. 精确乘性局部极大点：超阈真赢家唯一

不要求极大点存在。若某 μ 确实对全部保持质量的乘性扰动是 Φτ 的局部极大点，则 E 上所有合法真赢家几乎处处相同。

证明：固定两个可测真赢家 R,L，以及某个 |v|≤1、∫v dμ=0。局部极大允许步长邻域依赖 v；在该邻域内 (2) 左侧≤0，δ=0，故
τ∫E|π_Rv−π_Lv|≤tI/(1−t)²。令 t↓0，积分为零。

取一个可数确定概率测度的半开有理盒 π-system；对每盒 B 使用
v_B=1_B−μ(B)/W，它保持质量且 |v_B|≤1。两个后验对常数都积分为1，故 π_Rv_B−π_Lv_B=π_R(B)−π_L(B)。逐盒删去 receiver 零集后，两个概率后验在全部确定集上相同，因此在同一个 full-measure receiver 集上 π_R=π_L。

对同中心嵌套 cube，若 R<L，真赢家并列给

\[
\frac{\mu(Q_R)}{\mu(Q_L)}=(R/L)^n<1,
\quad \pi_R(Q_R)=1,\quad\pi_L(Q_R)=(R/L)^n.
\tag{5}
\]

与两个后验相同矛盾。连续 compact 窗口可取最小/最大真赢家；最大赢家的 {Rmax≥s}={M_[s,b]=M} 同样 Borel，所以若存在多个赢家，已被这对 selector 捕获。finite family 则取目录中的最小/最大 index。该证明保留原闭 cube 及边界原子，不将不同物理位置合成 free posterior。

这是必要的 receiver 刚性，比仅有 source S 常数更强；并不保证一个唯一赢家输入是极大点，也不保证全尺度的极值存在。

## 4. 任意有限来源 partition 的后验差预算

取任意固定有限 μ-可测 partition C1,…,CN，置
a_i(x)=π_R,x(C_i)−π_L,x(C_i)。对同一来源 partition 的辅助独立 Rademacher 符号 ξ_i，令 vξ(y)=ξ_i 于 C_i。此随机化只选择合法 source test；不假设 Y 的来源坐标独立，不按 receiver 重选符号。

对 A=Σa_iξ_i，直接展开给 E A²=Σa_i²、E A⁴≤3(Σa_i²)²。Hölder 从 E A²≤(E|A|)^(2/3)(E A⁴)^(1/3) 得

\[
\mathbb E_\xi|A|\ge\frac1{\sqrt3}
\Big(\sum_i a_i^2\Big)^{1/2}.
\]

对 (3) 求辅助随机期望，用 hinge 的凸性与正 Tonelli，得到一般完整来源式

\[
\boxed{\tau\int_E
\left(\frac t{\sqrt3}\|\pi_R(C_i)-\pi_L(C_i)\|_{\ell^2}-\delta\right)_+dx
\le2\epsilon W+\frac{t^2I}{(1-t)^2}.}
\tag{6}
\]

真赢家 δ=0 时，特别有

\[
\tau\int_E\|\pi_R(C_i)-\pi_L(C_i)\|_{\ell^2}dx
\le\sqrt3 B_\epsilon(t)W.
\tag{7}
\]

无精确/近极大前提的输入只享有 (2)，不能使用 (3)、(6)、(7) 的 ε费用。

## 5. 嵌套 cube 把差预算转成真实 source collision

一般 finite partition 公式 (6) 对所有 μ 有效。以下 collision 简化针对完整有限正原子源；首先合并物理上 coincident 的 atom，不能把两个相同 y 的标签赋予不同独立 source sign。partition 取这些不同物理 atom。

设 R(x)≤L(x)，R 是原真赢家，L 是原正响应候选。写 θ=n log(L/R)，q=e^-δ，β=μ(Q_R)/μ(Q_L)=e^(δ−θ)≤1。完整嵌套来源后验准确满足

\[
\pi_L(y)=\beta\pi_R(y)\ (y\in Q_R),\qquad
\chi_R(x)=\sum_{y\in Q_R}\pi_R(y)^2,
\]
\[
\|\pi_R-\pi_L\|_{\ell^2}\ge(1-\beta)\sqrt{\chi_R}.
\tag{8}
\]

于是 (6) 给真实 receiver-overlap 约束

\[
\boxed{\tau\int_E
\left[\frac t{\sqrt3}(1-e^{\delta-\theta})\sqrt{\chi_R}-\delta\right]_+dx
\le2\epsilon W+\frac{t^2I}{(1-t)^2}.}
\tag{9}
\]

对两真赢家 δ=0，(7)–(8) 给

\[
\tau\int_E(1-e^{-\theta})\sqrt{\chi_R}dx
\le\sqrt3 B_\epsilon(t)W.
\tag{10}
\]

若某子集 A 上 θ≥c>0 且 Q_R 捕获的不同物理 atom 数≤N，则 χ_R≥1/N，从而

\[
\tau|A|\le\frac{\sqrt{3N}}{1-e^{-c}} B_\epsilon(t)W.
\tag{11}
\]

这是完整原有限源、Lebesgue receiver、所有原尺度都参加 M 后的约束：大尺度分离的并列 plateau 在近极大点只能发生于足够弥散的捕获后验，或占极小体积。(9) 同样保留非零真实 response gap，可避免只考虑精确 ties 的不稳定性。

**没有把 (10) 当作一般 source-square 上界。** χ_R 是 receiver 捕获后验对物理 atom 的 collision；它不同于 ∫S²dμ，也不同于两个 receiver 的 μ(Qx∩Qx')/m_xm_x'。一般 nonatomic μ 上细化 partition 的 ℓ²差可以趋零，不能免费把 (8) 的 atom collision 或1/N用于全部输入。

## 6. 已执行的新三轮保存-cell 后处理

详见 nearflat_cube_geometry_probe_20261007_registration/results/receipt 与该 prefix 的三个 posterior_cells gzip。父任务明确协调此 probe：仅使用既有 log_max_variation_probe_20261007_n1/n2/n3_cell_profiles.json.gz 的完整 cell、原响应、六个 ±t 扰动响应；无 oracle、arrival 或 MC 重跑。新阈值 τ'=原τ/2 在新执行前登记，原结果和原阈值未修改。每 cell 取原同 max 的最小/最大尺度；无 tie 时两者相同。

|n|保存 cell 数|新 E cell 数|新 E 的 exact-tie cell 数|新 tie 体积|I=τ'\|E\||A=τ'∫E\|δp\||
|---|---:|---:|---:|---:|---:|---:|
|1|10|7|2|1|1/2|1/16|
|2|256|69|12|3/4|17/32|9/256|
|3|10648|1775|32|5/16|1343/2048|35/4096|

全部阈值、winner、capture posterior、β、collision 与 RHS 用 Fraction。三 positive t 固定为1/8、1/32、1/128。物理 coincident 来源按 signed mass / mass 合并，三个输入的 v均 |v|≤1且∫vμ=0。保存的 posterior-min/max、δp 与 collision 可逐 cell 重建，没有来源采样。

9/9 对称 switch 下界严格通过；9/9 实际 pair Φτ' 增量由独立有理 log 区间严格认证为正。解析 RHS=tA−t²I/(1−t)² 本身正的有5/9：n1/n2在1/32、1/128正；n3只在1/128正。负 RHS仍是有效下界，不因实际 pair增量正就将其改成正。

log区间来自 dyadic range reduction 后 log r=2Σz^(2j+1)/(2j+1)，z=(r−1)/(r+1)≤1/3；80项截断余项被 2z^161/[161(1−z²)] 支配。所有区间端点为确切有理数，符号相减使用 lower/upper 正确方向。160位 Decimal另作诊断，未当作 interval certificate。程序初次 n1/n2成功，n3仅在把巨大有理端点转字符串时触及 Python4300digits 上限；已单独保存 repair receipt，提高 serialization ceiling后重复保存-cell算术，已存 n1/n2 payload逐项相等复用，没有改变数学输入/参数/判据。

三个 A都严格正，因此对当前固定 v，(2) 给所有足够小 positive t 的 pair增量严格正。例如 t<min(1/2,A/(4I)) 时 tA−t²I/(1−t)²>0。由于两扰动均保持 W，至少一个方向提高 Φτ'；故这三个输入在新 τ' 下**严格不是乘性 localmax**。这条排除依据解析 A>0与精确代数，不把有限3步长检查冒充任意小邻域论证。

这仍是 finite-family 工具验证，不认证 continuous [a,2a] 的 oracle、不声称三个输入 nearoptimal、不确定 K，也不给维数阶数证据。

## 7. 距离 sqrt(n) 几何付款的准确缺口

(3)/(6)/(9) 是 nearflat之外可用于真几何的更强必要条件，清楚保留“一个全局 source test，receiver 选尺度”的量词差别。它排除具有真实可区分并列/近并列后验的大体积区域，而不先要求全输入 source平方或 pairGaussian。

剩余不是再证明 S≈未知K；需要证明大 K 的近极大输入必然产生足够多、足够可由同一 source test/固定 partition 区分的原近赢家，且其 (6) 左边有 dimension-effective 下界。当前没有这种强制 nearwinner overlap 的几何结论。单纯 R/L 的分离不保证 q≈1；单纯 nearwinner多也不保证后验 ℓ²差不小。一般来源可以弥散，source partition不具维数无关有效 atom count。不能据有限原子验证套 (11) 到 arbitrary μ，不能把 continuum细化损失隐藏在 n^{o(1)}里。

即使以后取得这个有限窗口费用，全尺度迁移和 actual FIRST/history仍需原合同；本稿保留所有原窗口尺度，没有每尺度另领完整 W，也没有免费全尺度结论。当前目标 sqrt(n)n^{o(1)} weak(1,1)仍未证明。
