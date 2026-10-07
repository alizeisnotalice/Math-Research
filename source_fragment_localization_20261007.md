# 真实高势 fragment 的组件局部化：单标签预算与重叠链边界

2026-10-07。父任务指定 6.1-sol high。本文给一份适用于每个原 L¹ 输入的固定空间分解，不增添覆盖假设：原 B={V_hi>H} 的支持集按窗口 b 的邻接组件拆为 ξ_j，任一原 receiver 的正 B 捕获只属于一个标签。因此把 global-fullness 换成 local-fullness 时没有组件个数费，来源抵扣总额仍只为 H W_B。这支付远离重复 clusters 的全局质量稀释障碍。连通长链仍能持续 fragment；本文给出两个真实有限赢家空间族，后者只在小固定 H 下否定“高势自一致性自动保证组件 fullness”，不反驳 H≈√n 的目标尾。

使用已读 [E03](</Users/zhengzhihao/.codex/skills/math-e03-tree-bellman-cross-layer-budget/SKILL.md>)、[G01](</Users/zhengzhihao/.codex/skills/math-g01-carleson-measure-box-and-tent-test/SKILL.md>) 的质量/packing 审计流程及 method/cube-interface/provenance。这里直接证明 countable source partition 与单标签求和，不套未验证树 Bellman、Carleson 嵌入或 private TC 接口。

查重发现 [kernel-core CP5](../anisotropic_cube_flow/kernel_core_peeling_next.md) 已证明 L¹ 支持邻接组件与查询单组件原理，且明确大组件可任意长。本文不把这个几何原理重新报为新成果；新增的是将它接到 [实际高势一次截断](general_tail_goodlambda_20261007.md) 的 J_full 切线预算，同时核验真实 B 下远离/重叠两类例子。旧 small-component halo receipt 不另领一次。

## 1. 量词与原对象

对每个 μ=fdy≥0、W<∞，固定原 μ_hi≤μ、合法可测原 winner R(x)∈[a,b]（a>0，b≤2a）、原 E_R⊂{M>τ}、I_R=τ|E_R|<∞。始终使用完整 μ 的

\[
 M(x)=h_{R(x)}*\mu(x),\quad m(x)=M(x)R(x)^n,
 \quad A_x=Q(x,R(x))\setminus Q(x,a).
\]

冻结

\[
 V(y)=\tau\int_{E_R}\frac{\mathbf1_{A_x}(y)}{m(x)}dx,
 \quad B=\{V>H\},\quad\xi=\mu_{\rm hi}|_B,
 \quad w=\xi(\mathbb R^n).
\]

w=0 时尾为0。以下 w>0。B 依赖真实全 profile，但不是按 receiver 重组。来源是原 ξ≪dy，未添加背景或重新选 ξ 的赢家。上一稿准确恒等为

\[
 F(H)=\int(V-H)_+d\mu_{\rm hi}
 =\tau\int_{E_R}\frac{\xi(A_x)}{m(x)}dx-Hw. \tag{1}
\]

下述分解对每个这样的 ξ 都存在，无“如果 packets 覆盖”的额外前提。

## 2. 固定可数空间组件与每 row 唯一标签

令 S=supp ξ，为有限 Radon 测度的闭支持，ξ(S)=w。取开放扩厚集

\[
 U=S+(-b/2,b/2)^n
   =\bigcup_{y\in S}\bigl(y+(-b/2,b/2)^n\bigr).
\]

U 的连通组件 U_j 是两两不交开放集。每个组件含一个有理点，故至多可数；用有理点固定枚举，无可测 argmax 问题。S⊂U，令

\[
 S_j=S\cap U_j,\qquad\xi_j=\xi|_{S_j},
 \quad w_j=\xi_j(\mathbb R^n),\qquad
 \xi=\sum_j\xi_j,\quad\sum_jw_j=w. \tag{2}
\]

非零组件的 w_j>0：它包含某个支持点的开放小邻域，而支持定义给正 ξ 质量。零项也可直接丢弃。组件可无界、非凸、任意长；不从连通性推小直径。

**查询单标签引理。** 对每个 x，ξ 不给 Q(x,R(x)) 的面质量，故正捕获都来自 S∩int Q(x,R(x))。任意两个内部来源 y,z 满足 ||y−z||∞<R(x)≤b，因此两开放盒 y+(-b/2,b/2)^n 与 z+(-b/2,b/2)^n 相交，处于同一 U 组件。于是每个原 query 的全部正 ξ 捕获属于同一 j；annular 子捕获也只有这个标签。

定义可测

\[
 c_j(x)=\xi_j(A_x),\qquad E_j=\{x\in E_R:c_j(x)>0\}.
\]

联合 cube 核可测及 Tonelli 保证 c_j 可测。各 E_j 两两不交，且在 c(x)=ξ(A_x)>0 的 row 上，c(x)=c_j(x)，j 唯一。c=0 的 row 没有 B 交通。特别地

\[
 \sum_j\tau|E_j|\le I_R. \tag{3}
\]

这不是任意固定空间网格都具有的性质。细网格可以让每个 query 完整捕获许多 cells，导致 receiver 被重复收费；当前组件分解牺牲小直径，换取真实唯一标签。

这里若 c_j(x)=0，单个来源点仍可能形式上位于 A_x，不能逐点误称没有 incidence；Tonelli 给 ∫ξ_j∫E_R\E_j1_Ax dx=∫E_R\E_j c_j dx=0，才保证下述 ξ_j-a.e. 恒等。对 ξ_j-a.e. y，所有对 V(y) 有正 dx 贡献的原 row 属于 E_j，所以原高势资格原样局部保留：

\[
 V(y)=V_j(y):=\tau\int_{E_j}
           \frac{\mathbf1_{A_x}(y)}{m(x)}dx>H,
 \qquad\int V_jd\xi_j=\tau\int_{E_j}\frac{c_j}{m}dx.
 \tag{4}
\]

没有更换完整分母或沿 source 分量重新定义 winner。由 (4) 得 H w_j<τ|E_j|；但这只是一份局部 first-moment 关系，不控制组件直径或各 query 的 c_j/w_j。

## 3. local-fullness 的一般无组件数费预算

固定 δ₀∈(0,1)、η∈(0,1]。在每 E_j 上按原 b_B=c_j/m 切

\[
 D_{0,j}=\{b_B\le\delta_0\},\quad
 D_{F,j}=\{b_B>\delta_0,\ c_j\ge\eta w_j\},\quad
 D_{G,j}=\{b_B>\delta_0,\ c_j<\eta w_j\}.
\]

令 I_{F,j}=τ|D_F,j|，J_{F,j}=τ∫D_F,j c_j/m，J_G,loc=Σ_j τ∫D_G,j c_j/m。低 share 支由 (3) 支付 Σ_j J_0,j≤δ₀ I_R。

在 D_F,j 对原 incidence y∈A_x 的捕获，m≥ηw_j 与 m>τR^n 给

\[
 \frac\tau m\le\min\left\{\frac\tau{\eta w_j},
                       \frac1{(2\|x-y\|_\infty)^n}\right\}.
\]

上一稿的完整层饼重排与 ℓ(r)=r (r≤1)、1+log r (r≥1) 原样给

\[
 J_{F,j}\le w_j\ell\left(\frac{I_{F,j}}{\eta w_j}\right)
 \le\frac H2 w_j+\eta^{-1}e^{-H/2}I_{F,j}. \tag{5}
\]

使用 (2)、(3) 非负求和，只付一次来源与 receiver：

\[
 \sum_jJ_{F,j}\le\frac H2w
                    +\eta^{-1}e^{-H/2}I_R.
\]

因此一般严格合同为

\[
 \boxed{F(H)\le
  [\delta_0+\eta^{-1}e^{-H/2}]I_R
           +J_{G,\rm loc}-\frac H2w.} \tag{6}
\]

这是从 (1) 重新划分的一份替代预算，不是在旧 global-fullness receipt 之上再相加。每个来源的负抵扣总额始终 H，不能旧拆分先用 H/2 再另给每组件 H。由于 w_j≤w，global-fullness 必为 local-fullness；故 local fragment 是上一稿 global fragment 的子集，J_G,loc≤J_G,global。新分解确实删除了“只因别处存在远离来源质量而 fullness 失败”的 rows。

在 local fragment 内还严格有

\[
 m(x)<\frac{\eta w_j}{\delta_0},\qquad
 R(x)^n<\frac{\eta w_j}{\delta_0\tau}. \tag{7}
\]

若 ηw_j/δ₀≤τa^n，则这个组件的 D_G,j 为空。一般只剩 w_j>δ₀τa^n/η 的较大组件，数量至多 ηw/(δ₀τa^n)。这依输入质量与尺度，不是无维数的费用常数；允许一个组件非常大。

式 (6) 与 (7) 已证，没有假定这些组件普遍被 full 捕获。还未证的是

\[
 J_{G,\rm loc}-\frac H2w\le\delta_1 I_R. \tag{8}
\]

下面两个真实空间族分别检验 (6) 删除什么、保留什么。

## 4. 远离重复 clusters：global fragment 全在，local fragment 全消失

固定 n=1、允许尺度 {1,2}、a=1、b=2、τ=1/4、d=1/8。对任意整数 J≥1，完整 L¹ 输入

\[
 f_J(y)=d^{-1}\sum_{j=0}^{J-1}
                 \mathbf1_{[4j-d/2,4j+d/2]}(y),\qquad W=J,
 \quad\mu_{\rm hi}=\mu.
\]

保留原 receiver 子集

\[
 E_R=\bigcup_{j=0}^{J-1}
 \{x:9/16<|x-4j|<15/16\}. \tag{9}
\]

在这些原 rows 上，Q_1 不捕获任何来源，Q_2 完整捕获恰第 j 个 packet；其它 packets 距离足够远。完整 averages 是0、1/2，所以 R=2 是真实唯一 finite winner，M=1/2>τ，m=1。也有 A_a f=0，保留任意 K 的 anchor-cold。没有添加背景。

annular 条件对整个来源 packet 恒真，故每个原来源 y 的

\[
 V(y)=\tau\,2(15/16-9/16)=3/16.
\]

取 H=1/8，则 B 在 μ-a.e. 来源上覆盖全部 μ，w=J。扩厚 U 的组件正是 J 个分离 intervals，每个 w_j=1。取 η=1/4、δ₀=1/8，J>4 时每 row 的 c=1<ηw，但 c=1≥ηw_j；b_B=1。global-fullness 全失败，local-fullness 全成立，local fragment 为空。

该例的真实尾与 receiver 质量是

\[
 I_R=3J/16,\qquad F(H)=J/16=I_R-Hw.
\]

这里 H 很小，(6) 的指数系数并不小；例子检验来源标签/预算结构，不用于声称尾吸收。它说明 global μ(B) 不收缩并不意味着需要另付 J 份 W：组件局部化在任意 J 上仍仅支付 Σ_j w_j=w。对任意具有真实高势 B 的紧支撑 base，可在查询互不影响的远离翻译 copies 上使用同一原理；不能把各 copy 当额外独立主输入费用。

## 5. 重叠链：真实 B 高势自一致，组件仍大且 fragment 持续

仍用上述 n=1、τ、d 与尺度 {1,2}，改完整 centers 为 z_j=3j/2：

\[
 f_J(y)=d^{-1}\sum_{j=0}^{J-1}
           \mathbf1_{[z_j-d/2,z_j+d/2]}(y),\qquad W=J.
\]

保留原 receiver 子集为相邻来源中点的小 intervals

\[
 E_R=\bigcup_{j=0}^{J-2}
 \left(\frac{z_j+z_{j+1}}2-\frac1{16},
       \frac{z_j+z_{j+1}}2+\frac1{16}\right). \tag{10}
\]

每 row 的 Q_1 不捕获来源：到最近来源支持的距离≥3/4−1/16−1/16=5/8>1/2。Q_2 完整捕获恰相邻两 packets，最远端距离≤7/8<1；其它来源支持距离>1。完整 averages 为0、1，故原 R=2 唯一、M=1>τ、m=2；同样 anchor-cold 对任意 K 成立。

每个 interval 长1/8、annular 指示捕获这两个来源全质量。来源两端 packets 只有一份 interval，内部 packets 有两份。因此真实原 profile 为

\[
 V(y)=\begin{cases}1/64,&\text{两端 packet},\\
                    1/32,&\text{内部 packets}.
       \end{cases} \tag{11}
\]

取 H=1/128，所有来源都满足 V>H，故真实 B 仍覆盖全 μ。相邻扩厚 intervals 的宽为2+d=17/8，大于 center 间距3/2，整个 U 只有一个组件，质量 w_1=J、直径随 J 增长。

在每个原 row 上 c=2、b_B=1。J>8 时 c=2<ηw_1，故 **全部 rows 都是 local fragment**；local-fullness 捕获份额0。这是完整真实输入、原 winner 与真正高势 B 的反例，否定一个无 H 限制的“高势必使该组件覆盖有常数 fullness”论断。它没有按 x 重组来源。

准确 first moment 与尾是

\[
 I_R=(J-1)/32,
 \quad\int Vd\mu=2/64+(J-2)/32=I_R,
 \quad F(H)=(3J-4)/128. \tag{12}
\]

两族的 E_R 都是原合法子集，不是完整 {M>τ}。式 (12) 不能当成完整 E 的一般弱比。H 是小固定常数；随着 n 的目标 cutoff≈√n 增加，这个低 profile 会退出 B。本例不反驳目标阶数，也未实现原 soft FIRST/fullfuture/CPGP/LCA/history。它只说明“B 是高集”的逻辑条件与空间组件唯一标签，本身不足以排除长组件 fragment；还必须实质使用高 cutoff 和完整实际几何。

把连通组件再细分为每个小 packet 会恢复这个例子的 c_j=w_j，但每 row 同时见两个标签。一般再细分后标签数可任意大，不再有 (3)；须核验真实 row multiplicity、分配权重或其他跨片预算。不能直接把 (5) 逐小 packet 求和并仍报同一 I_R。这里恰有2标签，不证明一般重叠族也有常数2。

## 6. 剩余接口与 actual gates

当前一般局部化没有 receiver 自适应 packets、重选赢家或重复扣 source 质量。它严格减少 global fragment，并将所有未付 row 归到一份较大、可能很长的真实 B 组件。原 B 的高势条件逐组件保留为 (4)；可能的下一步应控制这些长组件内的重叠，而不是再次要求一层固定 packets 普适覆盖。

若要用空间树逐层切组件，必须明确原 annular incidences 在每层如何分配、每个 row 的总权重及每个 source 的 H 负项消耗；单层 Σw_child=w_parent 不自动给多层总质量≤w。E03/G01 的树望远镜或 packing 不能替这些字段。重叠链已有真实持续残余，但不是大 H 反证；仍需证明 (8) 或在原输入中构造大 H 的真实对抗例。

对带原 g(x,y)∈[0,1] 的 actual V^g，可以从该势的 B^g 定义 ξ、组件，再以 c_j^g=∫A_x g dξ_j 使用同一单标签/fullness 预算，因为正支持仍在原 Q_R。这保留所有已经积分进 g 的 FIRST、early/fullfuture、CPGP、LCA、all-history。实际 weighted local-fullness 可能更少，不能用未加门的 c_j 替它再支付原 g 交通；纯 hard 的 first-moment 保留率亦不能无证明迁移。本文不改变原 geom 主账。

## 7. 完整 S 的 near-flat 约束：主要尾不能靠 rare-source 消失

本节只接 [log-max 变分 §5](log_max_variation_20261007.md) 的**完整**势；不赋予受限 V_hi 同一性质。令 E={M>τ}、A_x=Q(x,R(x))、μ_hi=μ，定义完整 S=τ∫E1_Q/m。前面所有 component/fullness 证明只用 A_x⊂Q，故这套选择也合法。假定当前完整 L¹ 输入满足一份已校准的 near-flat 合同

\[
 d:=\frac1W\int|S-K|d\mu,
 \qquad 0<H<K.
\]

这里 K 是同一 input 类的 log-max 上确界常数，或只是给定 near-flat 比较值；不预设它≈√n。root 的近优化归约在其相同合法输入类上给 d≤3√(εZ)，Z=1+nlog(b/a)。若极值类包含原子，不能不经原正宽化/近似接口把其 near-flat 输入自动当成本文 L¹ 组件实例。

由 |S−K|≥K−H 在 {S≤H} 上，以及 positive-part 的1-Lipschitz性，严格有

\[
 \frac{\mu\{S>H\}}W\ge1-\frac d{K-H},\qquad
 \left|\frac{F_S(H)}W-(K-H)\right|\le d,
 \qquad\left|\frac IW-K\right|\le d. \tag{13}
\]

负的质量下界只作0下界；无需新数值。令 q=δ₀+η⁻¹e⁻ᴴ⸍²，(6) 在完整 S 上反给一个必要的 fragment 交通下界：

\[
 \frac{J_{G,\rm loc}}W
 \ge(1-q)K-\frac H2-(1+q)d
                -\frac{Hd}{2(K-H)}. \tag{14}
\]

推导仅是 J_G≥F_S−qI+(H/2)w，分别代入 (13)；没有新 paid 费。当 K≫H、d≪K 且 q<1 固定时，B 的来源质量几乎全部保留，fragment 必须承载主要尾。故不能期待 nearmax 的来源高集自动稀少，亦不能据 nearflat 获得免费 mass shrink。

如果将来对这批完整输入真正证明尾 F_S≤δI（δ<1），由 (13) 才能得 (1−δ)K≤H+(1+δ)d，并在同一类合法近优化极限中限制 K。当前 (14) 只定位必须支付的真实交通，没有证明这份尾合同。对于 anchor/high-source/FIRST/LCA 受限 V_hi，profile 的门及其 input 变分项不同，不能照搬 (13) 或宣布其 B 质量同样近全。

## 8. 预登记空间证书与失败保留

原 [registration](source_fragment_localization_20261007_registration.md) 后执行 [guard](source_fragment_localization_20261007_guard.py)；在远离族开 interval 的下端误要求严格距离，首轮失败。精确端点是9/16−1/16=1/2，来源面零质量且内部严格。保留 [failure](source_fragment_localization_20261007_failure.json)，没有改原输入或阈值。

再单独 [v2 登记](source_fragment_localization_20261007_v2_registration.md)，仅将这一端点检查改为非严格，执行 [guard v2](source_fragment_localization_20261007_guard_v2.py) 一次。[v2 results](source_fragment_localization_20261007_v2_results.json)、[v2 receipt](source_fragment_localization_20261007_v2_receipt.json) 保存脚本/登记 SHA256。J=12,24,48，75 个精确 Fraction 断言 PASS；远离族组件数12/24/48，重叠链组件数1/1/1。

这些是两个新的一维原有限赢家空间族的解析构件证书：核完整捕获、唯一 winner、annular 门、profile 与真实 B、组件/fullness/fragment、一次质量抵扣。没有 MC、背景质量、原 soft gates 样本或一般阶数估计。原 failed run 未改成 PASS，旧守卫未重跑。

本文新增的严格一般结果是替代预算 (6) 及 local fragment 缩减；远离重复簇不再是该接口的障碍。长重叠组件的 (8) 仍未付，且小 H 真例子表明不能靠组件存在或高集名称自动消除它。一般 sqrt(n)n^{o(1)} 与原 geom 主核保持未解决。
