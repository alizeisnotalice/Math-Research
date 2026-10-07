# Hard truewinner 的包络熵：宽径向混合、角坐标及面 trace 缺口

2026-10-07。只新增本稿，不改其它文件或主账。按父任务 6.1-sol high 要求工作；工具没有公开实际运行模型 ID。本轮没有新已证阶数，没有数值或 toy 守卫。

**结论：有限 window 给出的 joint KL≤log(Z/m) 正确，但这是既有 envelope entropy 身份。** 新审计明确：参考噪声是 cube-cone 方向与一个宽 O(n) 的径向混合，只有条件角坐标具有产品结构；原 hard winner 的尺度导数则含相对 bulk prior 奇异的 surface trace。尚未证明一个原几何 score 的 P-均值支配 m，故没有 Gibbs 平方根付款。

读取 [D02](/Users/zhengzhihao/.codex/skills/math-d02-gibbs-entropy-duality/SKILL.md) 及 method/provenance/cube-interface；只使用自含的 KL 非负/Gibbs 概率代数，不引用未读的一般 Donsker–Varadhan 定理。全文检索现有 entropy 稿，已读旧 [high_column_direction_entropy](../high_marked_cover/high_column_direction_entropy.md)、[source_operator_entropy_localization_audit](../high_marked_cover/source_operator_entropy_localization_audit.md)、[entropy_characteristic_route](../signed_band_probability/entropy_characteristic_route.md)、[selected_source_chain](../../cube_general_20261002/fresh_general_cover/conditional_cube_transport/channel_entropy/selected_source_chain.md)，及新 [radial_shared_uniform_budget](radial_shared_uniform_budget_20261007.md)。旧 HE2/HE3/HE6 已包含一列/global entropy 和角集中；不重复记为新成果，不重做旧 entropy 假设反例。

## 1. 原 truewinner、weak joint 与参考概率

采用 side-length 约定

\[
 h_L(z)=L^{-n}\mathbf1_{\|z\|_\infty\le L/2},\qquad
 a\le L\le b\le2a,\quad a>0.
 \tag{1}
\]

原完整正来源 μ 的质量为 0<W<∞，可具有任意坐标依赖。令 H(x)=max_L(h_L*μ)(x)、E={H>τ}、τ>0，选择原 measurable truewinner L(x)，g(x)=τ/H(x)·1_E。finite family 使用原 tie rule。连续 family 以原 L1 来源 μ=fdy 为例，(x,L)↦h_L*f(x) 连续：平移/尺度改变的 cube 对称差体积趋零，f 的绝对连续积分控制差。compact window 的最大值存在；最小 maximizer 可测，因为 {L(x)≤c}={H_[a,c](x)=H_[a,b](x)}。这里不引入 restricted-source winner。

定义原 selected joint 及其质量

\[
 J(dy,dx)=μ(dy)g(x)h_{L(x)}(x-y)dx,\qquad
 m=\frac{τ|E|}{W},\quad J(\mathbb R^{2n})=mW.
 \tag{2}
\]

E 非空时 m>0。设 r(z)=2||z||∞、ℓ_0(z)=max(a,r(z))，正包络为

\[
 k_*(z)=\ell_0(z)^{-n}\mathbf1_{r(z)\le b},\qquad
 Z=\int k_*=1+n\log(b/a).
 \tag{3}
\]

令 Q 在 (y,z) 坐标为 (μ/W)(dy) k_*(z)dz/Z，并令 x=y+z。**Y 与 noise Z_noise 独立，Y 的坐标仍为原 μ 的任意联合 law。** Q 不是 μ 的产品坐标参考，也不是原 receiver posterior。

H≤k_* *μ，正 Fubini 给 τ|E|≤ZW。这只是已知 bounded-window strong envelope。P=J/(mW) 为概率，Q 是一份完整来源概率，未按尺度或方向重新领 W。

## 2. 真正 joint entropy 与 source Palm

在原 incidence 上 L(x)≥ℓ_0(x−y)，定义 d(y,x)=n log[L(x)/ℓ_0(x−y)]≥0。零 incidence 的 RN 取 0。准确

\[
 \frac{dP}{dQ}(y,x)=\frac Zm g(x)e^{-d(y,x)}
       \mathbf1_{x\in E,\ r(x-y)\le L(x)},
 \tag{4}
\]

所以 RN≤Z/m。因为 m≤Z、0≤d≤n log(b/a)，且 H≤Wa^{-n} 给 g≥τa^n/W>0 在 E 上，涉及的 logarithms均可积。直接取 log 得既有完整 joint 身份

\[
 D(P\|Q)+\mathbb E_P[d-\log g]=\log(Z/m),\qquad
 0\le D(P\|Q)\le\log(Z/m).
 \tag{5}
\]

令真正 source column

\[
 S(y)=\int_E g(x)h_{L(x)}(x-y)dx,\quad0\le S\le Z,
 \quad\int S\,dμ=mW.
 \tag{6}
\]

Palm source 为

\[
 \frac{dP_Y}{d(μ/W)}=\frac{S(y)}m.
 \tag{7}
\]

它不是 μ/W；source coordinates不会因 Q 的噪声角坐标而成为 independent。S=0 的来源没有 P mass。条件 receiver P(dx|y)=g h_{L(x)}(x−y)dx/S(y)，其相对 prior k_*/Z 的 KL 为

\[
 D_y=\log[Z/S(y)]+\mathbb E_{P(\cdot|y)}[\log g-d],
 \quad
 D(P\|Q)=D(P_Y\|μ/W)+\int D_y P_Y(dy).
 \tag{8}
\]

全部 KL 非负，完整 entropy 预算同时支付 source Palm change 和 conditional receiver change，不能分别各领 log(Z/m)。P_X=Lebesgue|E|/|E| 是真实 weak calibration；P(dy|x)=μ|Q(x,L(x))/μ(Q(x,L(x)))。这个原来源 posterior 一般不具产品结构。

## 3. 包络 prior 的准确 cube-cone 分解

除 noise 为 0 或 max-face tie 的 Lebesgue 零测集合，写

\[
 z=\frac R2\Omega,\quad R=2\|z\|_\infty>0,
 \quad\|\Omega\|_\infty=1.
 \tag{9}
\]

cone 概率 ς_n 如下：face index I uniform于 {1,…,n}，face sign ε uniform于 ±1，Ω_I=ε，其余 n−1 个 Ω_i iid Uniform([-1,1])。体积公式为 dz=nR^{n−1}dR dς_n(Ω)。由 (3)，Q 下 Y、R、Ω 彼此独立，R 的 law为

\[
 q_R(dR)=\frac nZ\left[
 \frac{R^{n-1}}{a^n}\mathbf1_{0<R<a}
 +\frac1R\mathbf1_{a<R<b}\right]dR.
 \tag{10}
\]

core 质量 1/Z，core conditional noise是 Uniform(Q_a)，其未条件化的 n个 Cartesian noise坐标确实 independent。outer 质量 B/Z，B=n log(b/a)，outer R logarithmically uniform；b=a 时 outer branch 不存在。

特别，T=n log(R/a) 的完整径向 density 为

\[
 q_T(dt)=\frac1Z\{e^t\mathbf1_{t<0}+\mathbf1_{0<t<B}\}dt,
 \quad Z=1+B.
 \tag{11}
\]

core T为负 Exp(1)，outer T为 Uniform([0,B])。其 MGF准确为（β>−1，β=0取连续值）

\[
 \mathbb E_Qe^{βT}
 =\frac1Z\left[\frac1{1+β}+
                 \frac{e^{βB}-1}{β}\right].
 \tag{12}
\]

若 window ratio固定大于1，B为 O(n)，outer conditional径向 variance B²/12为 O(n²)。因此包络 noise不是一份整体的 Cartesian product law；宽径向混合不能按 n个 bounded independent increments 当成 √n concentration。这里 R 是相对 source 的 norm radius，通常不是 selected L(x)；二者差正是 (5) 的 depth d。

真正保留的 independent angular坐标只在 Q 的条件 (Y,R,I,ε) 下成立。把 deterministic selected label L(Y+Z_noise) 附加到 Q，不改变 KL，却会条件化这些角变量；它不是额外 independent label。

链式 KL还可精确拆为 source、radial、cone 三项：

\[
 D(P\|Q)=D(P_Y\|μ/W)
 +\mathbb E_{P_Y}D(P_{T|Y}\|q_T)
 +\mathbb E_{P_{Y,T}}D(P_{\Omega|Y,T}\|ς_n).
 \tag{13}
\]

这保留 arbitrary μ，不估计其 differential entropy；每项都在同一总 log(Z/m) 内。P 下 radial/source/direction可相关，不能把 Q 条件独立当作 P 独立。

## 4. Gibbs 的角工具已有，但 m 的支配桥未证

旧 HE3 的合法 centered angular score为

\[
 V(\Omega)=\sum_{i\ne I}(|\Omega_i|-1/2),\qquad
 \mathbb E_Qe^{tV}=
 \left[\frac{\sinh(t/2)}{t/2}\right]^{n-1}
 \le e^{(n-1)t^2/8}.
 \tag{14}
\]

它与 source/radius 无关，故 arbitrary μ 不影响此 MGF。D02的 Gibbs 证明是令 dQ_t=e^{tV}dQ/E_Qe^{tV}，由 D(P||Q_t)≥0 得 tE_PV≤D(P||Q)+log E_Qe^{tV}。正负 t分别应用并优化，给既有

\[
 |\mathbb E_PV|\le\sqrt{\frac{n-1}{2}\log(Z/m)}.
 \tag{15}
\]

允许 coefficients依赖先给定的 (y,R,I,ε)，若其 bounded-square norm一致，conditional MGF同样可用；若 coefficients/测试方向在看完整 Ω 后选出，原 product-MGF证明未涵盖它。尤其不能令来自 original source posterior 的向量免费独立于该方向。

若能额外证明一条**原 truewinner 几何**不等式 E_PV≥c(m−1)−C，其中 c>0、C≥0 为维数无关常数（或另一真正 centered score具有相同 MGF并支配 m），才会推出

\[
 m\le1+C/c+c^{-1}\sqrt{\tfrac{n-1}{2}\log(Z/m)}.
 \tag{16}
\]

这将给想要的 Gibbs 路线阶，但这份 mean-domination桥尚未得到。式(16)仅是 conditional代数，不是原上界的新估计，不据此运行数值。entropy身份没有迫使 (14) 的 centered mean为正，更没有迫使它增长如 m。

一个看似能支配 m 的循环候选是 source score S(y)：E_PS=∫S²dμ/(mW)≥m，而 E_QS=m，Q-MGF为 ∫e^{tS}dμ/W。它依赖未知 source column，高均值本身就是未付量；仅用 S≤Z给线性 envelope，不能登记为 √n MGF。source likelihood log(S/m) 的 P均值恰是 source KL，其 t=1 Q normalizer为1，只重写 (8)。joint log-RN score亦只重写 (5)。这些都不是新几何势。

## 5. Hard dilation 的 singular surface 项

原 hard kernel的 distributional尺度导数准确为

\[
 \partial_{\log L}h_L
 =-nh_L+\mathcal B_L,
 \quad
 \mathcal B_L(dz)=\frac{L^{1-n}}2
 \sum_{i=1}^n[\delta_{L/2}(dz_i)+\delta_{-L/2}(dz_i)]
 \prod_{k\ne i}\mathbf1_{|z_k|\le L/2}dz_k.
 \tag{17}
\]

每对 faces质量1，故 ∫B_L=n，∫∂logL h_L=0。没有遗漏 factor L/2 或两面。B_L相对原 bulk k_*(z)dz奇异：它支撑于固定 L 的 cube faces，后者 Lebesgue 零测。

对 smooth来源，在有合法 trace的 interior stationary hard maximum，

\[
 L\partial_L(h_L*μ)(x)=0
 \quad\Longrightarrow\quad (\mathcal B_L*μ)(x)=n(h_L*μ)(x).
 \tag{18}
\]

这是 surface posterior和bulk response之间的接触关系，不是 bulk accepted law下一个有限 centered angular score的均值身份。general L1来源的 surface convolution仅有 a.e.尺度/coarea或 distributional意义；original selected L(x)可以落在 exceptional trace，不能未经 regularization/trace迁移就在该点代入普通导数。endpoint winner也只给 one-sided条件。

若在带 L标签的参考上处理 (17)，其 face measure相对 dz dL仍支撑在 codimension-one集合 |z_i|=L/2；需要另一个 surface prior/coarea预算。原 selected label L(x)的 graph不能把这份新measure自动变成 (4) 的 bulk prior。arbitrary μ的产品化也不能补这个缺口。

**固定 L 或 dz dL 下的奇异性，不推出自适应 selected graph 的 face事件仍是零概率。** 例如 μ=δ_0，连续 window 内 a<2||x||∞<b 的真实 hard winner为 L(x)=2||x||∞：捕获原点后的响应 L^{-n}随 L下降，因此最小捕获尺度恰为 winner，selected graph正落在 face。此例只提醒自适应 graph的量词区别，不作为原 L1/全部 actual门的反例。这里的正确结论仅是 fixed-L bulk KL不能自动控制 selected trace/coarea；迁移必须按真实自适应图像另行证明。

## 6. 低 bulk KL 为何不能直接支付 trace

固定一个 L_0∈[a,b]，在一面内部选 subface |z_k|≤L_0/4（k≠1）。取厚度 ε 的 inward strip z_1∈[L_0/2−ε,L_0/2]，ε足够小时位于 Q_b。Q_noise密度 k_*/Z≥b^{-n}/Z，该 strip概率至少

\[
 c_{n,a,b,L_0}\,ε,
 \quad c_{n,a,b,L_0}=\frac{(L_0/2)^{n-1}}{Zb^n}>0.
 \tag{19}
\]

用可测 score V_ε=ε^{-1}1_strip模拟单位 surface trace，则任何固定 t>0 都有

\[
 \mathbb E_Qe^{tV_ε}\ge c_{n,a,b,L_0}ε e^{t/ε}
 \longrightarrow\infty.
 \tag{20}
\]

这是固定 n下精确的 positive-exponential-moment障碍，不是 fitted variance/δ tradeoff，也不产生任何 2/3阶推断。把尺度 derivative再除bulk density构造 likelihood-score，会保留这种随 trace thickness增长的高度；单靠 (5) 没有一个 uniform finite MGF可用。

(19)–(20)只审计 Gibbs直接作用于 raw face trace的接口，未给完整 actual weak反例，也未声称原 surface项不可由别的几何支付。可以更换 prior去含 surface、平滑尺度并清算其新normalizer，或证明原winner限制下的source/face占用；每条都需要另一个真实费用，不能继承 Q的logn KL并删除奇异性。

## 7. 本轮可操作缺口及完成状态

joint/source Palm/KL身份是旧 envelope构件。新审计给出 (11)–(12) 的宽径向prior与 (17)–(20) 的精确 hard-face缺口。真正下一步需要以下之一：

* 一个原 winner/fullfuture/实际门确实强迫的 bulk geometric score，其 P-均值支配 m，且在保持 arbitrary μ、radial mixture与post-selection dependence后有可核 MGF；
* 一个独立的 surface/coarea或 source-tail预算，将 (18) 的 trace转成可积source费用，再明确清算相对参考的改变。

仅证明 small KL、重新命名 source S或将prior换成 products不会完成这一步。pure-hard truewinner工具还没有原 soft FIRST/early/CPGP/LCA全部资格；也不能直接把本稿的 J认作当前 R_angle。没有新的主账 paid费，一般 √n n^{o(1)}目标仍未证明。现应回到原 actual source-tail/列占用的真实probe或解析几何，而非继续对已有 angular MGF堆toy检查。
