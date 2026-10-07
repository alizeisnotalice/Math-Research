# 原高源尾：真实高势集合的一次质量抵扣与未付 fragment

2026-10-07。父任务指定 6.1-sol high。本文证明一个严格的尾子支：从原 V_hi 自身定义一份全局来源集合 B，在原查询对 B 有固定 fullness 的 receiver 上，来源重排可用截断尾中的负项抵扣一次 W_B，余费系数为 η⁻¹e⁻ᴴ⸍²。剩下的小比例捕获 fragment 仍未证明 good-lambda；没有闭合一般 sqrt(n)n^{o(1)} 或原 geom。

使用已读 [E04 SKILL](</Users/zhengzhihao/.codex/skills/math-e04-snell-envelope-optimal-stopping/SKILL.md>) 及其 method/cube-interface 检查停止迁移边界。下文是确定性 source–receiver 交换，不调用 Snell 或可选停止。既有 A01/E01 的完整来源归一化保留。

查重：[single-capture §2–4](single_capture_source_rearrangement_20261007.md)、[完整 packet square](general_source_packet_square_20261007.md)、[direct mass-band §3](direct_capture_mass_band_square_20261007.md)、[aggregate 剩余尾合同](aggregate_mass_band_budget_20261007.md)、[source-tail](source_occupation_tail_geometry_20261007.md)。径向 log 重排、min-cap 核和截断尾充分方向均有旧构件，不能报成新机制。本文新连接是：**B 由真实整份原高势定义，receiver 积分前只冻结这一份 B；其自依赖合法，负 H W_B 只用一次。** 不再假定通用 mass-band pair Gaussian 衰减或固定 packets 常数覆盖。

## 1. 冻结完整输入、原赢家与原余项

设 μ=fdy≥0、W=∫f<∞，μ_hi≤μ 是上一稿已经冻结的高源限制。原尺度窗口为 a≤R(x)≤b≤2a，a>0；原真实 winner 为 R，完整平均和捕获质量为

\[
 M(x)=h_{R(x)}*\mu(x),\qquad m(x)=M(x)R(x)^n.
\]

E_R 是上一稿保留的原 receiver 子集；只要求 M>τ>0，并令 I_R=τ|E_R|<∞。原 anchor-cold、低源删除等门不重选。中心 cube Q(x,L) 的边长为 L，h_L=L^{-n}1_Q；边界按原约定，L¹ 来源对固定查询边界零质量。

单 anchor 的正剩余核准确是

\[
 (h_R-h_a)_+(x-y)=R(x)^{-n}
          \mathbf1_{Q(x,R(x))\setminus Q(x,a)}(y).
\]

记 A_x=Q(x,R(x))\setminus Q(x,a)。本次研究的实际纯 hard 高源余势与尾为

\[
 V(y)=\tau\int_{E_R}\frac{\mathbf1_{A_x}(y)}{m(x)}dx,
 \qquad F(H)=\int(V-H)_+d\mu_{\rm hi}. \tag{1}
\]

它就是上一稿的 V_hi，不增加背景、不更换分母。联合核可测时 Tonelli 给

\[
 J_{\rm hi}:=\int Vd\mu_{\rm hi}
 =\tau\int_{E_R}\frac{\mu_{\rm hi}(A_x)}{m(x)}dx\le I_R. \tag{2}
\]

又因 τ/M<1，V(y)≤∫sup_{a≤L≤b}h_L(x-y)dx=1+nlog(b/a)。故所有定义有限；finite winner 的可测性直接成立，continuous winner 沿原合法可测 selector。H≥这个 O(n) cap 的尾当然为零，不能把它称作平方根进展。

## 2. 真正高势来源集合与精确一次抵扣

固定 H>0，定义

\[
 B=\{y:V(y)>H\},\quad \xi=\mu_{\rm hi}|_B,
 \quad w=\xi(\mathbb R^n).
\]

B 可以依赖整份 input、原 R、E_R 与 H。它在 source 积分前是一份全局可测集合，不能因当前 x 改变，也不在重选 μ 的 winner 后重新判定。若 w=0，F(H)=0；以下设 w>0。定义

\[
 c(x)=\xi(A_x),\qquad b_B(x)=c(x)/m(x)\in[0,1].
\]

准确的 source–receiver 交换是

\[
 \boxed{F(H)=\tau\int_{E_R} b_B(x)dx-Hw.} \tag{3}
\]

其来源是 ∫_B Vdμ_hi−Hμ_hi(B)，没有改成任意来源 supremum，也没有条件 weak 平均。由 (2) 只有定义性的 w≤I_R/H；不存在由此得到 w/W≪1 的一般结论。高势有可能覆盖几乎整份来源。即使 w≤I_R/H，I_R/W 本身正是尚未付的弱比，不能循环以它支付自身。

取 δ₀∈(0,1)、η∈(0,1]，在原 receiver 上切三个不交集合

\[
 \begin{aligned}
 D_0&=\{b_B\le\delta_0\},\\
 D_F&=\{b_B>\delta_0,\ c\ge\eta w\},\\
 D_G&=\{b_B>\delta_0,\ c<\eta w\}.
 \end{aligned} \tag{4}
\]

所有集合都与 E_R 相交。令 I_F=τ|D_F|、J_j=τ∫D_j b_Bdx。直接有 J_0≤δ₀ I_R。D_F 对这一份真实 B 有 fullness；D_G 则只捕获 B 的小比例，却在完整原查询中占 >δ₀ 的质量份额。这个事实不能叫几何 diffuse，也不能通过细分 source 标签制造。

## 3. J_full 的完整重排与已付尾子支

对每个 y∈B，令

\[
 P_F(y)=\tau\int_{D_F}\frac{\mathbf1_{A_x}(y)}{m(x)}dx.
\]

在实际 incidence y∈A_x⊂Q(x,R(x)) 上，D(x,y)=(2||x−y||∞)^n≤R(x)^n。原阈值与 fullness 分别给 m(x)>τR(x)^n 和 m(x)≥c(x)≥ηw。因此

\[
 \frac\tau{m(x)}\le
    \min\left\{\frac\tau{\eta w},\frac1{D(x,y)}\right\}. \tag{5}
\]

D=0 时把第二项解释为 +∞，min 仍有限。固定 y 后，|{x:D(x,y)≤v}|=v。设 v_F=|D_F|，φ(v)=min(τ/(ηw),1/v)。φ 非增，层饼给任何 |A|=v_F 的集合

\[
 \int_A\phi(D(x,y))dx
 =\int_0^{\tau/(\eta w)}|A\cap\{\phi(D)>t\}|dt
 \le\int_0^{\tau/(\eta w)}
             \min(v_F,|\{\phi(D)>t\}|)dt
 =\int_0^{v_F}\phi(v)dv.
\]

把原 incidence 指示上放至 D_F 只增积分，得到

\[
 P_F(y)\le\ell\left(\frac{I_F}{\eta w}\right),\qquad
 \ell(r)=\begin{cases}r,&0\le r\le1,\\1+\log r,&r\ge1.\end{cases}
\]

这既不要求 B 在空间小、来源坐标独立，也不要求 D_F 为 cube。积分一份 ξ 得

\[
 \boxed{J_F=\int_B P_Fd\mu_{\rm hi}
       \le w\ell\left(\frac{I_F}{\eta w}\right).} \tag{6}
\]

ℓ 是全局凹函数；在 T=e^{H/2}≥1 的切线为

\[
 \ell(r)\le\ell(T)+\ell'(T)(r-T)
          =\frac H2+e^{-H/2}r.
\]

所以严格的已付部分为

\[
 \boxed{J_F\le\frac H2w+\eta^{-1}e^{-H/2}I_F.} \tag{7}
\]

代入唯一的 (3)，

\[
 \boxed{F(H)\le
    [\delta_0+\eta^{-1}e^{-H/2}]I_R
                 +J_G-\frac H2w.} \tag{8}
\]

式 (8) 已证；J_G 仍未付。此处 H w 的一半抵扣 fullness，只余一半供 fragment，不能把两部分各抵扣完整 H w。若另分别控制 J_F−Hw 和 J_G−Hw，不能把两界直接相加充作 F(H) 的上界：这会多减一次 Hw，而 (3) 只有一个负项。

旧 packet 稿用 τ/m≤2τ/(ηM_i+τD) 得 2log(1+I_i/(ηM_i))；direct mass-band 稿已经使用 (5) 类 min-cap。这些是同一径向来源几何的相关构件。本文 (6) 的 sharp ℓ 不改变一般阶数，真正新增的是针对真实高 V 集合的 (3)、(7)、(8)，而非重新报一份 packet 平方预算。

## 4. 未付 fragment 的原几何约束与最小合同

在 D_G 上有

\[
 m(x)=c(x)/b_B(x)<\eta w/\delta_0,
 \qquad R(x)^n<\eta w/(\delta_0\tau). \tag{9}
\]

故它同时保留原小完整捕获质量门与半径上限
R≤b、R<(ηw/(δ₀τ))^{1/n}；允许原 R=b，不能把第一个不等式也写成严格。若 ηw/δ₀≤τa^n，则与原 m>τR^n≥τa^n 矛盾，D_G 为空。这是可验证的特殊 paid branch；一般 w 大时不够。

由 (8)，一个比原全尾更具体的充分缺项是

\[
 \boxed{J_G-\frac H2w\le\delta_1 I_R.} \tag{10}
\]

若 (10) 证明，则全尾 ≤(δ₀+δ₁+η⁻¹e⁻ᴴ⸍²)I_R。这里 δ₁ 是待证 allowance，不能登记为已付费。也可定义 V_G(y)=τ∫D_G1_Ax(y)/m(x)dx，则

\[
 J_G-\frac H2w
 =\int_B(V_G-H/2)d\mu_{\rm hi}
 \le\int(V_G-H/2)_+d\mu_{\rm hi}. \tag{11}
\]

式 (11) 只是正部上放；为 V_G 重新取高集再无限递归，会把 cutoff 逐次减半，并不能仅凭定义得到收缩。B 已冻结，来源质量抵扣有限；任何新分层必须证明真正的几何收缩或跨层不重付合同。

原范围下，D_G 可以把 B 分成许多空间互相远离的片，每个原 query 只见一个小片；(9) 不限制片数，也不控制同一 source 随不同原查询反复进入 fragment 的次数。这不是已认证的反例，只是 (9) 所尚未包含的量。固定 packet 覆盖已失败的结果不能反向当成这个实际高势 B 的尾反例；需从 B={V>H} 的自一致性证明或构造真实输入。

## 5. restricted-source maximal：合法包含与循环边界

令 M_ξ(x)=sup_{L∈[a,b]}h_L*ξ(x)，为这份 B 的 restricted maximal。D_G 上原 annular ξ 响应为 c/R^n>δ₀M，因而

\[
 D_G\subset\{M_\xi>\delta_0M>\delta_0\tau\}. \tag{12}
\]

若原 E_R 保留 M>K A_aμ，则 ξ≤μ 又给

\[
 D_G\subset\{M_\xi>\delta_0K A_a\xi,
                       \ M_\xi>\delta_0\tau\}. \tag{13}
\]

这是一个真正保留的 cold 事件包含，使用完整原分母，不按 x 重组 source。它没有说 ξ 自己的 winner 等于 R，更没有把原 annular profile 变成 ξ 新赢家的 profile。

任意 t∈(0,1) 的完整 b_B>t 集合均在 {M_ξ>tτ}。已知 bounded-window positive envelope 质量 C_n=1+nlog(b/a) 只给

\[
 \tau|\{b_B>t\}|\le\min(I_R,C_nw/t),
 \qquad\tau\int b_Bdx
 \le\begin{cases}
 I_R,&I_R\le C_nw,\\
 C_nw[1+\log(I_R/(C_nw))],&I_R>C_nw.
 \end{cases} \tag{14}
\]

第二式由 layer-cake 积分第一式，b_B≤1。这仍带 O(n) 常数。若在 (13) 调用未知一般 cold weak 常数 C_cold(δ₀K)，会得到 |D_G|≤C_cold(δ₀K)w/(δ₀τ)：只是把待解决的 cold weak 问题移给 ξ。还没有契约保证它比原常数小，δ₀K 的冷门甚至变弱，不能据此声称 good-lambda。

也不能用 w≤I_R/H 后把 C_nw/δ₀ 当成小费：这只得 C_n I_R/(δ₀H)，H≈√n 时仍大。原 B 的来源质量占比可以接近1；没有免费 W_B/W 收缩。

## 6. winner 最优性实际提供什么，尚未提供什么

拆完整 μ=ξ+ζ（ζ≥0），记 A_Lν=h_L*ν、M_ν=sup_{L∈[a,b]}A_Lν。原 true winner 保证 M=M_μ，并对每个允许 L 有

\[
 A_L\xi-A_R\xi\le A_R\zeta-A_L\zeta. \tag{15}
\]

它准确说明 ξ 在新尺度的收益可以被完整互补来源的平均损失抵销。不能从完整 μ 的 winner 推出 ξ 的 winner，更不能忽略右边。另有

\[
 0\le M_\mu-M_\zeta\le A_R\xi,
\]

其中右界因 M_ζ≥A_Rζ，左界因ζ≤μ。上述式子都是点态强比较；它们没有为 ∫_D_G c/m 的同源重复次数提供已付空间积分。原 cold 门涉及 A_aμ，而 fullness 与 fragment 涉及实际 B，两个方向没有现成的乘法收缩。

E04 的 Snell 构件要求共同滤过、适应可积奖励与条件期望递推。当前固定 x 的尺度平均只是完整输入的确定性函数；不同 x 的“首次”不是一个已供给的同源鞅停时，winner 最大值也不是可选抽样结论。单个中心的嵌套 cube 链不直接控制全部中心的 source occupation。因此本文不从 (15) 或 deterministic first threshold 推出可选停止费。

原 FIRST/fullfuture、CPGP、唯一 LCA 与 all-history 在本文均冻结。若对带原 g(x,y)∈[0,1] 的 actual traffic 使用相同方法，必须从它自身定义 V^g 和 B^g，并令 c^g(x)=∫B^g1_Ax(y)g(x,y)dμ_hi(y)。则 c^g≤m、低比例/fullness 重排 (3)–(8) 仍逐项成立；门 c^g≥ηw 是实际 weighted fullness。但它可能覆盖更少，不能删 g 后重领交通，亦不能由纯 hard 的 J_hi≥(1−ν−K⁻¹)I_R 自动得到 actual 子交通的同一下界。原主账的吸收常数要独立校准。本稿没有将这份纯 hard allowance 记入 geom paid 总账。

## 7. 三轮标量证书及收口

先登记后执行一次：[registration](general_tail_goodlambda_20261007_registration.md)、[guard](general_tail_goodlambda_20261007_guard.py)、[results](general_tail_goodlambda_20261007_results.json)、[receipt](general_tail_goodlambda_20261007_receipt.json)。n=16,64,256；H=2(√n+2)log2；η=1/4、δ₀=1/8。59 个 Fraction 断言 PASS，full fee η⁻¹e⁻ᴴ⸍² 分别为 1/16、1/256、1/65536。log2 用有理 atanh 区间，核切线多个参数点、原 O(n) cap 下的 H、标量分类和 fragment mass cap。

额外核验的 δ₁=1/4、ν=1/4、K=√n 只表示：**若** fragment 合同 (10) 另被证明，这些常数满足上一稿的尾吸收条件。没有验证 (10)，没有原空间积分、gate 样本、MC 或新阶数实验，旧数值没有重跑。

当前严格所得为 (8)：真实高势 source 的 low-share 支可吸收，fullness 支由一次 W_B 抵扣加指数小费支付。当前最小未证项是 (10)，连同原完整 μ 的 winner/cold/fullfuture 门，在小捕获比例却大原份额的真实 fragment 上控制同源总 occupation。它比任意 pair 衰减或强 square 要求更弱、字段明确，但仍未成为已证 good-lambda；一般 geom 主核保持未付。
