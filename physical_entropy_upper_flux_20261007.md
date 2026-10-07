# 物理尺度熵的上通量：标量变差、实际截门与局部正部缺口

2026-10-07。本轮只新增本稿，不改主账，不重跑旧数值。结论是**没有取得新的任意正 L¹ 输入的 √n polylog n 上通量合同**。更前面的实际回代缺口可以精确定位：标量熵正变差保留空间取消；不同接收点的赢家尺度和原接受门会选择空间局部通量。一个完整正 L¹ 窄盒输入证明，删门后要求接收点局部正变差的 source-once 预算一般至少为线性 nW。因此标量 √n 变差上界即便成立，也不能通过这个删门步骤支付 actual 余项。

这不是 actual FIRST/geom 反例：下面窄盒只检验一个一般放大接口，没有认证原 FIRST、CP/GP、远对、lowS、lowcoin、诊断 K、角门及自适应 core 等完整资格。本稿不给任何新已付项。

## 1. 原输入、已知合同和真实面导数

固定完整同一输入 f≥0、f∈L¹(Rⁿ)、W=∫f、λ>0，a>0。全边长 L 的中心立方体概率核记 h_L，u_L=h_L*f，

\[
\Phi_\lambda(t)=t-\lambda\log(1+t/\lambda),\qquad
E_\lambda(L)=\int_{\mathbb R^n}\Phi_\lambda(u_L(x))\,dx.
\tag{1}
\]

令 s=log(L/a)∈[0,log2]。对第 i 个坐标，B_{i,L} 是该坐标在 ±L/2 两端各质量 1/2、其余坐标为 side-L 均匀概率核的卷积算子。每个 B_i 保正且保总质量。真实尺度导数为

\[
\partial_s u_L=\sum_{i=1}^n(B_{i,L}f-u_L),\qquad
\int B_{i,L}f=\int u_L=W.
\tag{2}
\]

先对光滑来源作面微分；一般 L¹ 来源用 L¹ 近似与导数范数 ≤2nW 延拓为 L¹ 绝对连续身份。因为 Φ 为 1-Lipschitz，可取链式微分并用 Fubini：

\[
j_s(x)=\partial_s\Phi_\lambda(u_L(x)),\quad
F_i(s)=\int\frac{u_L(x)}{\lambda+u_L(x)}
       (B_{i,L}f(x)-u_L(x))\,dx,
\tag{3}
\]
\[
E'_\lambda(s)=\int j_s(x)dx=\sum_iF_i(s)
=\lambda\int\frac{\sum_i(u_L-B_{i,L}f)}{\lambda+u_L}\,dx.
\tag{4}
\]

这些身份保留完整 f、原非线性 Φ 与 Lebesgue 接收测度。|F_i|≤W；正部面总质量给已知 V⁺E≤n log2·W。已有[物理尺度稿](physical_scale_bounded_entropy_variation_20261007.md)还给完整正 L¹ 同源输入的 Ω(√nW) 必要下界，以及一个只依赖 ∥f∥₂² 的二次 Fourier 变差界。本轮不重复下界或数值，也不将二次界代替 (4)。

## 2. 标量 V⁺ 为什么尚不能回代 moving 赢家

标量熵正变差是

\[
V^+E_\lambda=\int_0^{\log2}\left[\int j_s(x)dx\right]_+ds.
\tag{5}
\]

它先在全部接收点积分，再取正部。不同接收点选择不同尺度时顺序不再相同。具体取任意可测 ρ_h(x),ρ_s(x)∈[a,2a]，设

\[
H(x,s)=\operatorname{sgn}(\rho_h(x)-\rho_s(x))
1_{\min(\rho_h,\rho_s)\le ae^s<\max(\rho_h,\rho_s)}.
\]

绝对连续链式身份精确给

\[
\int\{\Phi_\lambda(u_{\rho_h(x)}(x))-
          \Phi_\lambda(u_{\rho_s(x)}(x))\}dx
=\int_0^{\log2}\int H(x,s)j_s(x)dxds.
\tag{6}
\]

绝对可积性由 ∫∫|j_s|≤2n log2·W 保证。即使 H 被限制为 0≤H≤1，任意这样的门所需的通量上界也是

\[
\sup_{0\le H\le1}\int\!\int Hj_s
=\int_0^{\log2}\int[j_s(x)]_+dxds
=\int_{\mathbb R^n}V_s^+\{\Phi_\lambda(u_{ae^s}(x))\}\,dx.
\tag{7}
\]

最后一个等式对几乎处处接收点的绝对连续版本成立。右侧是**接收点局部正变差**，一般大于 (5)。不能交换空间积分与正部来用 (5) 支配 (7)。如果门只依赖 s，0≤H(s)≤1，则确实可由 (5) 支付；实际赢家和原来源/历史接受函数不满足这个限制。

当前[hard-average 转接稿](hard_average_transfer_endpoint_20261007.md)式 (3) 的 actual 泛函含 h_{R_h(x)}(x−z)、带原历史门的 Q_C^b(x,y,z,θ) 和原共同来源对。它甚至不是 (6) 的一个已经证明的熵差表示。fullfuture cap 是原自由响应的行约束，不使接受门与 x、y、z 或历史独立。新增角门后的后验也不等于截门前的产品概率。

因此需要先取得一条保留全部实际资格的转接，例如

\[
\mathcal R_{\rm current}\le A_nW+c_nV^+E_\lambda
                         +\varepsilon\lambda|E|,
\tag{8}
\]

其中所有损失须另有同源已付预算；或者直接证明适用于原接受函数的 gated 通量估计。式 (8) **未证**，只是明确缺项。不能由门 ≤1 就把它改成任意 H 的外类、再免费使用 (5)。本稿亦不将 arbitrary-H 外类的失败升级为原实际门的失败。

## 3. 完整正 L¹ 窄盒：局部正变差 ≥c nW

这一解析检查只排除 (7) 的 √n polylog n 一般预算；它不是新的标量熵下界，也不再研究特殊族的渐近赢家。

先取 a=1，任意固定 λ>0，n≥1，定义

\[
\epsilon=\frac1{100n^2},\quad m=4\lambda2^n,\quad
f(x)=m\epsilon^{-n}1_{[-\epsilon/2,\epsilon/2]^n}(x).
\tag{9}
\]

这是完整非负 L¹ 输入，W=m；全程固定源、λ，没有受限来源的再归一化。记 r=2∥x∥∞。对 1+ε≤r≤2−ε 的接收点，side-1 接收 cube 与源盒沿一个最大坐标不交（端点相切只有零 Lebesgue 质量），所以 u_1(x)=0。取 L_x=r+ε≤2，则整个源盒被接收 cube 包含，故

\[
u_{L_x}(x)=\frac{m}{(r+\epsilon)^n}\ge4\lambda.
\tag{10}
\]

对 t≥4λ，log(1+t/λ)/(t/λ) 非增且 log5<2，因此 Φλ(t)≥t/2。点态正变差至少为从 1 到 L_x 的净增量：

\[
V^+_{L\in[1,2]}\{\Phi_\lambda(u_L(x))\}
\ge\Phi_\lambda(u_{L_x}(x))
\ge\frac{m}{2(r+\epsilon)^n}.
\tag{11}
\]

此盒源的响应关于 L>0 为连续分段光滑函数，亦可直接使用绝对连续正变差；单个 L_x 不需要组成一个对所有 x 相同的分割，因为这里只检验局部 V⁺。

full-side radial r 的壳体积为 rⁿ，故 Lebesgue 壳 Jacobian 为 nr^{n−1}dr（没有 2^{−n}）。又 r≥1 且 nε≤1/100，给

\[
\left(\frac r{r+\epsilon}\right)^n
\ge(1+\epsilon)^{-n}\ge e^{-n\epsilon}>\frac12.
\tag{12}
\]

积分 (11) 得完整空间下界

\[
\begin{aligned}
\int V^+_{[1,2]}\{\Phi_\lambda(u_L(x))\}dx
&\ge\frac{mn}{2}\int_{1+\epsilon}^{2-\epsilon}
           \left(\frac r{r+\epsilon}\right)^n\frac{dr}{r}\\
&\ge\frac{mn}{4}\log\frac{2-\epsilon}{1+\epsilon}
\ge\boxed{\frac{\log(3/2)}4\,nW}.
\end{aligned}
\tag{13}
\]

最后一步由 ε≤1/100<1/5，因而 (2−ε)/(1+ε)≥3/2。一般 a 可将输入作 f_a(x)=f(x/a)，两边共同乘 aⁿ；固定 λ 不变。W 随 n 允许增长，这属于“任意完整正 L¹ 输入”的量词，且 (13) 已除以同一 W。

式 (13) 排除 dimension-free、polylog 以及 √n polylog n 的一般 receiver-local V⁺/W 上界；(5) 保留的正负空间取消没有被检验或否定。尤其 (13) 不排除所求标量 √n polylog n 界。没有声称实际接受门实现 H=1_{j_s>0}，没有为这份源补造 CP/GP/FIRST 或原赢家带资格，故也不推出 actual geom 不可付。

## 4. 面平均与 Efron–Stein 尝试的精确欠项

将来源概率记为 μ/W，取 Y∼μ/W、独立 U∼Uniform([-1/2,1/2]ⁿ)，接收点 X=Y+LU。无条件 U 是产品随机量，(3) 可视为对原函数 Φ′λ(u_L(Y+LU)) 的第 i 个 seed 坐标“端点平均减内部平均”，再对完整 Y 积分。

Efron–Stein 可控制内部产品概率下的方差，却不免费控制这个边界 trace：一般有界 seed 函数可以只支撑在宽 δ 的端点条带上，内部方差 O(δ)，两端 trace 为 1。因此不存在从任意 seed 函数的内部方差到端点差的无参数统一不等式。这里的任意函数不是自动可实现的原 Φ′λ(u_L)，所以这只准确定位方法需要补充的原响应 trace 估计，不作为熵上界反例。

进一步，固定接收点 X 后，Y 和 U 的后验由共同完整 f 决定，通常不再产品；原来源对、共同 θ 和实际接受门的条件化不能用无条件 seed 独立性替换。将 hard 来源换成 soft 来源也不能修复 (3) 的面 trace 身份。

作为审计标尺，若能直接证明原非线性、原 Lebesgue 积分下的**集体**估计

\[
\int_0^{\log2}\sum_i|F_i(s)|^2ds
\le C(\log(2+n))^A W^2,
\tag{14}
\]

且 C,A 对 λ、a 和完整 f 一致，则两次 Cauchy–Schwarz 给

\[
V^+E_\lambda\le\sqrt{n\log2}\,
\left(\int\sum_i|F_i|^2ds\right)^{1/2}
\le\sqrt{C\log2}\sqrt n(\log(2+n))^{A/2}W.
\tag{15}
\]

式 (14) **没有证明，也没有登记为成立的新阶数候选**。逐坐标 |F_i|≤W 只给右侧 n log2·W²，回到线性阶；逐坐标各领取 W 不是 source-once 集体预算。当前没有能从 Efron–Stein 消去这额外 n 的合法 trace 能量证书。即使 (14) 成立，§2 的实际桥接仍需独立证明。

## 5. 本轮交付范围

已核清原 Φ 面通量身份、moving 端点的空间权重身份、任意门与 local 正部的精确关系，并给 (13) 的完整原 L¹ 解析论证。未新增一般 √n 上界，未修改主账，未付 current actual 余项。

本轮无新增实验：没有已论证可行的新分层/阶数机制需要预注册数值判断；窄盒检查有完整解析证明，重复特殊族数值不能补齐原实际截门的转接。下一项真正需要的数学输入是保留实际门的 source-once 转接/commutator 预算；不能将 (13) 所排除的任意门外类作为替代合同。
