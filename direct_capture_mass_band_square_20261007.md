# 完整捕获质量 band 的直接平方函数：无需 packets

2026-10-07。独审父任务提出的直接表示。结论严格成立：把 **原完整捕获质量** 分 dyadic bands，每个来源的单 band profile≤2，全部 band 的平方函数≤2I₁。它无需 packet、fullness 或 winner 极大性。跨 band 总 profile 仍未付；原幅度带和 bounded window 只恢复 O(n) 基线，没有平方根结论。

沿用已读 A01/E01 的完整来源与 occupation 范围守卫；具体结论完全由捕获距离、体积和 Tonelli 自证，不调用外部定理。未运行数值或改主总账。

## 1. 查重与正确简化范围

已全文核对 [gated occupancy §1–7](gated_occupancy_bridge.md)、[single capture §1–5](single_capture_source_rearrangement_20261007.md)，检索本课题目录的 capture-mass/mass-band/质量分层对象；另核 [stopped diffusion supplement](../anisotropic_cube_flow/stopped_diffusion_supplement.tex.txt) 的“单簇捕获份额”段。

旧 gated §3 是硬核 cube-cone 径向身份及全窗 O(n)W；§6 是 complete mass≤W 导出的最大半径与固定来源窄尺度窗；§7 是有限硬原子的份额阈值及 O(log N) 交通费。它们都使用相同的 capture/半径/volume 原理，但分层变量与本稿完整 m(x) 不同，未找到其写出本稿 (5)–(7) 的无 packet 来源平方函数。旧 single capture §4 的 exact-capture-set diagonal 相关，但有 background βR^n+m_A 分母及 spike 来源。旧 supplement 的固定 packet 半径 cap 加 Tonelli 也已包含同一 elementary 机制，不能把这种机制重新宣传为新的 cube theorem。

因此本稿是对现有路线的 **直接简化和接口整理**：若只想要完整原 profile 的 band-diagonal，使用完整 m(x) 足够，不需要先固定 packet、验证 fullness、再分 share band。新 [packet share-band 稿](general_packet_share_band_square_20261007.md) 的内容仍适用于 packet-restricted 子来源与真实 share；本稿保留所有来源，因此处理的是更简单且不同的变量。

## 2. 原来源、原 receiver 与 mass bands

设 μ 是完整有限正 Borel 测度，W=μ(R^n)>0；τ>0。给原可测 receiver E 和原可测 full-side R(x)>0，令

\[
 Q_x=x+[-R(x)/2,R(x)/2]^n,
 \quad m(x)=\mu(Q_x),\quad m(x)>\tau R(x)^n,
 \quad I_1=\tau|E|<\infty. \tag{1}
\]

这包含真实 finite winner 或合法 continuous selector 下的 strict 原事件；证明本身只用 (1)，不会改选 R。闭 cube 规范、原子来源均可保留。μ 有限保证 0<m(x)≤W，无无限质量 band。

任取 B>0，并令 k∈Z，

\[
 A_k=2^kB,\qquad
 D_k=\{x\in E:A_k<m(x)\le2A_k\},
 \quad I_k=\tau|D_k|. \tag{2}
\]

下端开、上端闭；m=2^kB 属于 k−1，所有正 m 恰属一个 band。D_k 是原 receiver 的划分，并非按来源 y 重组选取的输出。令

\[
 T_k(y)=\tau\int_{D_k}\frac{\mathbf1_{Q_x}(y)}{m(x)}dx,
 \qquad S(y)=\sum_kT_k(y)
       =\tau\int_E\frac{\mathbf1_{Q_x}(y)}{m(x)}dx. \tag{3}
\]

所有 countable 求和以非负 Tonelli 理解。B 仅选择质量档位的 origin，不增加源质量或改变 incidence。

## 3. 点态 single-band cap 与平方函数

对 x∈D_k，由 (1)–(2)，

\[
 R(x)^n<\frac{m(x)}\tau\le\frac{2A_k}\tau,
 \qquad\frac\tau{m(x)}<\frac\tau{A_k}. \tag{4}
\]

若同一原来源 y 被捕获，则 ||x−y||∞≤R(x)/2，故所有这样的 x 都在固定来源中心的 cube

\[
 Q\left(y,(2A_k/\tau)^{1/n}\right),
\]

其 Lebesgue 体积是2A_k/τ。因此直接得到

\[
 \boxed{0\le T_k(y)\le
       (\tau/A_k)(2A_k/\tau)=2.} \tag{5}
\]

这个 cube 只上包已经发生的原捕获 receiver，不把原 R 替成一个新查询。它不依 n、source support、atom count、density cap、packet 几何或 fullness。

同一完整 μ 下 first moment **准确**为

\[
 \int T_k(y)d\mu(y)
 =\tau\int_{D_k}\frac{\mu(Q_x)}{m(x)}dx
 =I_k. \tag{6}
\]

所以 T_k²≤2T_k 与 Σ_kI_k=I₁ 给

\[
 \boxed{\int\sum_kT_k(y)^2d\mu(y)
       =\sum_k\int T_k^2d\mu\le2I_1.} \tag{7}
\]

没有要求 μ∈L²(dx)，也没有把 source μ 冒当 product measure。这是 dimension-free **diagonal/square-function** 预算，其右侧是仍待证明弱输出量 I₁，不是已付 O(W) 交通费。

可保留一个无需额外工具的较紧常数，但并非主改进：对 V=(2||x−y||∞)^n，捕获给 V≤2A_k/τ，且 τ/m≤min(τ/A_k,1/V)。利用 |{V≤v}|=v，

\[
 T_k(y)\le\int_0^{2A_k/\tau}
                   \min(\tau/A_k,1/v)dv
       =1+\log2. \tag{8}
\]

故 (7) 的2可换成1+log2。一般质量比 ρ>1 的同一证明给 single-band cap≤1+logρ；仅用矩形 coefficient×volume 给ρ。本文保留 (5)–(7) 的常数2作为最简接口，不另为常数优化启动数值。

## 4. 原完整 S 的跨 mass-band 交叉仍须几何

精确 source square 为

\[
 \int S^2d\mu
 =\sum_k\int T_k^2d\mu
       +2\sum_{k<l}\int T_kT_l\,d\mu. \tag{9}
\]

D_k 不交不意味着 T_k 的来源支撑不交。原 y 可被不同质量档的许多 receiver 捕获；不能把 (7) 写成 ∫S²≤2I₁。

令有正 I_k 的实际 Gram

\[
 C_{kl}=\frac{\int T_kT_l\,d\mu}{\sqrt{I_kI_l}},
 \qquad\alpha_k=\sqrt{I_k/I_1}.
\]

I_k=0 的 band μ-a.e. 为零，可删除。在有限截断上 C 是真实 source Gram、C_kk≤2、||α||²=1，而

\[
 \frac{\int S^2d\mu}{I_1}=\alpha^TC\alpha \tag{10}
\]

按非负截断极限理解。此 PSD 仅属于完整 source profiles，不能把它移植到带 receiver-dependent LCA 远对门的另一核。bounded diagonal 本身不控制跨档 Rayleigh 值。若另证 (10)≤K_n，则原 ∫S dμ=I₁ 与 Cauchy 才给 I₁≤K_nW。将希望的 K_n=O(sqrt(n)n^{o(1)}) 再写一遍，不是证明它。

若仅有 L 个非零质量 bands，则 Cauchy 给

\[
 \int S^2d\mu\le L\sum_k\int T_k^2d\mu
                         \le2LI_1. \tag{11}
\]

任意 μ、任意 R 下质量可趋零，L 可无限；即使窗口 R≥a 且 (1) 给 m>τa^n，若没有原响应上幅度带，质量范围 W/(τa^n) 也可任意大。

## 5. 原 amplitude band/window 恢复 n 基线

若进一步保留原幅度 τ<m(x)/R(x)^n≤2τ、R(x)∈[a,2a]，则

\[
 \tau a^n<m(x)\le2\tau(2a)^n=2^{n+1}\tau a^n.
\]

选择 B=τa^n 后，只可能有 k=0,…,n 的 n+1 个档。由 (11)，

\[
 \int S^2d\mu\le2(n+1)I_1,
 \qquad I_1\le2(n+1)W. \tag{12}
\]

这仅为已知 O(n) bounded-window 基线的另一写法。旧 positive hard envelope 给更直接的 S(y)≤1+nlog2，从而 ∫S²≤(1+nlog2)I₁；不能把 (12) 登记为更好的已有费用。需要省掉的仍是同-source 跨 mass-band 的 n，而不是证明 packet coverage 后再创造 n 档。

**粗化档位不能假省平方根。** 把质量 band 的比值改为 ρ>1，令 d=logρ，(8) 的同一径向 max-kernel 积分给每块 T≤1+d，全部 block 平方函数≤(1+d)I₁。当前质量跨度是 Δ=(n+1)log2，一般覆盖这个跨度需 L=ceil(Δ/d) 块。再用跨块 Cauchy 的系数是

\[
 (1+d)\lceil\Delta/d\rceil. \tag{13}
\]

该表达式≥Δ；d≈sqrt(n) 时，每块成本≈sqrt(n)、块数≈sqrt(n)，乘积仍 O(n)。L=1 时也需 d≥Δ，single-block 费用本身已为 O(n)。实际输入若只占少数档当然可以省，但那须证明实际质量占用限制，不能从改分块获得。这里没有新增阶数实验或用常数优化替代 geometry。

与 packet share-band 的区别：后者固定 i 的 band 起点随 M_i 平移，需 fullness 才取得 p_i 的上下端；本稿所有 receiver 直接按同一个完整 m(x) 分档，既不付 η^(-1) 也不需要 source 分配，但没有定位哪个 packet 或哪个 actual gated 子来源产生跨档合作。

## 6. 实际门可以保留，付款校准不能省略

若真实 coefficient a(x,y)∈[0,1]，可定义同一原 D_k 上

\[
 T_k^a(y)=\tau\int_{D_k}\frac{\mathbf1_{Q_x}(y)}{m(x)}a(x,y)dx.
\]

仍有 T_k^a≤T_k≤2，且 Σ_k∫(T_k^a)²dμ≤2Σ_k∫T_k^a dμ≤2I₁。这里保留的 FIRST/fullfuture/CPGP/LCA/history gate 必须是原真实已经积分后系数，不是一个任意事件代理。其 first moment 一般 **小于** I₁，若缺原阈值下界，不能从它自动恢复 τ|E|。换成固定硬子来源 ν_h≤μ 的版本也只有相应的受限 first moment，不能把它重新校准成完整 I₁。

本稿因此提供一个更简洁的一般 source-square 起点及明确 cross-mass defect，不新增 actual geom paid fee。未来真正要证的是保留原门的跨 mass-band 控制，或更弱的原来源截断-tail 合同；完整 massband diagonal 已无需繁复 packet 细化。没有数值重跑、没有把 (7) 或 Gram 存在误报为平方根预算。
