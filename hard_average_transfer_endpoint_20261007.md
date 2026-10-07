# Moving hard average 的实际 endpoint 转接：原泛函、signed collar 与截门损失

2026-10-07。只新增本稿；不修改总账、不重跑已有数值。本轮没有取得新的 actual W 空间付款。可核成果是把需要的转接写成原标签上的精确泛函，并把两项可能误用的取消具体展开。新 frozen 移动物理尺度 Gaussian 网沿用已经核查的工具，不重复其证明或有限矩共同目标障碍。

## 1. 本轮输入与确切待估泛函

原 tex 5522–5585 给早首跳的单份来源测度与完整 continuation，6096–6204 给规范 contact 表示及它未付的重返项，4873–5070 为 NB--L/M/S 的准确方向。读取 `gated_occupancy_bridge.md` §1–2、`actual_gate_contract_audit_20261007.md` §2、§5，及当前 joint 账。J05 已在上轮读其 SKILL 和必要 references；本轮仍按其真实强度、原标签、完整后续和可预测性合同工作，不补造新停止规则。

原 f≥0、μ=f dx、W=∫f 是完整输入；原分配 ν_C=fη_C dx、Σ_Cη_C≤1，硬源 ν_h=μ_hi|出生集合。原可测 σ(x),L_s(x),R_h(x)∈原允许目录，L_s,R_h∈[a,b]、a≤b≤2a，q∈[3λ/4,λ]，u=h_Rh*μ∈(2λ,4λ]。原完整 FIRST 与所有 future s≥σ,L∈[a,b] 保持
\[
P_{\sigma,L_s}μ(x)=q,\qquad P_{s,L}μ(x)≤q.
\tag{1}
\]
全部 GOOD/score/owntrace、出生、原格、共享 θ 与唯一 distinct-child LCA、失败币、时钟、strict CP/GP、远对和重捕获保留；当前 R_dagger 的 lowS、lowcoin、诊断双 cutoff、来源路径核心补集也保留。辅助 Bernstein K 不等于实际 continuation 活跃数。

为容纳依实际 i,t,r,w 的门，记 b_C(x,y,z,θ;t,i,r,history)∈[0,1] 为这些门的完整接受函数，已经平均未来细历史时使用其 backward projection，不能把未来事件当可预测指示。定义带门的原子核
\[
\begin{split}
Q_C^b(x,y,z,θ)
={}&\sum_i\int_0^{\min(\sigma(x),T_n,1/2)}(1-t)^{n-1}\int
1_{y+r e_i\notin C}\,G_{1-t,L_s(x)}(r)\\
&\quad\times [h_{L_s(x)}^{\otimes n}*K_{\sigma(x),t,L_s(x)}]^b(x-y-r e_i),drdt,
\end{split}
\tag{2}
\]
其中上标 b 表示完整路径门的正子核，不是删掉 continuation；故 0≤Q_C^b≤Q_C^early≤P_{σ,L_s}。若所有额外门只依端点，则可将 b 直接乘在标准 continuation 外。原目标按同一 θ,LCA 先合并后为
\[
\boxed{\mathcal R_\dagger
=\mathbb E_θ\int_E\int\int h_{R_h(x)}(x-z)\,
q^{-1}\sum_C Q_C^b(x,y,z,θ)\,dν_C(y)dν_h(z)dx.}
\tag{3}
\]
若原接口仅给正支配，使用 ≤ 而非等式。C 来源分配及 b 不得复制；原门在 y,z 上，不能改判于 w=y+r e_i。完整 Q/P 才是≤1 的接受比，固定 t,r 的密度比不是概率。式 (3) 与现有 source-once g 形式相同，只把本轮需要的 history 重新显式恢复。

删掉硬行和其它门可以用完整 u≤C_hq，得到 actual 子交通的软空间支配；但必须先有一份与 L_s(x) 无关的带初始标签正参考退出测度。NB--L 对固定参考 live 源的轴向退出常数2提供这一步的特定方向，不给不同 L 的各自 live 演化共同 W。早期首跳参考源总质量≤2W，完整 continuation 仍需支付；这不是新的小源。

## 2. 全 hard-average maximal 是更强问题，不能免费回代

固定 c，定义使用共同 ν 的
\[
\mathfrak M_{hH}^{c}ν(x)
=\sup_{a≤L≤b,\tau>0}(h_L^{\otimes n}*H_\tau^{c,L}ν)(x).
\tag{4}
\]
完整 frozen H 含总率 n/c 的 holding 分量，因此逐点
\[
h_L*H_\tau^{c,L}ν≥e^{-n\tau/c}h_L*ν,
\qquad
\boxed{\mathfrak M_{hH}^{c}ν≥M_{cube,[a,b]}ν.}
\tag{5}
\]
取 τ↓0 即得，既不依赖 source atom 也不要求 ν 产品。故若直接提出一般 (4) 的 √npolylog 弱预算，它至少包含原共同窗口 cube 的一般预算；若提出“不添费用继承 O(log n)”，更不能由 frozen 定理逻辑推出。式 (5) 是转接强度检验，不构造满足全部实际门的反例。

实际 (3) 的首跳退出与原 fullfuture、远对、CP/GP 等门可能缩小这个更强外类，所以不能拿 (5) 宣称 actual 不可付。反过来，放宽到 (4) 后再以原门已删除为由忽略其中的硬 maximal，也不能成立。

已有 Gaussian 1/√n 尺度网只控制 H_τ^{c,L}ν。h_L 移动使每个节点输入变成 h_L*ν；每份质量虽为 W，但不是共同函数。给定 uniform seed v 后，Gaussian 中心是 Lv，随尺度移动，不能用 centered Gaussian 网。普通正硬核比较需要 1/n 尺度格才保常数，这是已核物理接口，不重复计算。

## 3. 保留 signed collar 的精确障碍分解

这里只在一个固定 (c,L) 下使用已有障碍，而非让不同参数免费共用它。记正 odometer U、域 Ω={U>0}、σ_fr=L_fr U=ν_b−μ_b、ν=ν_b+ν_g、μ_total=μ_b+ν_g≤κ。定义
\[
R_\tau^{fr}\sigma_{fr}=H_\tau\sigma_{fr}
-\tau^{-1}\int_0^\tau H_s\sigma_{fr}\,ds.
\]
卷积可交换且 U∈L¹，故严格有
\[
\boxed{h_L*H_\tauν
=h_L*H_\tauμ_{total}
 +\frac{h_L*U-h_L*H_\tau U}{\tau}
 +h_L*R_\tau^{fr}\sigma_{fr}.}
\tag{6}
\]
第一项满足 \(h_L*H_\tauμ_{total}≤κ\)，因为 \(μ_{total}≤κ\) 且两算子均正、保常数。这里使用 \(ν=μ_{total}+σ_{fr}\)，并不假定 \(H_\tauμ_{total}=μ_{total}\)。平均项与 R 项相加恰为 \(h_L*H_\tauσ_{fr}\)。
不加 h 时在 Ω^c 的中项≤0。加 h 后，只有在 (h_L*U)(x)=0 时能据正性继续去掉该中项；例如 Q(x,L) 与 Ω 不交时成立。一般 x∉Ω 不够，因为 x 的硬方体可以触及 Ω。这不是空间边界零测争议，而是正质量平均。

若用全部 collar Ω+Q_b 支付，现有 κ|Ω|≤W 不控制 κ|Ω+Q_b|。其精确源成本是所选障碍支持的平行体，不能把 |Ω| 原值沿用；当前 trimmed 来源核心已经支付满足指定体积证书的支，不可因此把任意障碍 collar 也算作已付。原 U 的有限占用仅用于定义，∫U 没有已证明的 O(W) 常数，不能将 (h_L*U)/τ 单独计成新正源费用。

对于固定参数，保留取消可写
\[
\frac{h_L*U-h_L*H_\tau U}{\tau}
=h_L*\frac1\tau\int_0^\tau H_s\sigma_{fr}\,ds.
\tag{7}
\]
这是合法 signed 身份，但把右边换成 σ_fr,+ 后，虽有 σ_fr,+≤ν，得到的正上界又含 moving hard-average maximal；没有从 (7) 消除 (5) 的困难。更不能为不同 L 单独重建 U_L 后声称所有 σ_fr,L,+ 合计仍≤ν：只逐 L≤ν 不给 source-once 总和。

## 4. Fullfuture 的尺度符号不能直接继承到截门交通

在原 selected σ 和普通可微物理内点 L_s，(1) 给完整响应的局部最大值；因此
\[
\partial_{\log L}(P_{\sigma,L}μ)(x)|_{L_s}=0,
\tag{8}
\]
端点或尖点只给相应单侧差分/导数符号，不能假定零。为不丢 hard face，先对任一合法 future 物理 L' 用有限差分
\[
D_x(y)=P_{\sigma,L'}(x-y)-P_{\sigma,L_s}(x-y),
\quad\int D_x(y)dμ(y)≤0.
\tag{9}
\]
对固定 x,z,θ，原 actual 接受 a_x(y)∈[0,1] 不必独立或对参数光滑；精确地
\[
\begin{split}
\int a_xD_x\,dμ
&=\int D_xdμ-\int(1-a_x)D_xdμ\\
&≤\int(1-a_x)(D_x)_-dμ
≤\int(D_x)_-dμ.
\end{split}
\tag{10}
\]
所以完整差分非正不保证 selected 差分非正。要用尺度取消支付 (3)，必须控制被原门排除的来源所承担的负部，或给接受门的可用 commutator 预算；仅有 0≤a≤1、唯一 LCA 和 lowcoin 不给该控制。这里没有构造新的 abstract toy，也没有将它作为 actual 反例。

对式 (2) 的 history 子核，参数变化还会改变 G_{1-t,L}、hard seed 和所有以后跳核；若冻结原 b 来比两个参数，只得原 test-family 差分。它不是重新计算原 FIRST/history 的差分。原 fullfuture 是原 μ 的自由响应约束，不能从中给重新条件化的退出/停止源继承 q。

## 5. 规范 ρ 与共同物理尺度：不丢重返的安全弱式

原 tex 固定 L，原 U_uniform 在整个路径中只抽一次。规范 ρ_0 与尚未接触的到达 ρ_v 保留 Y,C,U_uniform、首跳状态和跳标；holding 接触测度 B_L 质量≤W。原 selected 图像密度则还包括首次停止后严格 v<γ(x) 的完整增广传播。因而
\[
A_{C,L}(γ(x);x,y)=B_L^{eligible}\text{ 的规范图像密度}
 +\text{此前接触后仍传播至原 }x\text{ 的非负项}
\tag{11}
\]
应理解为该已有 augmented-history 分解，不把任意测度在图上求普通点值。不把 clock-contact 在 K_{s,s} 的恒等原子上重计，不把原 sourcepair 门转到停止点。

固定 L 的 B_L 部分可直接用质量支付其保持的实际正接受函数≤1；剩余重返项仍需图像 occupation 或统一逃逸。将每个物理 L 的 B_L 逐个相加复制 W；随机先选 L_j 的接触律是 p_jB_{L_j}，恢复 selected 输出要付逆权。停止后输入已经条件化原 uniform seed，与重新独立抽 h_L 的输入类不同，不能把 K*κ 免费写成 h_L*K*κ。这些弱式来自旧 actual contact，不冒称本轮新支付。

## 6. 本轮判断与下一条具体所欠估计

尚未找到对全部当前 actual R_dagger 的新可收费分层。本轮无新的维数阶数候选，故不注册抽象数值、不重跑旧下界、contact、frozen 网或 core 守卫。下界目录 A/B 的任意联合正权与增长源标签继续允许；本稿没有以峰型数或源深度给 hard 平均外类降维。

有意义且未证的下一接口，须在原 (3) 内给出下列三种之一的**实际空间估计**：保留符号的 hard collar 高频/平均正部预算；fullfuture 差分的 actual 截门负部/commutator 预算；规范原 stopping-history 的重返 occupation 预算，并覆盖共同移动 L 的参数迹。其费用须由同一原 μ 的 source-once 总质量给出，不能重启 FIRST、不能以逐输出捕获源替代输入。仅证明完整行导数非正、再列一个 lowK 或源质量条件，尚未做到这一步。

当前主账仍保留 R_dagger；本稿没有新增已付费，也没有证明一般 √n n^o(1) 目标。
