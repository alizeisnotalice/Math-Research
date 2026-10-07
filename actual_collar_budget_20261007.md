# 实际重返 occupation 的一份来源：初始接触候选与 early 资格降级

2026-10-07。只新增本稿，不修改主账。本轮直接研究 `hard_average_transfer_endpoint_20261007.md` (3) 的实际 history。选取最可能给固定来源的重返项；没有得到新的困难 actual 子支付款。曾形成一项来源预算候选，随后在原 early 截止处严格降级为空支／旧 late-first-jump 工具，故不运行组成 toy 数值。

## 1. 原实际对象及工具合同

保留同一完整 f、W、ν_C=fη_C、Ση_C≤1，原 y,z、细格、出生、共享 θ 与唯一 LCA、失败币、重捕获、strict CP/GP、FIRST fullfuture、owntrace/score/时钟和当前 R_dagger 的全部门。原第一跳时间 τ₁<T_n=(12+log(n+2))/n，第一跳必须从原 C(y) 退出，且原 continuation 完整保留。原接受函数在初始 y,z 与实际历史上，不转到观察位置或重建 FIRST。

沿已读 J05 的真实强度和补偿合同，并使用已证明 hazard-anchor 的正比较；不把技能当来源预算定理。读取原 tex 6096–6204 的完整 contact 与重返定义、NB 4873–5070、early/late 5522–5585，及 hard-average 稿的确切目标。下界总表 A/B 的任意正联合权、增长标签、XOR/Cantor/Gibbs 继续允许；没有产品源、有限标签或源深度假设。

固定物理 L，抽一次 U_L∈Q_L，初始 Y~ν_C，原物理状态 X_s=Y+Z_s，观察 O_s=X_s+U_L。首跳资格按 X_s，接触时钟按 O_s；U_L 在整个历史只抽一次。

## 2. 改造接触图但不改原终端资格

令 E_0 为当前待研究的实际输出子集，原 σ(x)≤S<1（可先作此有限窗口再取极限），且 σ(x)≥s_0>0。把 σ 在全空间作固定 Borel 延拓，再定义
\[
\gamma(v)=\min\{S,\max(s_0,\sigma_{ext}(v))\}.
\tag{1}
\]
在 E_0 上 γ(x)=原σ(x)，故规范首次接触图不改变该支的终端观察时刻；全空间 γ≥s_0。γ 只用来分解原自由历史，未替代 FIRST、fullfuture 或实际门。接触后原状态和全部后续传播保留。

将该支按首次接触是否在任何首跳之前拆分。首跳前 O_s 恒为 v=Y+U_L；因此“首次接触在首跳前”恰为初始位置的 clock contact 时间 γ(v)，且必须 γ(v)<τ₁。其带原 Y,C 的正源为
\[
d\kappa_{0,L}(v,t,y,C)
=(1-\gamma(v))^n h_L(v-y)\,dν_C(y)\,dv\,\delta_{\gamma(v)}(dt).
\tag{2}
\]
没有自由选择图像密度版本的问题：这是 no-jump 路径的规范 Jacobian。time 条件源一般不是 Bochner L¹(dt;L¹(dx))，不能直接引用只处理该输入类的 hazard 公式。

但这一个源确有共同物理尺度支配：以 S_h(v-y)=sup_{L∈[a,b]}h_L(v-y) 代替 h_L，保留同一 y,C 与 time=γ(v)；已有 cube-cone 恒等式给 ∫S_h=N_h=1+nlog(b/a)。于是单份 joint 参考源满足
\[
\|\overline\kappa_0\|≤(1-s_0)^nN_hW.
\tag{3}
\]
这一步没有将停止点当新原输入，实际门仍检查 y,z；只是正参考的接触来源预算。

停止后观察位置的增量是原 K_{s,t,L}，**不再另加 h_L**：U_L 已在 v 中。重新独立平均 U 会改变条件停止源，禁止如此操作。

## 3. 无额外 hard-average 的部分面包络与一个有效但昂贵的源费

固定 source 时刻 t<1，K_{s,t,L}=[rδ₀+(1−r)G_{1-t,L}]^{⊗n}，r=(1−s)/(1−t)。保留所有 masks A；每个 A 的正核在其 |A| 维坐标面上为偶、逐坐标递减的概率密度，其物理尺度 envelope 质量≤1+|A|log(b/a)。holding A=∅ 是共同 δ₀，质量1。

以原 Bernoulli 峰值 c_|A| 控制所有 s 后，得到一个固定 t 的正**测度**包络 Λ_t，
\[
K_{s,t,L}≤Λ_t\quad(t≤s≤1,a≤L≤b),\qquad
\|Λ_t\|≤N_hD_n,
\quad D_n≤1+(\pi/2)\sqrt n.
\tag{4}
\]
不把部分面当 full bounded density，不对 t 再免费取 supremum。对 (2) 的 time graph，可按空间 v 的 L¹ marginal disintegrate 出 time/原来源标签；部分面卷积使用该规范切片。式 (4) 是逐 t 的 kernel-measure 正支配，随后 Tonelli 总质量预算成立。holding 传播不能被指定普通 delta 点值；对于严格重返 t<γ(x)，同位置 holding 的接触图贡献实际为零，正上包仍可保留它。

删门只在带原标签的正支配后，以原硬行 C_hq 抵消一次 q，给出候选分支的真实上界
\[
\mathcal R_{initial-contact,return}^{act}(E_0)
≤C_h(1-s_0)^nN_h^2D_nW.
\tag{5}
\]
这是两个不同物理接口的真实费用：初始 uniform seed 的来源 envelope N_h，后续 noise 的 partial-face envelope N_hD_n。不能隐藏其中任何一个，也没有用 weak-L¹ 可加性。

若进一步恢复“接触后必须首跳早于 T_n”，可先以 G_{1-t,L}≤2G_{1-t,b} 支配这个首跳的单轴 mark，保留其初始 y,C，再取后续 (4)。在 s_0<T_n 时，参考来源 mass≤
\[
2N_hW[(1-s_0)^n-(1-T_n)^n],
\]
从而
\[
\boxed{\mathcal R_{initial-contact,return}^{act}(E_0)
≤2C_hN_h^2D_n[(1-s_0)^n-(1-T_n)^n]W.}
\tag{6}
\]
概率差来自原 no-jump survival；不是把 endpoint 后验当独立 first-jump 先验。

## 4. 原资格使表面平方根候选降级

式 (5) 的朴素选择 s_0=2log(n+2)/n 会给 (1−s_0)^n≤(n+2)⁻²，形式上 N_h²D_n=O(n^{5/2}) 变成 O(√n)。然而当 log(n+2)≥12，s_0≥T_n：首次接触时间 γ(v)≥s_0 后才发生 τ₁，必有 τ₁≥T_n，违反实际 early。该支此时为空，不能登记为非空的一般空间进展。

保持 s_0≤T_n，且只按最大允许的 no-jump 衰减估算，e^{-nT_n}=e^{-12}/(n+2) 仅是一阶 n⁻¹；两个 N_h 与 D_n 仍留下约 n^{3/2}W。不能借“首跳在该时钟之后”绕过原 early 截止。

这里两种费用的薄窗阶数必须区分。令 Δ=T_n−s_0≥0。固定常数范围 nΔ=O(1) 下，均值定理给
\[
(1-s_0)^n-(1-T_n)^n
≤nΔ(1-s_0)^{n-1}=O(Δ),
\]
因为 \((1-s_0)^{n-1}=O(n^{-1})\)。因此式 (6) 含 \(N_h^2D_n=O(n^{5/2})\)，Δ=O(1/n) 一般只给 \(O(n^{3/2})W\)，不能宣称平方根；必须取 Δ=O(1/n²) 才由这个上界得到 \(O(\sqrt n)W\)。对于 Δ=C/n、固定 C>0，生存差实际上为 Θ(n⁻¹)，所以这里也不能凭一般符号 O(1/n) 暗藏额外一阶增益。

另一方面，此候选支始终蕴含原 τ₁≥s_0。无需 contact 分解，原 FJL 的单轴共同退出源只有一个 N_h，直接给
\[
\mathcal R_{\tau_1≥s_0}^{act}
≤2C_hN_hD_n(1-s_0)^nW
=O(\sqrt n)W
\quad(s_0=T_n−O(1/n)).
\tag{7}
\]
式 (7) 在 Δ=O(1/n) 即成立，因 \(N_hD_n=O(n^{3/2})\) 且 \((1-s_0)^n=O(n^{-1})\)。它比 (6) 更便宜，且只是原 late-first-jump 小来源公式向左延伸常数 hazard 长度。故本轮不把 (6) 的更薄特例称新 occupation 支、不重复领取旧费用、不为它做三轮测试。

## 5. 已首跳后的接触源仍是具体剩余障碍

其余首次接触已在至少一次跳跃之后。pending/ineligible/eligible 状态需按原首跳判定，不能只用接触位置重新判格。对应来源由 survivor arrival ρ_{v,L} 与 jump-entry 构成；其 marginal 固定 L 的 mass≤W，但 L 改变了此前传播、uniform seed、未接触条件及条件原来源。没有得到一份共同 \(\overline\kappa(dt,dv,dy,C,history)\) 且 mass=O(W/√n) 或其他足够小预算的正支配。

NB--L 只比较同一 live 源的一次轴向退出，不能以它支配这些不同尺度的 live 演化；原 fullfuture cap 测的是同一原 μ 的自由传播，不给条件停止源继承 q。若只保留各 L 的 mass≤W、重新领一份 ν_L 后使用 (4)，依然复制来源预算。若用 (1) 的 γ≥s_0 对已首跳接触也乘 (1−s_0)^n，则错误：该因子只属于初始 no-jump 接触，已跳接触的最后保持期是 [(1−γ(v))/(1−a)]^n，a 为实际到达时间，可以任意接近 γ(v)。

这准确定位了未能完成的交换，而不是只给一个新行条件：需要控制实际 survivor-arrival 的 moving-L **联合来源**，或其原 history 重返 occupation；不能把前述 initial-contact mass 因子扩展给它。目前没有该空间证明。

本轮不新增主账费；当前 actual R_dagger 保持。候选在原资格下已降级，无实质未验证的新阶数策略，因此没有预注册或运行 abstract toy/系数数值。目标一般 √n n^o(1) 仍未闭合。
