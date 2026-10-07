# 共同实际赢家的正缺额交换：精确接口与空间收费失败处

2026-10-07。独立解析尝试；只新增本稿，不改联合账、历史数据或别人文件。结论：得到一个保留原实际交通、同一完整来源且不添森林层数的正交换恒等式，但没有新的平方根空间费用。另一条利用 strict CP 小捕获比例的“捕获—未捕获”交换也有精确表达；其小比例换来了失去来源中心支撑，不能按旧 hard column 结算。两条路线均未覆盖支付 R_dagger。

## 1. 冻结对象与本轮边界

已读最新 `current_joint_budget_20261007.md` 的同价自适应核心更新，以及 `hard_average_transfer_endpoint_20261007.md` 全文；另核对 `gated_occupancy_bridge.md`、`actual_gate_contract_audit_20261007.md`、`actual_parent_occupancy_route_20261007.md`、`far_source_geometry_bridge_20261007.md` 和 `auxiliary_exit_column_bridge_20261007.md`。本稿不用 ordered/frozen hard-average、障碍 collar 或已付 fixed-delay inner 分支。这里只用原共同响应、正 Tonelli 和显式正耦合；不调用外部输运存在性定理。

固定完整 μ=f dx、W=μ(R^n)，原共同 q、FIRST σ(x),L_s(x)，原 R_h(x)∈允许目录 J⊂[a,b]，原事件 E 与 hardband 2λ<u(x)≤4λ，并显式保留 **3λ/4≤q≤λ**。该阈值关系直接来自原提示词《中心立方体_geom余项最终估计_Pro提示词_20261006.txt》的“冻结实际对象”段（本轮重新核对第 101–105 行），也与 `actual_gate_contract_audit_20261007.md` §2 的共同 FIRST 合同一致。故 q≤λ<u/2；本稿不足二倍行放大依赖这一关系，不能仅由 hardband 推出。完整 future 合同仍为 P_{σ,L_s}μ=q、P_{s,L}μ≤q 对所有合法 future s≥σ,L∈[a,b]。保留完整 early/history/continuation、初始 y 标签、出生、唯一 LCA、失败币、strict CP/GP、far、重捕获，以及最新 lowS、lowcoin、K 双 cutoff 和 adaptive-core 补集。参考来源不重启 FIRST，不重新归一化。

在原 source-once 合并之后，记当前真实接受权为 g_dagger,θ(x,z)∈[0,1]，并设

\[
 d\tau_{x,\theta}(z)=h_{R_h(x)}(x-z)g_{\dagger,\theta}(x,z)d\nu_h(z),
 \quad T_\theta(x)=\tau_{x,\theta}(\mathbb R^n)\le u(x).
 \tag{1}
\]

这里 g 已包含全部原软端积分、实际历史和所有余项门；不是 arbitrary g 的新模型。E 外取 T=0，定义表示交通 R_τ=E_θ∫_E T_θ dx。若原合并接口精确，则 R_dagger=R_τ；若上游只给正支配，则仅有 R_dagger≤R_τ。后续耦合恒等式对 R_τ 精确成立，不能反称原余项等式。

## 2. 同一原输入上的正缺额交换

把原 hard 尺度 R=R_h(x) 用作辅助 future 物理参数，保持原 σ，不改变 L_s。它属于合法 [a,b]，故

\[
 F_R(x)=P_{\sigma(x),R}*\mu(x)\le q.
 \tag{2}
\]

定义固定原输入上的参考正测度

\[
 d\xi_x(z')=[h_R(x-z')-P_{\sigma,R}(x-z')]_+d\mu(z'),
 \quad D(x)=\xi_x(\mathbb R^n).
 \tag{3}
\]

由正负部分与 (2)，

\[
 D(x)\ge u(x)-F_R(x)\ge u(x)-q>u(x)/2>0.
 \tag{4}
\]

现在显式给出交换，不要求被实际 g 选择的来源自身承担同样的 deficit：

\[
 d\Pi_{x,\theta}(z,z')=d\tau_{x,\theta}(z)\,D(x)^{-1}d\xi_x(z').
 \tag{5}
\]

它的第一边缘恰为原 τ，第二边缘恰为 (T_θ/D)ξ，且

\[
 0\le T_\theta(x)/D(x)\le u(x)/(u(x)-q)<2.
 \tag{6}
\]

因此原实际接受质量可以不改第一边缘、以不足二倍的行放大，转到**同一个完整 μ 的来源标签** z'。原 z 的所有门和历史留在 Π 第一坐标；z' 不被宣称具有原 z 的出生、LCA、FIRST 或 early 资格。Π 是收费耦合，不是新的实际路径。此接口对原有限 hard 目录同样有效：只需 R 合法和 (2)，未免费把 hard winner 扩为连续窗。它没有借用 hard 的所有尺度最大性，也因此还没有获得该最大性可能供应的空间增益。

正 Tonelli 给一个精确的 source-once 交换列：

\[
 \mathcal R_\tau
 =\int C_\dagger(z')d\mu(z'),\qquad
 C_\dagger(z')=\mathbb E_\theta\int_E
 \frac{T_\theta(x)}{D(x)}
 [h_{R_h(x)}-P_{\sigma(x),R_h(x)}]_+(x-z')dx.
 \tag{7}
\]

式 (7) 对定义的表示交通 R_τ 是精确等式，不添 J，不按输出重取 μ；全部历史仍通过 T_θ 和原选择器出现。只有精确原接口下可将左端改为 R_dagger；正支配接口下只能写 R_dagger≤R_τ=∫C_dagger dμ。

## 3. 正交换为什么没有立即支付空间列

原辅助噪声为 N_{σ,R}=[(1−σ)δ_0+σw_R]^⊗n，P_{σ,R}=h_R*N_{σ,R}。在闭 Q(x,R) 内，P_{σ,R}(x-z')=R^(−n)N_{σ,R}{z'+ζ∈Q(x,R)}，而 Q 外 h_R=0、P≥0。因此精确地

\[
 [h_R-P_{\sigma,R}]_+(x-z')
 =h_R(x-z')e_{\sigma,R}(x,z'),
 \tag{8}
\]

其中 e 是辅助退出概率。它使用同一原来源，不是 soft 首跳退出的 2W 参考测度。式 (8) 在边界上使用相同闭硬核规范；对本轮 L¹ 输入无来源边界原子问题。

在当前剩余 σ>1/n 上，已有原核解析下界 e>1/5；σ=2/n 可用 e>1/3。因而对 b=2a，删掉实际选择合同后的参考列满足

\[
 \frac13[1+n\log2]
 \le\int\sup_{a\le R\le2a}
 [h_R-P_{2/n,R}]_+(x-z')dx
 \le1+n\log2.
 \tag{9}
\]

该结论来自同一来源 hard-cone 列的解析公式，已经在 `auxiliary_exit_column_bridge_20261007.md` 核查。故将 (6) 替换为常数 2、再删 E/FIRST/选择器/余项门做 Schur，只能得到 O(n)W，不能取得 √n 费用。

这里严格区分两件事：(9) 否定“交换以后正缺额本身有 √n 列”的无门路线；它**不反驳** (7) 的实际加权列，因为后者保留 T_θ/D、hardband 和全部联合合同。用达到 (9) 的单来源或任意逐点尺度选择来冒充 R_dagger 输入，会丢失上述门。现有 e>1/5 已经逐来源处理了 selected gate 退出概率的下界，所以本稿不再把那一步重列为缺口。

真正仍欠的估计是同一 μ 的

\[
 \int C_\dagger(z')d\mu(z')
 \le A_nW+\varepsilon\lambda|E|,
 \quad A_n=O(\sqrt n\,n^{o(1)}),
 \tag{10}
\]

并覆盖全部当前 R_dagger。正耦合只给 (7)，没有给 (10)。将接收点依赖的密度 (T_θ/D)[h_R−P]_+ 当作一份固定新输入，会正好省略待付列；不能这样处理。

## 4. strict CP 的捕获—未捕获质量交换及其失去的支撑

另一尝试保持原 θ、唯一实际父组 P，并先看 R_h≤L_s 序。记 h_P=ν_h(P) 为原固定完整出生质量，m_P(x)=ν_h(P∩Q(x,R_h))。strict CP 与 hard 行给

\[
 \alpha_P(x)=m_P(x)/h_P<C_hv_j/\sqrt n,
 \tag{11}
\]

在右侧<1 的一般目标维度成立。将当前实际交通依其唯一 P 正分配为 τ_{x,θ,P}，T_P=τ_P(R^n)。软后验总量≤1 给 T_P≤m_P/R_h^n，且 Σ_P T_P=T_θ，不能让每个 P 再获得完整 T_θ。

当 m_P<h_P 时，令

\[
 d\chi_{x,P}(z')=\frac{1_{P\setminus Q(x,R_h)}(z')}{h_P-m_P(x)}d\nu_h(z'),
 \quad d\Pi_P=d\tau_P\otimes d\chi_{x,P}.
 \tag{12}
\]

第一边缘仍是原 τ_P；第二边缘的 ν_h 密度满足

\[
 \frac{T_P(x)}{h_P-m_P(x)}1_{P\setminus Q(x,R_h)}(z')
 \le\frac{\alpha_P(x)}{1-\alpha_P(x)}R_h^{-n}
 1_{P\setminus Q(x,R_h)}(z').
 \tag{13}
\]

这确实产生了一个行上的 O(1/√n) 因子；但它把收费来源换到了**原 hard 盒外**。右边不是 h_R(x−z')，无法调用来源中心的 N_h hard column。对固定 P，即使 T_P 的原输出必在 P+Q_b，粗空间列只是

\[
 \int_{P+Q_b}a^{-n}1_{z'\notin Q(x,R_h(x))}dx
 \le |P+Q_b|/a^n,
 \tag{14}
\]

父组边长大时没有维数多项式上界。完整 fixed/adaptive core 工具已经支付有合格证书的来源；R_dagger 恰保留证书不适用或原 Rh 大于保存认证半径的部分。不能将那个工具的好支费用重复用于 (14) 的任意剩余 P。

此外 (12) 的参考来源跨父层重用：唯一 LCA 只保证原来源对 (y,z) 的分配唯一，不保证交换后的 z' 在各 P 的收费列中只出现一次。若逐 P 只用 ν_h(P)，Σ_P h_P 最多仍有 J；原 lowcoin 的 source-fixed 容量才有相应已证多层预算，不能把 (11) 的输出捕获比例当成新的 prior capacity。R_h≥L_s 的另一序只给原 soft 加权响应比例小，根本不是一个硬捕获集合；不可把它代入 (12) 冒造未捕获质量。

这条交换失败的精确位置是：小捕获比例在转移第一边缘后，换成了没有来源中心空间支撑的第二边缘，并可能跨层复制参考质量。不是 (11) 算术错误，也不是 complete h 可以替成 captured m。

## 5. 判断与保存范围

本轮新增可核结果仅为 (5)–(7) 的显式行耦合和 (12)–(13) 的捕获补质量耦合；没有新的阶数假设、已付分层或全覆盖费用，故不登记新下界构造数值，也不重跑已有三轮 shell、packet、fixed-delay 或 frozen 工具。相应 (9) 是已有原核的全 n 解析阻碍，不靠数值拟合；它只拒绝删门的收费办法。

若继续沿共同实际赢家交换，需要对 (7) 保留的 weighted column 给真正的 spatial certificate，或使 (12) 的参考来源选择具有同一 μ 的 source-once 支撑/容量证书并支付不合格补集。现有行最大性、τ≤u、fullfuture q 和 strict CP 小比例没有提供这些证书。原 hard 连续最大性可以给额外尺度关系，但不能改变 (13) 的盒外支撑；原有限目录还须单独保留合法尺度范围。

没有修改主账，R_dagger 尚未支付，一般 geom 的 √n n^{o(1)} 目标未完成。
