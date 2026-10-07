# 实际重父组的捕获分数与空间占用：查重及最小缺口

2026-10-07。结论：未从共同 strict CP/GP、ℓ¹-far 和原 FIRST/winner 得到新的剩余空间费用。本稿保留两个尺度顺序下的正确行约束，列出旧共享覆盖工具的一个可选反链应用；不将它接入主账，不把既有 CA/source-cone 表示重报为新进展。

## 1. 实际门排除了上一压力模型

`joint_low_future_contract_probe_20261007.md` 的原核/FIRST/futurecap/非集中合同仍成立，但不是当前实际余项样本。均匀二元顶点输入的 ℓ¹ 直径为 nd=ka；原 tex 2354–2385 已支付所有距离≤H_mark 的来源对，三轮此直径远小于 H_mark。其 W=1 又使 h(P)≤1，v_j≤1，而 q min(L_s,R_h)ⁿ√n≥√n/3>1，所以 strict CP 也不能成立。实际余项为零，无需再运行 φ。

原 tex 2500–2610 的两盒 CW true-FIRST 构型不能未经实际门核验就当新残余样本。该处展示的大 far-soft 质量与小 near-hard 质量的跨盒交通，可由 SC-F 的固定来源容量引理支付；反向跨盒可检查 SC-R。这里不由跨盒支付推断整个 CW 输入在全部输出处的实际余项为零，同盒及其它输出仍需原资格核验。仅重新构造低 Z、futurecap、far 和 CP/GP 数值，若省略原容量币与平衡分数，仍未触及目标。

## 2. 本稿保留的固定来源合同

固定完整原输入、共同 q、原 FIRST σ,L_s、hardwinner R_h、出生集合 B、共享源森林种子及唯一 LCA 父组 P。记

\[
 M(P)=\mu_{\rm hi}(P),\quad h(P)=\mu_{\rm hi}(P\cap B),
 \quad S_P(x)=\int_P P_{\sigma,L_s}(x-y)d\mu_{\rm hi}(y),
\]
\[
 H_P(x)=R_h^{-n}\mu_{\rm hi}(P\cap B\cap Q(x,R_h)).
\]

始终 M≥h。h 是完整出生合格质量，不是当前捕获质量。原 soft 行给 S_P≤q；原 hard 行给 H_P≤C_hq。保留 ℓ¹-far>H_mark、原 SC-F/R 或 HR 容量币、平衡权重 (1−p_Cy)(1−r_Cz)、出生、原全部历史/first-exit/full-continuation，以及最新低未来包络门。

第 j 个预定 soft 带的 v_j=(1−β₀s_j)^{n/2}∈(0,1]，β₀=3/(4e)。目前未付父组同时满足

\[
 v_jh>q\min(L_s,R_h)^n\sqrt n,
 \qquad h>q(L_sR_h)^{n/2}. \tag{1}
\]

这些条件及其参数都沿原 CP/GP 使用，不重算 FIRST 或来源。

## 3. 双序的正确小响应约束

若 R_h≤L_s，实际 hard 捕获比例

\[
 \gamma_h(x,P)=\mu_{\rm hi}(P\cap B\cap Q(x,R_h))/h
 \le C_hqR_h^n/h
 < C_h\min\{v_j/\sqrt n,(R_h/L_s)^{n/2}\}. \tag{2}
\]

若 L_s≤R_h，则应使用完整 soft 核的归一化响应

\[
 \gamma_s(x,P)=L_s^nS_P(x)/M
 \le qL_s^n/M
 <\frac hM\min\{v_j/\sqrt n,(L_s/R_h)^{n/2}\}. \tag{3}
\]

(3) 是带原 soft 权重的响应比例，不是 soft source 被某 hard 盒捕获的质量。保留 h/M，不能换成 1 后宣称取得更强来源容量。等尺度时两式都可使用。

它们只说明较小尺度端的父组平均响应很小；并不说明每份原来源的空间列小。`large_parent_capture_route.md` 的 LP1–LP7 已核 hard 序的取消：将 (2) 乘回 h/R_hⁿ只恢复 H_P≤C_hq；分割当前 capture 为新硬子测度会依赖 x。soft 序把 (3) 乘回 M/L_sⁿ只恢复 S_P≤q。把同一完整 soft P 配给多个固定 hard patch，或按每个输出重新选择 patch，会重复领取来源质量；已有 SC-F/R 不允许这样使用容量。

## 4. 已查重的来源一次占用接口

已读原 tex 2608–3015、4766–4869，及 `high_marked_receipts/large_parent_capture_route.md`、`carleson_audit.md`、CP threshold 和 GP ledger。CA1–CA6 已利用 hard 来源 z 的互不相交 sibling rings，将原软行总和压到 q，并用一次原 hard 列得到 N_hW；没有 J 的额外损失，但仍是 O(n)W。现有 source-cone Γ_v* 表示也已精确记录，不能把换坐标本身当新费用。

最小尚缺量仍是原实际占用：对原 hard 来源 z，保留真正来源对及历史门后的接受权 α_θ(x,z)∈[0,1]，其中来源级别/带/唯一 LCA 的和已经计入，研究

\[
 \mathcal R_{\rm rem}
 =\int_B\mathbb E_\theta\int_E
 h_{R_h(x)}(x-z)\alpha_\theta(x,z)\,dx\,d\mu_{\rm hi}(z). \tag{4}
\]

等式按原门条件平均的 RN 接受权理解；若仅有子核正支配，则用对应上界版本。原 y 标签、两个不同尺度及所有门在 α 内，不能用 reference 扩大后的总质量反推原边际。来源只计一次。需要证明 (4)≤√n·n^{o(1)}W，或给其可支付互补分解。

(2)–(3) 是先在同一个 x 对完整父组来源积分后的约束；(4) 需要固定 z 后对不同 x 的真实 Lebesgue 占用。没有证明允许交换这两个约束并产生 v_j/√n 增益。将 α 删为 1 只回到已知 N_hW；将现有 Γ 的平方或 radial 高占用票据断言为小，又是在增加未证接口。本稿不另命名这种未证量为引理，不创建 abstract extremizer 或特殊核实验来代替它。

## 5. 旧共享覆盖的可选直接应用：输入密度父组

这项只应用旧 tex 4100–4135 的 q|U| 支付及 BA 最大祖先反链思想；不是新核/空间列理论，也未登记为本轮新增主账。固定 T>0 和森林种子，在输出之前标记所有满足

\[
 q|P+Q(0,b)|\le T M(P) \tag{5}
\]

的父组。P 为边长 ℓ_P 的轴盒时，其闭包 Minkowski 体积恰 (ℓ_P+b)ⁿ；半开边界不改变该 Lebesgue 体积。选所有最粗的标记祖先 A。森林有限层，所以每个标记 P 有唯一最大标记祖先；这些 A 构成可数不交源反链，Σ_A M(A)≤W_hi。标记条件不要求跨层单调。

任何 LCA P 被标记的实际源对，其 hard z∈P⊂A，R_h≤b，因而原输出 x∈A+Q(0,b)。固定种子的 U_θ=∪_A(A+Q(0,b)) 支持这份子交通。先将唯一 LCA、唯一输出带及原 pair 门合并，再对整个子交通使用原行密度≤q⁻¹·q·C_hq=C_hq；不是逐 parent 分别付一份行上限。于是

\[
 \mathcal R_{\rm marked,\theta}
 \le C_hq|U_\theta|
 \le C_hT\sum_A M(A)\le C_hT W_{\rm hi}. \tag{6}
\]

再平均原种子无新倍数。可进一步在原带种子交通中删除整个 U_θ，费用仍按 (6)；不得把所有种子的 U_θ 并集当作一个免费输出集合。取 T=√n 给合法的 √nW 可选子支；来源反链避免逐层的 J 或逐带 K_s 重复收费。

这个明确合同不要求原 Q 捕获整个 A，也未把 h 改成 M。CP/GP 仍使用原完整 h；(5) 是另一个完整 M 的来源固定判据。它不能支付所有父组：补集为 q(ℓ_P+b)ⁿ>√n M(P)，若删除整个 U_θ，则原 LCA 的所有标记祖先也都被排除。该密度缺口与 strict CP/GP 可以共存，不推出实际输出稀少；ℓ¹-far 本来已排除了靠近小父组的旧支付情形。没有证明剩余实际交通全部落在 (5) 的覆盖中。

## 6. 本轮终态

没有新一般空间列证明，没有重跑 φ 或任何冻结数据，没有拟合阶数。上述可选 corollary 的解析范围与旧共享覆盖工具相同，不为旧工具再造三轮 toy 压力，也不将其主账费用再领一次。真正核压力的本轮进展是实际 near/CP 门排除了上一已完成模型，明确了后续测试必须先过哪些原门。

若继续研究，必须在 (4) 内利用完整 FIRST/futurecap、ℓ¹-far、strict CP/GP、容量平衡和实际历史的联合几何；仅增加父组行质量条件或删除占用门都不足。最小空间接口仍未证，一般 √n 目标未完成。
