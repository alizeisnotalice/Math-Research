# 唯一 LCA 容量失败的可付覆盖分层与 sibling 占用障碍

2026-10-07。只写本稿及同前缀守卫，不改旧原件。本轮给一项作用于原 actual 失败子交通的真实空间支付：覆盖概率 `ℓ≥δ` 的失败部分可由已有 source-fixed 收据预算放大支付。剩余 `ℓ<δ` 严格保留低覆盖及两条完整质量条件；没有从这些条件取得整个余项根号费用。

## 1. 查重与必须恢复的原 actual 失败权

已读原 tex 2591–2758 的 SC-F/R、HR-COST/PARTIAL；`high_marked_cover/next_spatial_capacity.md` SC1–5；`high_marked_receipts/route.md` HR3–5、`carleson_audit.md` CA1–8；`GENERAL_CAPACITY_TREE_LEDGER.md` §§2–6、`GENERAL_CURRENT_REDUCTION.md` §32；现有 `gated_occupancy_bridge.md`、`actual_parent_occupancy_route_20261007.md` 与前轮 actual 门审计。CA 已用 fixed-hard-source 的互不相交 sibling rings 去掉伪 J，得到旧 `N_hW`；source-cone 只改坐标。小空间祖先、轻质量祖先和固定密度父组也已有各自范围，不重新登记为本轮新工具。

本轮作用于**当前未付 actual 子交通**，保留完整 μ、μhi、q、σ/Ls/Rh、FIRST/futurecap、hardband、GOOD/score/owntrace、出生、原森林、共同重捕获、严格 CP/GP、实际 early firstjump exit 与 full continuation，以及所有已付删除补集／新 lowS 合同。

固定原几何平移 θ 的 source forest。令 Dθ 为 sameTop 且 y,z 处于不同最细格、从而存在合法 distinct-child 唯一 LCA 的 pair 域；原 actual far 子交通位于此域。本文公式中的 `1_sameTop` 实际指 `1_Dθ`，不包含同一最细格而无 distinct-child LCA 的 pair。父组 P 与孩子 Cy,Cz 由原 `(y,z)` 的唯一 LCA 决定。域外置 G=ℓ=w=0，关系 w=1−ℓ 只在 Dθ 上使用。写原完整 source masses

\[
m(C)=\mu_{hi}(C),\qquad h(C)=\mu_{hi}(C\cap B),\qquad B=\{t_C\le b\},
\]
\[
p_C=\min\{1,m(C)/(\sqrt n\,h(P\setminus C))\},\quad
r_C=\min\{1,h(C)/(n\,m(P\setminus C))\}.
\tag{1}
\]

零分母取一。原 forward/reverse **辅助币**对各 source child/direction 独立、在所有输出上固定；它们不是 FIRST/history/source posterior 的坐标。几何 θ 的概率平均与这两个币的积分分别处理。

取原 actual 基础交通的可测版本

\[
d\Pi_0=q^{-1}1_E(x)P_{\sigma,L_s}(x-y)h_{R_h}(x-z)
\,dx\,d\mu_{hi}(y)\,d(1_B\mu_{hi})(z).
\tag{2}
\]

纯角向资格也可留在基础核；下文正上界可去掉。合并唯一分配后，原失败子交通有表示

\[
\boxed{dR_{cur}=1_{sameTop}\,G(x,y,z,\theta)
(1-p_{C_y})(1-r_{C_z})\,d\Pi_0\,d\Pr(\theta),
\quad0\le G\le1.}
\tag{3}
\]

`G` 保存其余全部门。实际已积分 early 子核比值 `Qearly/P≤1` 可以位于 G，但不重复乘该比值。如果较晚门依赖原 HR 币／历史，用在“双币失败”条件下的门接受率定义 G；只有双币失败的先验质量仍为 `(1−p)(1−r)`，不假定门与币独立。w=0 的位置取 G=0。

重要范围：最新 prompt 146–153 只提供某 actual 交通或正支配的 `A≤1`。若先做正扩大而删除了原 failure 权，只知 `A≤1` 的扩大对象**不能**用本轮证明。必须在原实际子交通中恢复 (3)，或证明其 `A≤w=(1−p)(1−r)`。当前本文按原 HR 失败历史及后续 actual 删除的子项使用，不把一个任意 `A≤1` 重新命名成失败交通。

## 2. 原 joint union 的真实空间费

定义 source-fixed pair quantities

\[
\ell_\theta(y,z)=p_{C_y}+r_{C_z}-p_{C_y}r_{C_z},\quad
w_\theta(y,z)=1-\ell_\theta(y,z).
\tag{4}
\]

条件于 θ，ℓ 是两条原收据联合包含此 pair 的准确概率；不同孩子／唯一 LCA 保证其它层收据不能包含这对。其空间预算来自完整正收据，而不是单行概率。

原 SC-F/R 与 HR-COST 已给

\[
\mathbb E_\theta\int1_{sameTop}\ell_\theta\,d\Pi_0
\le J A_nW_{hi},
\quad A_n=N_h/\sqrt n+C_hJ_s/n,
\tag{5}
\]
\[
N_h=1+n\log(b/a),\quad C_h=16/3,
\quad J_s=D_{full}(1+2n\log4),
\quad D_{full}\le1+(\pi/2)\sqrt n.
\]

证明沿用旧 source-only receipts：每层 `Σ_C p_C h(P\C)≤W_hi/√n`、`Σ_C r_C m(P\C)≤W_hi/n`，先用原 soft/hard row，后用真实 Lebesgue column；`ℓ≤p+r` 扩为收据和。每层完整来源计一次，既有 J=O(log(n+2)) 层显式收费。没有按父组数量、后验历史、接收点数或候选 pair 再复制 W。

G≤1、实际 early 子核≤P、后来 source/output 删除只缩小 (5) 左侧。不能用 reference 2W 来源继承原 q cap；(5) 始终在原 kernel 与原来源下证明。

## 3. 当前失败交通的可付 high-union 分支

固定任意预先规定 `0<δ≤1`，只在 (3) 的当前实际剩余内分割 `ℓ≥δ` 与 `ℓ<δ`。对 high 支，严格点态

\[
1_{\ell\ge\delta}w
\le\frac{1-\delta}{\delta}\,\ell.
\tag{6}
\]

因此，保留原门直至正估计，得到真实 source-W 空间费

\[
\boxed{R_{cur,\ell\ge\delta}
\le\frac{1-\delta}{\delta}J A_nW_{hi}
\le\frac{1-\delta}{\delta}J A_nW.}
\tag{7}
\]

等号 `ℓ=δ` 归 paid；δ=1 的 high 支 w=0，(7) 正好为零。取

\[
k_n=\lceil\log_2(n+2)\rceil,\qquad\delta_n=k_n^{-4},\quad n\ge512,
\]

这是预先固定且可精确复现的有理 coin 阈值，与 `log^-4(n+2)` 同阶；它不改变前轮 lowS 的原 `ε_n=log^-4(n+2)`。则 (7) 为 `O(√nlog^5(n+2))W`，可置于现有 `O(√nlog^6)` 主预算。它是旧空间 receipt 的一个新明确 paid subbranch／合同强化，未发明新核列理论或证明原余项自动 high。

原 union 已付 `ℓdΠ0` 与当前 high-failure `G w1_{ℓ≥δ}dΠ0` 互补，不能声称原 union 已支付这一整份失败项。按本轮增量处理时，保留旧费用一次，再为当前 high 子支付 (7)。若重组完整原 HR union 项，可用

\[
\ell+Gw1_{\ell\ge\delta}\le\ell/\delta
\]

将同一个旧 `J A_n` receipt budget 合并成 `J A_n/δ`，不再另加一次旧与新增费用。此合并只适用于原 HR 部分的已核账；不能因数值相同就重排未知上游三个尺度项。

### 主账系数与 sameTop 的范围

(7) 对原 `sameTop` 失败交通直接成立，**没有额外乘二**。若从最早的整个 R_HMP 开始，旧 HR 的 sameTop 几何覆盖只给≥1/2，需沿旧式

\[
R_{HMP}\le2J A_nW_{hi}+2R_{bal}
\]

再代入高支支付；这里的二仅是原 sameTop cover，一次使用。旧 MM/HR 总账外层64/63因而给128/63，不把这个历史系数搬回当前 prompt。

当前最新条件账是

\[
\lambda|E|\le2B_{paid}+(4096/49)R_{heavy}.
\]

若当前 (3) 是 Rheavy 的剩余子项，本轮 fee 在这一账中乘 **4096/49 一次**。若同时保留目前已核的 `49/8192 λ|E|` 吸收，即

\[
R_{heavy}\le B_{other}+B_{coin}+(49/8192)\lambda|E|+R_{lowcoin},
\]

则移项后

\[
\lambda|E|\le4B_{paid}
 +(8192/49)(B_{other}+B_{coin}+R_{lowcoin}).
\tag{8}
\]

`Bcoin=((1−δ)/δ)J A_nW` 已含 Ch 于 A_n 中；不额外乘 Ch、sameTop 二或旧128/63。最新三个尺度整账仍沿 prompt 接口前提，本稿不以 (8) 认证其缺失上游原件。

## 4. low-union 的完整质量合同及 sibling 链检验

在合法域 Dθ 内，补集 `ℓ<δ` 等价于 `w>1−δ`；由于 `ℓ≥max(p,r)`，必有

\[
p_{C_y}<\delta,\quad r_{C_z}<\delta,
\]
\[
\boxed{m(C_y)<\delta\sqrt n\,h(P\setminus C_y),
\qquad h(C_z)<\delta n\,m(P\setminus C_z).}
\tag{9}
\]

实际正交通的 endpoints 保证 relevant m(Cy),h(Cz)>0，分母也为正；零分母时对应币为一而不进 low 支。严格 (9) 要求原**完整来源质量**，不换为每行当前 captured subsource，也不从它推出两个孩子质量可比。`p<δ,r<δ` 只是必要条件，不能反过来当成 `ℓ<δ`。

唯一 LCA 能否再给可求和额外增益？对固定 hard z，原 sibling rings

\[
S_j(z)=P_{j+1}(z)\setminus P_j(z)
\]

互不相交。定义带完整 actual 门与本轮 low selector 的原 soft 接受响应 a_j(x,z)≥0，则

\[
\sum_j a_j(x,z)\le q,
\qquad
R_{lowcoin}\le\int_B\mathbb E_\theta\int_E
h_{R_h}(x-z)\sum_j a_j(x,z)/q\,dx\,d\mu_{hi}(z)
\le N_hW_{hi}.
\tag{10}
\]

这是 CA1–3 已有 source-once O(n)W，没有新 J，但也没有新根号费。`ℓ<δ` 的 failure 权趋近一；它不是乘回 (10) 的小 δ 系数。

一个明确 source-tree arithmetic 障碍：每层有 b=512 个等质量、全 born 的 occupied children（其余孩子质量零）。沿固定 hard z 的 child 向上，parent 总质量为 b 倍 child，每层 relevant coins 都恰

\[
p=\frac1{\sqrt n(b-1)},\qquad r=\frac1{n(b-1)},\qquad
\ell=p+r-pr,\quad w=(1-p)(1-r).
\tag{11}
\]

取任意有限深度 D，最内 child 质量 b^-D；第 j 层 child 质量 b^(j−1−D)，新增 b−1 个同质量 sibling，故 rings 总质量恰 `1−b^-D`，完整来源 W=1。n≥512 的真实 dyadic n-cube 可有至少512个 occupied children；但此处只使用有限质量树，不赋予其实际输出资格。三轮 n=1024,4096,16384 时 (11) 严格满足 `ℓ<k_n^-4`，且 w 沿 ancestor index 恒定，未出现额外可求和衰减。即使链加长，`w(1−b^-D)` 仍接近一。它只检查低币与树关系本身没有给出所需小系数，**不认证 FIRST、far、CP/GP、early 或 actual 空间反例**。

更一般在未 clipped 的 distinct-child 情形，

\[
p_{C_y}r_{C_z}
=\frac{m(C_y)h(C_z)}{n\sqrt n\,h(P\setminus C_y)m(P\setminus C_z)}
\le\frac1{n\sqrt n}.
\tag{12}
\]

binary／全 born 两孩子时可等号；多孩子可以使乘积更小。(12) 是 source-only 的上界，不能换成 ℓ 的下界；任意多小孩子允许两个 coins同时很低，其个数不能免费忽略。也不能从不同 rings 的全局质量增长，推 original posterior a_j 的相同比例增长；kernel、两个尺度、输出及 fullgate 会重权，原行只给总和≤q。

对 ℓ作 dyadic inverse-probability 层分割是合法的，但每层 (7) 的 `δ^-1` 会增长；(5) 只付 `∫ℓdΠ0`，不控制 `∫ℓ^-1ℓdΠ0` 的全尾。反复增大币或每层重领收据预算，需支付相应增加的 repeated source masses；CA、原容量证明和 source-cone 身份均没有这种额外可求和结论。

本轮最小剩余是原完整实际合同再交 `ℓ<δ_n` 的 source-pair 空间积分。需要从 full future／far／strict CP-GP／实际 history 取得新的空间增益；source-tree 代数与原两个收据费用尚未提供它。

## 5. 三轮精确守卫与结论范围

主命题 (7) 来自旧已核空间收据和严格点态比较；新守卫只核完整有限 source-tree 的质量容量、joint union/failure、原失败权恢复、同一份来源及 threshold 端点，不模拟 actual FIRST 或重跑旧数据。δn=k_n^-4 使用整数 bit_length 与 Fraction 精确计算；不需要认证任何浮点或自然对数端点。

脚本 `capacity_pair_occupancy_exact_guard_20261007.py` 已终态执行，结果在 `capacity_pair_occupancy_exact_guard_results_20261007.json`：n=1024、4096、16384，√n=32、64、128，k_n=11、13、15；共 **55,674 项 exact Fraction 检查全部 PASS**，无随机种子、无 live handle。

每轮三种完整质量组件：32叶 binary/full-born 树、16叶 skew/partial-born 四叉树、精确压缩的512-child 星形及四层链。前两种树逐 ordered distinct pair 查唯一 LCA，使用相关于 pair 的任意接受率 G，并逐完整层核 (5) 使用的 source capacity 不等式；第三种按两类 parity orbit 精确求和，未用采样替代大树。每种组件分别测试主有理阈值、实际 ℓ 等号端点和 δ=1，总计27组 summed fee／互补 split 检查。512-child 链在三轮均 strict-low，完整 rings 总质量 `1−512^-4`、failure-weighted ring 总质量>99/100；这是有限质量树中没有额外祖先衰减的范围检查。

这些记录没有对真实 Lebesgue spatial column、FIRST、birth clocks、CP/GP、nonconcentration、far、early exit 或 full history 作数值认证。尤其压缩链不是 actual 余项反例，不宣称能满足原剩余所有门。空间费 (7) 的证明使用原 SC/HR 已证收据预算和点态 (6)，而非这些有限代数数据。旧／新账系数的独立有理核验另见 root 保存的 `capacity_ledger_scalar_review_20261007.json`；它也只核算术。

本稿的增量是把同一原旧 receipt 的覆盖预算放大，支付当前未付 failure 内的 `ℓ≥k_n^-4` 子支；与原 union/上游 paid 分支依次互补记账。最终 lowcoin∩lowS∩其余 actual 门的空间费仍未闭合，本轮没有证明完整一般 geom 的 `O(√n n^o(1))` 上界。
