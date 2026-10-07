# 完整 future 与 source capacity 的未归一化 Bernstein-tail 空间桥

2026-10-07。新的可付子支：在原 actual lowS∩lowcoin 交通上保留所有原标签与门，以原核的正 Bernstein 分解作诊断标签；按 source-fixed forward 容量 p_C 选择尾阈值。每个孩子的未归一化 tail 响应≤p_Cq，随后固定来源容量和真实 hard column 支付 `J N_h W_hi/√n`。这不是重复 coin amplification，不证明完整剩余空间费。

## 1. 查重与接口范围

已读原 `概率接口阶段证明.tex` 730–865 的 FIRST-MGF、FACE、MIXED、RANK；`signed_band_probability/entropy_characteristic_route.md` 及其独审；`high_marked_cover/fractional_actual_pair_budget.md`；旧 SC-F/R、CA 与本轮之前的 capacity/actual-parent 稿。

旧 FIRST 已给每行的 count tail，并说明条件来源 C 会出现 posterior 质量 `π_C^-1`；旧 entropy characteristic 保留群导数及边界，核体积权恰取消正 n；fractional LP 的实际 demand/source capacity 仍未付。本稿不重证这些通用结果，不把群导数设零或重报 source-cone/CA 的 `N_hW`。

本稿使用原完整 μ、q、停止软度 σ(x)、Ls(x)、硬赢家 Rh(x)、完整 future cap

\[
P_{s,L_s(x)}*\mu(x)\le q\qquad(s\in[\sigma(x),1]),
\tag{1}
\]

以及 source-only 森林种子 θ、合法 distinct-child 唯一 LCA、原 born hard 子源 `ν_h=1_B μ_hi`。所有 FIRST/GOOD/score/owntrace/birth/重捕获/strict CP-GP/far/nonconcentration/短壳/实际 early firstjump/full continuation、新 lowS 与 lowcoin 门保持。在本轮证明起点仍为原 failure 权 `G(1−p_Cy)(1−r_Cz)`，不能对已丢失 failure 的任意扩大 A≤1 追认原失败合同。

## 2. 保留实际历史的正诊断标签

对原核写

\[
P_{\sigma,L}(x-y)=\sum_{A\subset[n]}P^A_{\sigma,L}(x-y),
\quad P^A_{\sigma,L}=L^{-n}\sigma^{|A|}(1-\sigma)^{n-|A|}
\prod_{i\in A}\phi((x_i-y_i)/L)\prod_{i\notin A}h((x_i-y_i)/L).
\tag{2}
\]

在原 actual 测度每个 `(x,y,z,θ,全部真实历史)` 标签上，额外按正权 `P^A/P` 分配，不改变其原质量；`P=0` 的原基础交通为零，可任意定义诊断标签。所有原门连同实际 firstjump/continuation 保留在每份分配中，正权和恰一。

记诊断 `K=|A|`。它是原核 Bernstein mixture 的标签，**不是原 early firstjump 坐标 i、continuation 活跃数、jump 次数或已揭示路径信息**。尤其下面的 K 上限不约束这些真实历史。若先把原历史条件平均成 G≤1，使用同一个 G 乘每个 `P^A` 也恢复原交通；不声称原历史与任何自然路径 mask 独立。

## 3. 源容量相关 tail 的准确响应及端点

固定 θ、原父组 P 和其中 soft 孩子 C，记

\[
m_C=\mu_{hi}(C),\quad H_{opp}=\nu_h(P\setminus C),\quad
p_C=\min\{1,m_C/(\sqrt n H_{opp})\}.
\tag{3}
\]

原零分母约定 p_C=1；正 traffic soft endpoint 使 m_C>0，故 relevant p_C>0。若 m_C=0，则对应原来源交通为零；可把其 tail 置空，不取 log0。

预先以完整来源定义

\[
u_C=\log(1/p_C),\quad v(x)=n\sigma(x),\quad
\kappa_C(x)=v+\sqrt{2vu_C}+2u_C/3.
\tag{4}
\]

tail 为 **K≥κ_C(x)**，整数 K 因而使用 `ceil κ`，等号进入 paid。若 κ>n 则为空。p_C=1 时 u_C=0、κ=nσ，界仍成立；σ=0 时 K=0，p_C<1 的 tail 为空，p_C=1 的 tail 全收；σ=1 时 K=n，p_C<1 的 tail 为空，p_C=1 的 tail 全收。实际余项本来 σ>1/n，端点约定只是正核工具的完整范围。

设未归一化群 tail 响应

\[
T_C(x)=\int_C\sum_{|A|\ge\kappa_C(x)}P^A_{\sigma(x),L_s(x)}(x-y)\,d\mu_{hi}(y).
\tag{5}
\]

对任意 t≥0，原 Bernstein tilt 恒等式和 positivity 给

\[
T_C(x)\le e^{-t\kappa_C}B_x(t)^n
\bigl(P_{s_t,L_s(x)}*(\mu_{hi}|C)\bigr)(x),
\quad B_x(t)=1-\sigma+\sigma e^t,\quad s_t=\sigma e^t/B_x(t).
\tag{6}
\]

这里只测试原共同物理 Ls；s_t≥σ，输入子测度≤原完整 μ，因此 (1) 对每个 C 的最后一因子给 q。没有重新定义 C 的 FIRST 或归一化成条件来源。

原通用 Bernstein 估计给

\[
\inf_{t\ge0}e^{-t\kappa_C}B_x(t)^n\le e^{-u_C}=p_C.
\tag{7}
\]

为自含核阈值：`log B^n−vt≤v(e^t−1−t)≤vt²/[2(1−t/3)]`，0≤t<3。若 v,u_C>0，取 r=κ−v、t=r/(v+r/3)；所得 exponent≤`−r²/[2(v+r/3)]≤−u_C`，因为 `r²−2u_C(v+r/3)=(2u_C/3)√(2vu_C)≥0`。u_C=0 取 t=0；v=0 用上述确定 count 端点，避免 t=3。κ>n 的空尾无需使用有穷 tilt。σ=1 的确定 count 也已单列。

因此关键**未归一化**预算严格为

\[
\boxed{T_C(x)\le p_Cq.}
\tag{8}
\]

没有 `π_C^-1`，也没有先按 C 条件化再遗忘条件质量。不同 C 的 tilt 参数可以不同，因为本稿下一步逐 C 以固定来源容量结算，不把不同参数的响应拼成同一份 row q。

## 4. 真正空间付款：固定 child 容量与一次 hard column

只在当前未付 actual lowS∩lowcoin 交通中付 (5) 的 tail。保留全部历史后去掉 G、failure≤1 及其余门作正估计；按原唯一 LCA 层及 soft child C 分组，hard endpoint 必在原 `P\C`，得到

\[
R_{tail,\theta}
\le\sum_{\mathrm{levels},C}\int_E\frac{T_C(x)}q
\bigl(h_{R_h(x)}*(\nu_h|P\setminus C)\bigr)(x)\,dx.
\tag{9}
\]

每个原 pair 只出现在其唯一合法 LCA；同最细格/无合法 distinct-child LCA 不进入本交通。原森林每层的不交 children 给 Σ_Cm_C≤W_hi。应用 (8) 及真实 hard radial column

\[
\int\sup_{a\le R\le b}h_R(x-z)\,dx=N_h=1+n\log(b/a)
\tag{10}
\]

后，

\[
R_{tail,\theta}\le N_h\sum_{\mathrm{levels},C}p_CH_{opp}
\le\boxed{\frac{J N_h}{\sqrt n}W_{hi}}.
\tag{11}
\]

最后一步是原 source-fixed forward capacity `p_CHopp≤m_C/√n`，包括 Hopp=0、clipped p=1 情况。平均原 θ 无新增费用，**没有 Ch**：在 soft tail≤p_Cq 后按 hard 的真实 column 收费，不再使用 hard row Chq。来源每层用一次，原 J=O(log(n+2)) 明收，费用 `O(√nlog(n+2))W`。

这是旧 forward capacity 与全 future 倾斜的新明确组合子支；不是另立更强核列定理。它与过去 paid 交通的交集先删除，只在当前补集收费一次。它不意味着原 union 已付了这些 failure tail；当前新增 fee (11) 与 coin/highS 等先前子支依次互补登记。当前保留49/8192吸收的账，回代 fee 乘8192/49一次；若仍在未移项的Rheavy账则乘4096/49一次，不用旧128/63，也不额外乘sameTop二。

## 5. 两个不可用交换及剩余量

对固定 θ,z，合法 sibling rings S_j(z) 不交；共同未来参数给 `Σ_j P_{s,L}*(μ_hi|S_j)≤q`。但若每环各选 `(s_j,L_j)`，该和不再是原同一参数卷积，完整 future 只提供逐项≤q；source disjointness 不允许把不同核选择的响应合回一份 q。用 H 全系数包络可以正放大，却要保留它真实 O(√n) shape mass 与物理 dilation fee。本稿 (9) 的逐 C 容量结算避免该交换。

若另加 source capacity 价 a_j，例如 hOpp/m_C，准确倾斜得到的是

\[
\frac{B(t)^n}q P_{s_t,L}*\Big(\sum_j a_j\mu_{hi}|S_j\Big)(x),
\tag{12}
\]

输入变成了加权测度；原 future cap 只对完整 μ 与≤μ的子测度可用。只能据 sup a_j 支配，不能免费沿用 q、不能交换 source-fixed prior capacity 与 posterior 响应比例。LP5 与 entropy characteristic 正是仍保留这种空间重用/群导数，未供应遗漏的容量界。

付款 (11) 后，严格剩余是同一原实际测度的正分配，再交所有原 lowS∩lowcoin 门及 `K<κ_C(x)`。这是诊断-label traffic 的 exact complement；不要把它改写成原 pair 没有尾质量、不要把 K 换成原 continuation 活跃数。随着 p_C→0，u_C增加，阈值可超过 n；此时 tail 为空而原 traffic 全在补集。由此尾付款不能覆盖所有极小容量来源对。

Sibling rings 不交仍只给旧 `N_hW` 的剩余基线；K 的诊断上限未给 fixed hard z 的 selected空间占用衰减。完整一般 geom 的 `O(√n n^{o(1)})` 目标尚未闭合。本稿不构造抽象/原核 toy 当 actual counterexample。

## 6. 三轮精确守卫与终态范围

只核新正分解/tilt/tail与完整来源容量结合的离散系数；不会数值重跑 φ、旧 FIRST 或 actual history，不认证任何连续 FIRST 输入样本。证明的空间付款由 (1)、(8)、原容量和真实 column (10) 给出，不由有限系数模型拟合阶数。

`actual_future_capacity_bridge_exact_guard_20261007.py` 已 exit0；同前缀 `actual_future_capacity_bridge_exact_guard_results_20261007.json` 保存三轮 n=1024、4096、16384，√n=32、64、128。18组正系数模型及15组端点，共 **24,171 项 exact integer/Fraction 检查全部 PASS_EXACT**，没有随机数或 live handle。

每轮 σ取 `2/n、1/(4√n)、1/2`，完整 Bernstein 系数取 `C_k=γ^k`，γ=1、1/2，接触响应归一成q=1。其全future响应恰 `[(1−s+sγ)/(1−σ+σγ)]^n`，s≥σ时≤1，有解析全参数 ceiling；未把有限future节点误当连续cap。该有限系数模型不是原φ核或真实cube输入，未认证原FIRST/hardband/CP-GP/actual历史。

守卫选整数 `U_C≥ln(1/p_C)`（由 `2^U≥1/p_C` 及 e≥2），再以整数√上包得到认证整数κ，故它不小于正文(4)的canonical κ。对这个稍高阈值，以精确有理最优tilt z≥1直接核 `B(z)^n z^-κ≤p_C`；无需用浮点log或sqrt决定tail端点。integer common denominator重构完整Binomial尾系数，并独立核递推端点；parity相关正接受率给gate-tail≤plain-tail。随后与512等质量born children的完整 `p_C Hopp=m_C/√n` 接合，验证未归一化群tail和归一hard-column费用≤1/√n，来源总质量始终一。σ=0/1、p=0/1及包含等号的端点均有明确证书。

正文canonical阈值的全部范围由(6)–(7)解析证明；守卫并不数值认证自然对数端点或原连续路径。真实hard column N_h及原森林层数J仍是已核上游解析输入。有限数据不替代(9)的actual分配或Lebesgue空间证明，不证明剩余标签有正交通。

本轮可登记的新补支仅(11)，在当前未付实际交通依次互补支付。新严格剩余为全部原门∩lowS∩lowcoin的正诊断分配再交 `K<κ_C(x)`；低容量使阈值可能超过n，故它仍可保留全部原子交通。本轮没有用K上限推early活跃数，也没有证明完整一般√n目标。
