# 共同软度—尺度正包络：同价来源筛选与 hard-winner 测试

2026-10-07。root 已全文独审 §1–5，通过 source-once 替换、参数范围与未付空间缺口核对；新守卫收据见 §6。本轮已读取 checkpoint 最后“逐原软来源低future筛选与必要几何”节。子任务 `get_goal` 返回 `goal:null`，故不冒称读取父任务 goal 对象；按明确交接目标继续一般 centered-cube geom `O(√n n^{o(1)})`。只编辑本稿及同前缀守卫，不修改其它作者文件。

本稿的主结果是一项严格的 source-once **替换**：把旧同参数高 Z 来源支扩成全参数固定正包络 selector，费用相同。它不证明新剩余为空或已付；低参数响应仍须与真实 source pair 几何耦合。

## 1. 查重与当前完整实际对象

本轮读取原 tex 730–865 的 FIRST-MGF/FACE/MIXED、2078–2195 的双物理尺度相关与合格列、3278–3348 的保留边界二维熵特征线、5040–5070 的 one-coordinate 正平均。旧 EC 的物理 shrink／soften 路径与硬面收费限制、旧 JM 的固定相关／移动双尺度差别都保留；不将这些已有工具重新称为主账进展。

当前是 `R_far,lowav` 的实际来源交通：完整 `μ,W`、共同 `q∈[3λ/4,λ]`、`σ(x)>1/n`、`L_s(x),R_h(x)∈[a,b]`、`a≤b≤2a`、`n≥512`，以及所有 FIRST/future cap、GOOD/score/owntrace、出生、森林、唯一 LCA、CP/GP、原分数、共同重捕获、删除门保持。保留 nonconcentration、hard fixed short shell、真实 firstjump exit、大 mark 和所选 vertices 的 far 门。

完整 soft 核与来源概率为
\[
p_x(y)=P_{\sigma(x),L_s(x)}(x-y)>0,\qquad
\rho_x(dy)=p_x(y)d\mu(y)/q,
\quad P_{\sigma,L_s}*\mu=q,
\]
\[
2\lambda<u(x)\le4\lambda,\qquad u(x)/q\le C_h=16/3.
\tag{1}
\]
在当前实际联合历史／hard source 门测度中，原软来源 `y` 边际被 `ρ_x` 支配。这是 `low_future_source_payment_20261007.md` §1 已核的 measure domination：先在原 kernel 下对任意 `μ|B` 使用正 first-exit/early 分解，再乘 `[0,1]` 门；不是参考 `2W` 放大之后恢复的性质。

完整未来 cap 仍对 `s≥σ(x),L∈[a,b]` 给 `P_{s,L}*μ(x)≤q`。原完整输入低平均输出条件也继续保留
\[
\overline P_{\sigma,L_s}^{(\alpha)}*\mu(x)/q<\varepsilon_n,
\quad\alpha=\lceil\sqrt n\rceil,
\quad\varepsilon_n=\log^{-4}(n+2).
\tag{2}
\]

## 2. 固定全参数正包络及可测性

使用已经证明的平均核
\[
\overline P_{s,L}^{(\alpha)}(v)
=\alpha\int_0^1t^{\alpha-1}
P_{s+(1-s)t,L}(v)dt.
\]
定义固定、与来源输入及 receiver 选择器无关的平移核
\[
\boxed{S_{n,a,b}(v)
=\sup_{\substack{0\le s\le1\\a\le L\le b}}
\overline P_{s,L}^{(\alpha)}(v).}
\tag{3}
\]
`s<σ(x)` 也在正包络中；这些参数**不称原 future 测试**，不声称其响应受原 cap `q` 控制。来源支付只需固定正包络的空间质量。

`s` dependence 是有限连续多项式；固定 `v` 的 `L` dependence 在 hard 面进入处采用闭核上跳值，右连续，窗口内可由 `L` 从右逼近，另加入 `b` 端点。原 `φ` 连续、physical normalization 连续。故 (3) 等于对 countable dense 软度、物理参数（含端点及 `b`）的 supremum。由此 `S(v)` Borel，可测原选择器与 `p>0` 给后续比值／selector Borel。只用可测逼近，不按点数收费，不把有限节点误认全连续 cap。

引用 `future_softness_moment_budget_20261007.md` (12) 的 **已证**全参数 column：
\[
\boxed{\int S(v)dv\le K_{n,\alpha},
\quad K_{n,\alpha}
=D_*+\log(b/a)[M_*+\tfrac12\sqrt{10n}D_*],}
\]
\[
D_*=3+\log(n/\alpha),\qquad M_*=2n/\alpha,
\qquad K_{n,\alpha}=O(\sqrt n\log n).
\tag{4}
\]
本稿不重证 Beta 峰值或 Fisher。注意这是 `sup` 放在**每个 fixed source 的核**内的 bound，足够支付选中来源；它不是 `sup` receiver 行响应仍 `≤q` 的声明。

## 3. actual source-pair 的 source-once paid majorant

对任意预先固定 `τ>0`，在当前实际交通中按原软来源标签切
\[
E_{\rm src,paid}(x)=\{y:S(x-y)/p_x(y)\ge\tau\}.
\tag{5}
\]
先积分原 hard 子概率、实际 history 与全部 `[0,1]` 门，再使用原 `y` 边际支配，得到
\[
\begin{aligned}
R_{\rm src,paid}
&\le\int_E\frac{u(x)}q\int_{\{S/p\ge\tau\}}p_x(y)d\mu(y)dx\\
&\le\frac{C_h}{\tau}\int\!\int S(x-y)d\mu(y)dx\\
&\le\boxed{\frac{C_hK_{n,\alpha}}\tau W.}
\end{aligned}
\tag{6}
\]
这是 coupled actual source-pair 的正 majorant；不假定乘门后软硬源独立，不把 `y` 换成 `z,w` 或 terminal point，不重新分配停止来源。原 `μ` 只使用一次，`q^{-1}` 只出现一次，(6) 没有额外 first-exit 因子二。

取 `τ=ε_n`，费 `O(√nlog^5(n+2))W`。由于
\[
S(x-y)\ge\overline P_{\sigma(x),L_s(x)}^{(\alpha)}(x-y),
\]
新 paid selector 包含旧 `Z_x(y)≥ε_n` 来源 selector。因此 (6) **替换** root 旧高 Z 来源费，不能额外再加同一笔 `C_hK/ε`；新增空间费用阶和吸收系数均不变。等号计入 paid。

旧“高完整平均输出”费仍保留，支付次序为其互补再用 (6)。`F̄/q≥ε` 只说明一个来源平均，不说明该输出的每个来源都 `S/p≥ε`；不能仅据 selector 包含旧高 Z，就宣称它包含整批旧高输出交通。

## 4. 同价留下的全物理／全 averaged-softness 条件

新剩余每份 actual 原软来源 `y` 满足
\[
\boxed{p_x(y)>S(x-y)/\varepsilon_n,}
\tag{7}
\]
从而对**所有** `s∈[0,1],L∈[a,b]` 同时有
\[
\overline P_{s,L}^{(\alpha)}(x-y)
<\varepsilon_n p_x(y).
\tag{8}
\]
特别地，`s=1` 的平均核就是全软 `P_{1,L}`，所以
\[
\boxed{P_{1,L}(x-y)<\varepsilon_n p_x(y)
\quad\hbox{for every physical }L,}
\tag{9}
\]
包含原真实 `R_h(x)`，包含 `L_s(x)`，也包含 full-soft kernel 自身最合适的 physical scale。这比仅在选中 `L_s` 的 low Z 强；不是增加一份昂贵全参数网格。

令 `D_opt(x,y)=log[p_x(y)/sup_{a≤L≤b}P_{1,L}(x−y)]`，则 `D_opt>log(1/ε_n)`；此来源对比不与 (4) 的包络常数 `D_*` 混用。原选中 `L_s` 的 Jensen lower bridge也继续成立，因为旧 `Z<ε` 是 (8) 的特例；不重新证明该桥。

对任意 fixed receiver 定义 `μ_{x,hot}=μ|{y:p_x(y)>S(x-y)/ε_n}`、`q_hot(x)=∫p_x dμ_{x,hot}≤q`。任意预先或 posterior 指定的正概率 test measure `η_x(ds,dL)` 都有
\[
\int\!\int\overline P_{s,L}^{(\alpha)}(x-y)
d\eta_x(s,L)d\mu_{x,hot}(y)
\le\varepsilon_n q_{\rm hot}(x).
\tag{10}
\]
可把 `η_x` 放在原 hard winner `R_h` 附近，或合法 joint future 路径；不需先猜哪个 scale 会救回来源。若 `η_x` 包含 `s<σ`，(10) 的证明是新 pointwise selector，不是原 cap。hard source 标签条件化后，只要 test 仍是 (8) 的正平均，也有相同逐来源小比值。

`μ_{x,hot}` 随 receiver 变化，**不是**可重启 FIRST 或重新领取 `W` 的 fixed source。不能将 (10) 代入任意新群 kernel，再声称其全空间质量是原 `W`。硬赢家 `m,u,R_h` 继续由完整 `μ` 定义；不换成 `μ_{x,hot}`。

## 5. 为什么尚未成为全部 remaining 的 paid majorant

需要支付的是所有原 source-pair gates 和 (7) 同时成立的
\[
\int u(x)\,d\Lambda_x^{\rm actual}(y,z,\zeta)dx,
\tag{11}
\]
其中 `Λ_x` 保留原初始 `y`、hard `z∈Q(x,R_h)`、早首跳／continuation、所有 far 与短壳资格，其 `y` 边际 `≤ρ_x`。在这一补集，(6) 使用的 domination 恰好方向相反：`p>S/ε`，不能继续拿 `S/ε` 当 paid 上界。

全物理测试 (9) 与原 hard winner 的有效接触是：同一 pair 的软来源在 `P_{1,R_h}` 中弱、硬来源在 `h_{R_h}` 中被捕获。但是这两条响应作用于**不同来源**，nonconcentration 只给原 hard 小盒捕获上限；没有已证机制将弱 soft 行转换成低 actual source-pair 交通或可付的来源重数。不能从全 future 正测试上限推“低测试则低原响应”，不能把参数平均的 soft direction 误当真实 continuation mask。

原 EC 3278–3348 的联合路径已说明：体积积分因子抵消后不能把 `n` 漂移免费变成来源费；原 JM 2078 后也已说明 fixed correlation 不能直接与随输出变化的两尺度 sup 交换。这些缺口在 (11) 中仍须解决，不用较弱 joint window average 重新包装。

本轮可严格登记的是 (6) 同价扩大 selector、(8)–(10) 的 uniform physical／averaged-softness 来源合同。该 selector 本身没有闭合空间目标。现已完成的 `joint_low_future_contract_probe_20261007.md` 对固定完整原核输入及正体积 receiver box，认证连续 rightmost FIRST／完整 future cap、hardband／唯一 hardwinner、nonconcentration、`β>v*` 与最强 `S/p<ε`；其出生、森林、CP/GP 及完整 actual 历史门仍 unknown。因此该压力输入不能称一般 geom 反例，也不能与此前四个有限候选混淆。新低 S 合同需要对完整 actual 联合几何继续检验。


## 6. 新 selector 三轮精确守卫收据

独占脚本 `joint_future_scale_exact_guard_20261007.py` 与 `joint_future_scale_exact_guard_20261007_results.json` 已实际运行，三轮 `n=512,1024,4096`，六个有限正核场景、十八个 threshold case、九个闭面 case 全部 `exact_PASS`。无随机种子，不重算已有 Beta 积分、`D/M/K`，也不模拟原 `φ` 或 actual FIRST。

有限 selector 组件采用 cyclic group 上完整 uniform 来源 `W=1`，三个物理标签、两个仿射软度端点，有限平移核均列归一。平均核对软度仍仿射，因此端点 sup 是该 **有限 surrogate** 全连续软度的精确 sup。它不是原 n 维 tensor 参数族，三个标签也不认证连续物理 cap。输入还含两个未归一化历史子权重、规范 hard source 概率及依赖 `(x,y,z,history)` 的正门，直接检查原 `y` 边际支配、actual 高 selector 收费链

\[
R_{\rm actual,paid}\le C_h\int p\,1_{S/p\ge\tau}d\mu dx
\le(C_h/\tau)\int S*d\mu dx
\le(C_h/\tau)C W,
\]

其中 `C` 是此 finite S 的 exact column mass，不冒充 (4) 的真实 `K`。同时检查旧同参数 selector 包含于新 selector、等号归 paid、剩余来源对全部 finite 参数测试的低响应。六场景分别保留 sibling-rescue 与 hard-persistent 模式，合计二十次旧低／新高标签及十四次仍低标签；这些只是说明 algebra 层扩大 selector 与保留补集都可发生，不是实际门资格样本或一般 coverage 结论。

闭面组件独立使用真实 hard centered-cube 核 `L^{-n}1_{2\|v\|_\infty\le L}`，在窗口 `[1,2]` 检查 face threshold `1,3/2,2`。内部／左端闭面用一列从右趋近的有理 `L`，exact 值单调趋近闭面值；右端 face 的 interior dense sup 为零而 `L=b` 值为 `2^{-n}>0`，所以显式加入 `b` 不能省略。有限证书核对所列有理逼近值；其真正极限等式由 `L^{-n}` 连续性解析给出，不把有限节点称连续 sup 认证。

首次执行仅在大有理数十进制序列化时触发 Python digit-limit；已改为长度分帧的 big-endian numerator/denominator SHA 与 bit 长度记录，随后完整重跑全部三轮通过。该记录修正不改变任何数学计算或比较。以上守卫不证明 (11) 的剩余空间费，也不认证真实大 N 合同的全部历史门。

root 已独立全文读取新脚本并核保存的十八条精确收费链、selector 包含及闭面 `b` 端点处理，通过；独审记录为 `root_joint_future_scale_saved_review_20261007.json`。最终数值范围限于 finite cyclic affine surrogate 的 selector／相关正门代数与真实 hard 核闭面组件，**不是原 φ 的数值样本，不认证实际 FIRST／完整历史门，不证明剩余空间费**。
