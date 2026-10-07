# 交叉壳预算：实际门控 queue、已付长迟延与未付接触交通

2026-10-07。解析推导；新增本文件，不修改 `gated_occupancy_bridge.md`、旧主账或数值记录。未启动新数值策略。

**本轮严格结果。** 保留 actual `g` 的径向交通有一个“出生—接触”queue 身份，且 queue 势逐点不超过一。由此得到实际接触测度的局部区间 packing `C_g(I)≤|I|+1`，但它没有支付全窗接触总量。更有用的是：真实对数迟延超过 `v` 的完整交通有源质量一次包络 `(1+B)e^{-v}H`；取 `v=log n` 得 `O(W)` 的通用互补子支。剩余是 `0≤β−s≤log n` 的短迟延 actual 交通，其平方能量仍未证。本稿给出短迟延 queue 的精确移位入流和能量测试身份，明确未知成本落在原实际接触测度与占用势变化项，未用自由 Carleson 或 HS 假设替代。

## 1. 冻结对象与技能范围

沿用 `gated_occupancy_bridge.md` 已核的 source-once 接口：原完整有限正输入 `μ=|f|dx`，质量 `W>0`；原硬高/出生子源 `ν_h≤μ`，质量 `H≤W`；原实际 `E`、`R(x)∈[a,b]`、`β(x)=n log(R(x)/a)∈[0,B]`，`B=n log(b/a)≤n log2`。在 `E` 外置 `g=0、R=a、β=0`。全部 FIRST/未来、共同捕获、出生/CP/GP、early-first-jump 完整 continuation 与源对门已经合并到实际 `0≤g_θ(x,z)≤1`；不在虚拟径向参数中重算这些门。

目标为
\[
R_g=\mathbb E_\theta\int_E\int h_{R(x)}(x-z)g_\theta(x,z)d\nu_h(z)dx.
\tag{1}
\]
`g` 可取原全部重交通或原实际深补支。若上游只给正支配，(1) 对原目标为上界；下述身份则对这个明确正支配交通成立。

已读 G01 与 F03 的 `SKILL.md`、method、cube-interface。仅使用其“先定义实际测度/层/交叉项、不能由根 packing 或形式平方免费推出嵌入”的范围守卫。本稿不调用 PDE/UR Carleson、径向 Schrödinger 或维数依赖 Schatten 定理，也没有认证新的 HS 核范数。

全体公式不要求 `R` 在连续窗口上最大。需要最大性质处另注明原允许尺度集 `\mathcal J`；它与森林层数 `J` 不同。

## 2. 每笔实际捕获的径向寿命：正测度定义

对 `x≠z` 令
\[
s(x,z)=n\log\frac{2\|x-z\|_\infty}{a},\qquad d(x,z)=\beta(x)-s(x,z).
\tag{2}
\]
核非零意味着 `d≥0`。来源原子对角线对固定来源的 `dx` 为零，再用 Tonelli 处理。定义完整实际交通的正测度
\[
dm(x,z,\theta)=W^{-1}\mathbf1_E(x)h_{R(x)}(x-z)g_\theta(x,z)
dx\,d\nu_h(z)d\Pr(\theta).
\tag{3}
\]
其总量为 `r_g=R_g/W`。由固定窗硬包络 `∫sup_{a≤R≤b}h_R=1+B`，先有
\[
\|m\|\le(1+B)H/W<\infty.
\tag{4}
\]
这只用于可积性与既有粗费，不是目标阶结算。

每笔交通在虚拟径向轴上有出生 `s` 和原赢家接触 `β`。定义
\[
C_g=\beta_\#m,
\qquad F_g(t)=\int\mathbf1_{\{s\le t<\beta\}}dm.
\tag{5}
\]
`C_g` 是真实门控接触测度，不是来源先验；其总质量仍是未知 `R_g/W`。这里采用死亡端点不含的右连续规范版本。区间长度为零的捕获同时出生和接触，没有运行占用，却仍保留两份相消的流量。

cube-cone Jacobian 给出生边缘
\[
s_\#m(dt)=\Theta_g(t)dt,
\]
\[
\Theta_g(t)=W^{-1}\mathbb E_\theta\int d\nu_h(z)\int d\varsigma_n(\omega)
g_\theta(x_{z,\omega}(t),z)e^{t-\beta(x_{z,\omega}(t))}
\mathbf1_{\{t\le\beta(x_{z,\omega}(t))\}},
\quad x_{z,\omega}(t)=z+\tfrac a2e^{t/n}\omega.
\tag{6}
\]
特别地出生边缘没有原子，`0≤Θ_g≤e^t H/W`（`t<0`）、`0≤Θ_g≤H/W≤1`（`0≤t≤B`），`t>B` 为零。

对任意 `q∈C_c^1(ℝ)`，每个寿命区间直接给
\[
-\int q'(t)F_g(t)dt
=\int[q(s)-q(\beta)]dm
=\int q(t)\Theta_g(t)dt-\int q(t)dC_g(t).
\]
所以精确分布恒等式为
\[
\boxed{dF_g=\Theta_g(t)dt-C_g,
\qquad\Theta_g(t)dt=dF_g+C_g.}
\tag{7}
\]
该证明只用正交通与区间积分，适用于有限 Borel 来源，不需对每个固定 `x` 假设半径累积质量绝对连续。对原 `L¹` 输入也可用下一节的 surface 导数得到同一身份。

## 3. 有界 queue 势、端点与接触 packing

对固定原 `(x,θ)` 定义虚拟径向响应
\[
A_t^g(x,\theta)=\int h_{ae^{t/n}}(x-z)g_\theta(x,z)d\nu_h(z).
\tag{8}
\]
`g` 固定在原实际赢家/FIRST 下，**不随虚拟 `t` 改变**。由 (5)，
\[
F_g(t)=W^{-1}\mathbb E_\theta\int_{E:\beta(x)>t}
e^{t-\beta(x)}A_t^g(x,\theta)dx.
\tag{9}
\]
以 `β≥t` 写同一式会给出死亡接触处的左极限版本；二者对 `dt` 几乎处处相等，分布身份相同。不能在原子接触时任意混用两个点值。

因 `A_t^g≤h_{ae^{t/n}}*μ`，该完整固定半径核积分为 `W`，得到
\[
0\le F_g(t)\le1\quad(t\in\mathbb R),\qquad
F_g(t)\le e^t\quad(t<0),\qquad F_g(t)=0\quad(t\ge B).
\tag{10}
\]
有限变差由 (4)、(7) 保证。`F_g(-∞)=0`，负核心满足
\[
\int_{-\infty}^0\Theta_g(t)dt=F_g(0-)\le H/W\le1.
\tag{11}
\]
若 `C_g` 在 `0` 或 `B` 有原子，它使 `F_g` 在该点向下跳跃，跳量正好是该原子质量。全轴积分给
\[
\int_{-\infty}^B\Theta_g(t)dt=C_g([0,B])=R_g/W.
\tag{12}
\]
正径向部分则为 `C_g([0,B])−F_g(0-)`；这一写法保留了两端接触原子，没有把端点原子塞进 `Θ_g dt`。

对任意有限 `α<γ`，(7) 的右连续版本给
\[
C_g((\alpha,\gamma])=\int_\alpha^\gamma\Theta_g(t)dt+F_g(\alpha)-F_g(\gamma)
\le(\gamma-\alpha)+1.
\tag{13}
\]
每个接触原子质量也不超过一。这是无 `J`、无原子数、无输入熵的实际局部 packing。但它允许全窗质量为 `B+1`，**没有**提供 polylog 接触总量或平方根费用。

如只登记 `0≤F≤1`、`0≤Θ≤1` 与 (7)，抽象长区间上的 `F=0、C=Θ dt=dt` 满足这些条件，仍可有线性总质量。这里不是原输入或 hardband 反例，只说明 queue 代数本身没有排除极短寿命流量的反复进入。

## 4. surface 导数和移动窗：与交叉壳的直接关系

令 `σ_t` 为 `x=(a/2)e^{t/n}ω`、`ω∼ς_n` 的边界概率测度。对 `L¹` 输入，在分布/几乎处处意义下
\[
\partial_t A_t^g=S_t^g-A_t^g,
\qquad S_t^g(x,\theta)=\int g_\theta(x,z)d(\sigma_t*\nu_h)(x;z),
\tag{14}
\]
其中最后记号表示对来源保留 `g(x,z)` 的表面卷积。可用 `(ae^{t/n})^nA_t^g` 的累积 cube 质量求导证明；`∂_t(e^tA_t^g)=e^tS_t^g`。普通 `g=1、ν_h=μ` 时，就是 `S_t=σ_t*μ`。

原未加权短迟延占用 `ψ_v` 因而有精确表示
\[
\psi_v(t)=W^{-1}\mathbb E_\theta\int_{E:\ t\le\beta(x)\le t+v}S_t^g(x,\theta)dx.
\tag{15}
\]
hard 放松 `Γ_v` 是把 `g` 换成 `1_E`、`ν_h` 换成完整 `μ` 后的同一式。平方是 **不同接收点** `x,x'` 在共同 `t` 上的交叉积分，不是现有同输出来源对 Schur/透镜界。

积分一个局部径向区间 `I` 时，每个固定 `x` 只涉及 `I∩[β(x)-v,β(x)]`。若实际最大赢家确在全连续 `[a,b]` 上最大，且这个交集位于 `[0,B]`，则 `A_t(x)≤u(x)`；由 (14) 得
\[
\int_I\Gamma_v(t)dt
\le[1+\min(v,|I|)]\,W^{-1}\int_{E:\beta(x)\in I+[0,v]}u(x)dx.
\tag{16}
\]
这给出了交叉壳向**接收交通测度**的局部转换，但右侧仍是未付 hardband 交通；不能用它反向证明来源预算。

原有限允许尺度集 `\mathcal J` 的赢家只给 `A_t≤u` 于 `ae^{t/n}∈\mathcal J`，不足以在整个交集积分 (16)。不使用连续最大时，累积质量单调只给 `A_t≤e^{β-t}u`（`t≤β`），相同计算可给较粗因子 `e^v`，仍落在该接收交通测度上。连续/有限尺度接口的区别没有消失。

## 5. 已付一般互补支：长迟延的精确固定列包络

现在对**原实际交通**按 `d=β-s` 分割。取任意 `v≥0`，长迟延为 `d>v`，其闭包络上界允许 `d≥v`。对每个固定来源与方向，合法 `β∈[0,B]` 且 `β-s≥v` 时的核/Jacobian 最大值是
\[
\sup_{\substack{0\le\beta\le B\\\beta-s\ge v}}e^{s-\beta}
=\begin{cases}
e^s,&s\le-v,\\
e^{-v},&-v<s\le B-v,\\
0,&s>B-v.
\end{cases}
\tag{17}
\]
积分整条 `s` 轴严格得到
\[
\int_{\mathbb R}\sup_{\substack{0\le\beta\le B\\\beta-s\ge v}}e^{s-\beta}ds
=(1+B)e^{-v}.
\tag{18}
\]
删除实际门只建立正支配，后取来源/方向/种子积分，故
\[
\boxed{R_{g,\mathrm{long}(v)}\le(1+B)e^{-v}H\le(1+B)e^{-v}W.}
\tag{19}
\]
这包括负径向核心，不额外复制来源。它不要求 hardband、连续最大、原子数、父层数或软尾截断，但**应用对象仍是保留全部门的当前 residual 子交通**。

对 `n≥512` 取 `v=log n`，`B≤nlog2`，得到
\[
R_{g,\mathrm{long}(\log n)}\le(\log2+1/n)W.
\tag{20}
\]
因此只需继续估计原互补短迟延 `0≤d≤log n`。其几何含义是
\[
2\|x-z\|_\infty\le R(x)\le n^{1/n}\,2\|x-z\|_\infty.
\tag{21}
\]
这是一条一般输入上的真实内区/迟延子支支付，不是整个余项的闭合。新版三个已付尺度分支是否已覆盖这批交通尚未独审；只能在当前残余项内部作替代分解，不把 (20) 与已有同一交通费用重复相加。

## 6. 短迟延的移位 queue：有偿入流已知，接触仍未知

定义短迟延出生密度
\[
\Theta_{g,v}(t)=W^{-1}\mathbb E_\theta\int d\nu_h(z)\int d\varsigma_n(\omega)
g_\theta(x_{z,\omega}(t),z)e^{t-\beta(x_{z,\omega}(t))}
\mathbf1_{\{0\le\beta(x_{z,\omega}(t))-t\le v\}}.
\tag{22}
\]
精确有 `e^{-v}ψ_v≤Θ_{g,v}≤ψ_v`。令截短的运行占用为
\[
F_{g,v}(t)=\int\mathbf1_{\{\max(s,\beta-v)\le t<\beta\}}dm,
\quad D_{g,v}=(\beta-v)_\#(\mathbf1_{\{d>v\}}m).
\tag{23}
\]
长寿命交通不在真正出生 `s` 进入此 queue，而在 `β-v` 注入；短寿命交通保留真正出生。区间端点直接给
\[
\boxed{dF_{g,v}=\Theta_{g,v}(t)dt+D_{g,v}-C_g.}
\tag{24}
\]
`0≤F_{g,v}≤F_g≤1`，支持在 `[-v,B]`，且移位入流已经有一次质量支付
\[
\|D_{g,v}\|=R_{g,\mathrm{long}(v)}/W
\le(1+B)e^{-v}H/W.
\tag{25}
\]
没有把这份入流当作免费新来源。对原 `L¹` 输入，还可以明确写
\[
D_{g,v}=(\beta-v)_\#\left[
W^{-1}\mathbb E_\theta e^{-v}A_{\beta-v}^g(x,\theta)\mathbf1_E(x)dx\right].
\tag{26}
\]
有限 Borel 来源需在 (26) 用 `||x-z||∞<R(x)e^{-v/n}/2` 的**开**内 cube 来表示 `d>v`，防止把 `d=v` 的出生交通同时放进入流；原 `L¹` 来源开闭无差别。

若 `β-v∈[0,B]` 且对应尺度在原 `\mathcal J`，实际 hard 最大条件还给该入流在这个分支上被 `e^{-v}` 倍移位的完整 hard 接触测度支配。`g` 下不能直接把完整 `u` 接触测度换成更小的 `C_g`。在尺度不允许或 `β-v<0` 处，不继承这一最大条件；全域有效的 (25) 已统一支付它们。

对任意非负 `η∈C_c^1(ℝ)`，(24) 的确切能量测试是
\[
\int\eta(t)\Theta_{g,v}(t)dt
=\int\eta\,dC_g-\int\eta\,dD_{g,v}-\int\eta'(t)F_{g,v}(t)dt.
\tag{27}
\]
取平滑概率核 `κ_δ` 并令 `η=κ_δ*Θ_{g,v}`（先截断后取极限），左侧是可核查的正自相关能量。右侧明确包含：真实接触—短出生交叉项 `∫η dC_g`、已付非负入流项、以及占用势的有符号变化 `−∫η'F_{g,v}`。不能删去最后一项，也不能以 `F≤1` 宣布接触交叉项来源一次 polylog；这样会把未知 `R_g/W` 再放回预算。

当前严格能量上界仍只是
\[
\int_0^B\Theta_{g,v}(t)^2dt
\le\int_0^B\Theta_{g,v}(t)dt
\le R_{g,\mathrm{short}(v)}/W\le1+B.
\tag{28}
\]
未得到维数统一改善。新推导把剩余量放在 (27) 的实际接触交叉与势变化上，而不是重新命名一个自由 Carleson 常数。

## 7. 保留同一 hardband 与合法赢家的解析压力守卫

以下是完整正 `L¹` **hard 放松**实例，检查 `Γ_v` 的必要 `v` 依赖，不认证 FIRST 或 actual geom 门。

取
\[
f(y)=c\left(1+\frac{\|y\|_2^2}{nM^2}\right)\mathbf1_{[-M,M]^n}(y),
\qquad \lambda=c/2,\quad M=n^2b,\quad n\ge2.
\tag{29a}
\]
完整来源质量为 `W=(4c/3)(2M)^n`。对所有
`x∈[-M+b/2,M-b/2]^n`，尺度 `[a,b]` cube 均包含于源盒，直接积分平方坐标给
\[
A_R(x)=c\left[1+\frac{\|x\|_2^2}{nM^2}+\frac{R^2}{12M^2}\right].
\tag{29b}
\]
它严格随 `R` 增大，故 `R(x)=b` 是唯一连续最大赢家，也是任意含 `b` 的原有限尺度集的唯一最大赢家，不依赖 tie 选择规则。接收内盒上
\[
c<U(x)\le c\left[2-\frac bM+\frac{b^2}{3M^2}\right]<2c,
\tag{29c}
\]
因此同一 `λ` 下全部属于 `2λ<U≤4λ`。其余点保留原尺度集的可测最大选择；这里只用内部点给下界，不另删源质量。

若 `z∈[-M+b,M-b]^n`，则对所有 `0≤t≤B` 和 cone 方向，`x=z+(a/2)e^{t/n}ω` 都落在上述内部接收盒。这个源内盒的完整质量比例准确为
\[
q_M=(1-b/M)^n\frac{3+(1-b/M)^2}{4}
\ge\frac34(1-1/n).
\tag{29d}
\]
因此对 `0<v≤B`，完整来源归一化的 hard 壳满足
\[
\Gamma_v(t)\ge q_M\quad(B-v\le t\le B),
\]
\[
\int_0^B\Gamma_v(t)^2dt\ge vq_M^2
\ge\frac9{16}v(1-1/n)^2.
\tag{29}
\]
这保留同一 `λ`、hardband 与实际最大赢家。它说明无权 `Γ_v` 的预算必须允许至少线性的迟延窗宽依赖，且“接近赢家壳因此自动有 `1/n` 小因子”不成立。它**没有反驳** `C(1+v)^p log^A n`（`p≥1`）或固定 `v` 的 polylog 能量；加权实际出生还保留 `e^{t-B}`，不能用 (29) 倒写成同阶一次交通下界。

本轮没有得到保持原 hardband/winner 后违反固定 `v` polylog 能量的解析实例，也没有证明该候选。不能把缺少证明当反证。未复用或重跑任何旧压力数据。

## 8. 可回代分解与确切下一接口

对当前原实际深补支取 `v=log n`，已有精确互补分解
\[
R_g=R_{g,\mathrm{short}(\log n)}+R_{g,\mathrm{long}(\log n)}
\le R_{g,\mathrm{short}(\log n)}+(\log2+1/n)W.
\tag{30}
\]
若 `g` 指全部重余项，同式仍成立；与既有浅父组分解不同时给同一交通重复收费。将 (30) 代入提供的主账只得到一个明确 `O(W)` 补支与原短迟延剩余，并没有一般端点结论。

具体下一步是利用原实际约束证明 (27) 中同一输入的接触—出生交叉与势变化合计给
\[
R_{g,\mathrm{short}(\log n)}
\le C\sqrt n\,n^{o(1)}W+\varepsilon\lambda|E|,
\qquad\varepsilon\le49/8192,
\tag{31}
\]
或直接证明对应短迟延 `Θ_{g,v}` 的 `W` 归一化 polylog 平方能量。只用 (13) 的一维局部 packing、(16) 的未知接收交通转换、(24) 的有界 queue 势均不足；需真正使用 hardband/共同未来/第一跳与来源质量的关联，或者给势变化的有偿控制。

**最终状态。** actual source-once queue、端点原子、局部 packing、全轴长迟延固定列费及短迟延移位入流均通过。一般短迟延交叉壳的平方根预算仍未证；原允许尺度集到连续最大版本的转换仍需独立结算。没有新增数值结论或更改一般主账状态。
