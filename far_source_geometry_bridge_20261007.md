# 实际远源与同一初始 cube：几何桥核查

2026-10-07，独立一般理论探索。只新增本稿；没有改动旧理论、脚本、数据或台账，没有运行数值实验。全文读了 `nonconcentrated_actual_far_interface_20261007.md` 和主 prompt《中心立方体_geom余项最终估计_Pro提示词_20261006.txt》，并核对旧 tex 的 first-exit、continuation 正分解及共同 full future ceiling 段落。最新上游三个尺度原件缺失仍作为依赖保留；这里不补造其定义。

**结论。** 得到精确的 L_s/R 几何分支、一个否定无条件 response 分离的局部配置、合法的固定参数边界交叉积分，以及完整 future cap 给出的同源 signed 行证书。它们尚未支付 `R_far,new` 的移动空间列费。没有新费用登记，也没有将局部配置宣称为完整 FIRST/GOOD/出生/CP/GP 反例。

## 1. 冻结实际对象

完整原正输入 μ、窗口 a≤b≤2a、共同 q∈[3λ/4,λ]、实际 σ,L_s,R 及所有原门均冻结。在 residual 上已经有

\[
\mu(B\cap Q(x,R))\le\eta_n m,\quad
h=2a/n,\quad \eta_n=49/(65536\sqrt n),\quad
m=\mu(Q(x,R)),\quad u=m/R^n.
\]

还保留 β=n log(R/a)≥v_*、0≤d=n log(R/(2\|x-z\|_\infty))≤v_*、σ>1/n、首跳真实退出且 mark 大于 a/n，以及已经登记的 selected-vertex-far 条件。记

\[
R_-=R e^{-v_*/n},\qquad
R_-/2\le\|x-z\|_\infty\le R/2.
\tag{1}
\]

原正 continuation realization 的 terminal aggregate 点为

\[
V=w+\sum_{j\in A}\xi_j e_j,\qquad
U=x-V\in[-L_s/2,L_s/2]^n.
\tag{2}
\]

U 是 (21) 中同一次初始 cube 偏移。给定 x 和这份 soft realization 后，它已经确定，不能再独立均匀重抽。原 soft 来源是 y；V 是历史参数产生的 terminal 点，不是免费替换的原来源。`A_f≤1` 中完整的 actual gates 保留；软、硬子来源也不重新归一。

## 2. L_s/R 的精确三分支

令 D=z−V，则 x−z=U−D。同一初始 cube 的 hard/terminal 响应差为

\[
\Delta_L(x;V,z)=h_L(x-V)-h_L(x-z)
=L^{-n}\mathbf1_{\{U-D\notin[-L/2,L/2]^n\}},
\quad L=L_s,
\tag{3}
\]

其中 (2) 已保证第一项存在。逐坐标 crossing 的精确条件是存在 i 使

\[
U_i-D_i>L/2\quad\text{或}\quad U_i-D_i<-L/2.
\tag{4}
\]

由此只有以下无条件结论：

- **L<R_-：** (1) 强迫 z 不在同一个初始 L cube，故 Δ_L=L^{-n}。等号 L=R_- 不能并入：d=v_* 时 z 可恰在闭 cube 面上。
- **L≥R：** 全部 hard z 都在同一个初始 L cube，故 Δ_L=0；距离 \|z−V\|∞>a/n 不改变这一结论。
- **R_-≤L<R：** 必须按实际坐标判断 (4)。该分支既允许 Δ_L=0，也允许 Δ_L=L^{-n}，单凭 far 与短壳不能确定。

短壳只控制 hard 源到 receiver 面的径向余量

\[
0\le R/2-\|x-z\|_\infty
\le\frac R2(1-e^{-v_*/n})\le\frac{av_*}{n}.
\tag{5}
\]

远离 V 的最大坐标可以是另一个坐标。甚至在同一个坐标上，两个点仍可位于初始 cube 内、远离彼此。因此不能将 (5) 当作 terminal-to-hard displacement 的法向正交叉下界。

## 3. 局部零 response 配置：允许不同比例 L_s/R

取 a=1,b=2,n≥512，固定 R=3/2、ε=1/(4n)，并令

\[
x=0,\qquad c=\frac R2(1-\varepsilon),\qquad
z=c e_1,\quad V=c e_2,\quad
y=V+\tfrac14e_3,\quad w=V.
\tag{6}
\]

首跳沿 e_3 的 mark 为 −1/4，真实大于 a/n；y 与 w 相距远大于原格侧长 a/n，所以首跳退出源格。选择没有后续 active coordinate 的正 continuation 分量 A=∅，terminal 正是 V。取 σ=2/n 和 0<t<σ/2，则这是原 first-jump/no-continuation 参数化中有正系数的局部配置，且 t<T_n。σ>1/n，并不在另一个 agent 核查的 σ 极接近 1 已付区。

这里 β=n log(3/2)>v_*，d=−n log(1−ε)<1<v_*。硬 z 到 y,w,V 三个点的 infinity 距离均至少 c>a/n。在小软度支中 A=∅，所选列表只有 y,w，terminal V=w；在其它支中三点列表的第三点也是 w。因此两种原列表的 far 几何都满足，没有另外增加顶点。

若 L_s=2，则 L_s/R=4/3，z 和 V 均在初始 cube 的严格内部，(3) 的差为零。若改取

\[
L_s=R(1-\varepsilon/2)<R,
\]

仍有 R_-<L_s<R 且 c<L_s/2，差仍为零。故这个失败不限于 L_s=R。初始 cube 包含、far、短壳、真实 first exit 均有严格余量；这些纯几何条件允许附近的正体积局部配置，而不是只有闭面等号。

此外 z 与 V 的坐标值只差置换。原 h_L 和全部同参数的对称张量软核 P_{s,L} 都满足坐标置换不变性，所以在 (6) 上

\[
h_L(x-z)=h_L(x-V),\qquad
P_{s,L}(x-z)=P_{s,L}(x-V)
\quad\text{对所有 }s,L.
\tag{7}
\]

这也否定“换成原对称张量平滑响应，就能仅由距离迫使非零 response 差”的纯几何命题。它不否定使用方向标签、整体同源平方、输入-dependent cancellation 或其它完整门后才成立的证书。

**范围限制。** (6) 没有提供实现完整 nonconc、FIRST、全部 future cap、GOOD、score、出生与森林条件的输入。它证明的是上述局部几何蕴含不成立，不能称为原实际 far 交通的反例，更不是一般弱型目标的反例。要用完整门排除该局部配置，仍须写出那一个门的具体不等式，而不是从 far 自动推断。

## 4. 边界 crossing 的合法积分与其收费缺口

先仅固定 L>0,V,z。用 Lebesgue 换元 x=V+U，初始 cube 密度下的 crossing 质量精确为

\[
\int h_L(x-V)\mathbf1_{\{z\notin Q(x,L)\}}dx
=1-\prod_{i=1}^n(1-|D_i|/L)_+
\le\min\{1,\|D\|_1/L\}.
\tag{8}
\]

这是固定空间核的公式；它不允许给定实际 x 后重新均匀化 U。far 给的是距离下界，不给 (8) 所需的小 \ell^1 上界。远距离时右边可为 1。

保留不同比例 R 后，固定核 column 的准确几何量是

\[
\int h_L(x-V)h_R(x-z)dx
=L^{-n}R^{-n}\prod_i\ell_i(L,R,D_i),
\tag{9}
\]

\[
\ell_i=\bigl[\min(L/2,D_i+R/2)
-\max(-L/2,D_i-R/2)\bigr]_+.
\]

若还限制短壳，积分区域须从 Q(0,L)∩(D+Q(0,R)) 删除 D+Q_open(0,R_-)；没有改用相同尺度。因 v_* 固定，完整 hard shell 的体积分数是 1−e^{−v_*}，接近 1，不能从 side thickness O(1/n) 直接声称体积 O(1/n)。

将这些固定 column 积分用于实际交通尚需支付选择：L_s,σ 依赖 x，V 的历史分布也随原 L_s,σ 变化；同时 z 的 hard 子源及 A_f 随同一个 x 被选择。现有逐 t continuation 正包络的 N_hD_n 正费仍在，(8) 没有无条件缩小它。这里没有以未证 conditional independence 或每个输出新一份来源来结算这笔费。

## 5. 一个真正来自共同 future cap 的同源 signed 行证书

共同 future cap 可以给出比局部距离更强的 **完整行平均**，但它尚不能直接绑定选出的 far pair。使用原 hard winner 的物理尺度 R 作为合法辅助 future 参数，保持原 soft winner L_s 不变。原 full cap 给

\[
F_R=P_{\sigma,R}*\mu(x)\le q.
\]

记原噪声概率

\[
N_{\sigma,R}=[(1-\sigma)\delta_0+\sigma w_R]^{\otimes n},
\qquad P_{\sigma,R}=h_R*N_{\sigma,R}.
\]

对每个完整 hard 来源 z∈Q(x,R)，定义辅助独立噪声下的退出概率

\[
e_x(z)=N_{\sigma,R}\{\zeta:z+\zeta\notin Q(x,R)\}\in[0,1].
\]

因为 inside μ 的 soft 响应≤完整 μ 的 soft 响应，有

\[
\int r_x(dz)(1-e_x(z))
=\frac1u\int_{Q(x,R)}P_{\sigma,R}(x-z)d\mu(z)
\le\frac{F_R}{u}\le\frac qu<\frac12.
\]

从而

\[
\boxed{\int e_x(z)r_x(dz)\ge1-q/u>1/2.}
\tag{10}
\]

等价的完整同源 signed 响应证书是

\[
\boxed{\int[h_R(x-z)-P_{\sigma,R}(x-z)]d\mu(z)
=u-F_R\ge u-q>u/2.}
\tag{11}
\]

所有混合来源在 (10)–(11) 中一起保留；没有把某个来源分量重选为 FIRST。辅助尺度 R 在合法 [a,b] 中，σ 也未被重选。这个新噪声用于行证书，不被宣称为原 soft path continuation 或原 terminal V 的同一次偏移。

它当前缺两座桥。首先，原实际 A_f 可以偏重于 e_x(z) 较小的 hard 子来源；(10) 的整体均值不能升级成每个 z 或 A_f-selected posterior 的相同下界。实际子源不归一、行门正支配，并不能修复这种下界方向。其次，对固定 R,σ，signed 核积分为零、正部质量至多 1；对 moving R(x),σ(x) 取正部再空间积分仍需 source-once 列证书。不能以 (11) 为理由免费调用其未知弱型费。

同样，future cap 只对原完整 μ 成立。将其写到 terminal V 的 pushforward、停止后的来源或任意 pair-selected 子来源，须先证明相应响应上限。辅助 endpoint 分布依赖同一原历史和输出，不能因核质量为 1 就改称与 μ 相同的来源测度。

## 6. 当前可接主账的状态与最小剩余桥

本稿没有新增已付子交通，因此不改变既有互补分解或吸收系数，不重复支付大 σ 已付区。几何三分支可直接作为未来证书的有效约束：(3) 仅在 L_s<R_- 强制正初始 response crossing，另两支保留实际坐标配置。

若要从共同 future cap 构造可支付空间对偶，最小缺口是下述二者之一：

1. 将实际 far 子交通与 (10)–(11) 的完整同源 signed 损失有偿耦合，保留 A_f、mixed 来源、原 L_s/R、first-exit 和完整 continuation；并支付未耦合补集。
2. 在 L_s<R_- 的真实 initial-cube crossing 子支上直接建立固定来源的 moving σ,L_s,R 列费；其余 L_s≥R_- 支仍须另付，不能按距离删除。

两者都需要 O(√n n^{o(1)})W 的空间成本，或当前剩余额度内的吸收小行项。这里只有行恒等式和固定参数积分，未取得该成本，也未证明它不可能。为避免忽略原门，本轮没有自行做数值压力；如果安排实例，必须先构造完整 μ 并核 full FIRST/future cap，再计算实际 A_f 与所需 signed 配对，单测 (6) 的几何不能替代该任务。
