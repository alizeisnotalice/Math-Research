# Logistic 正链与 actual FIRST：共同 C 来源、条件历史及未付占用

2026-10-07。只新增本稿及同前缀表示守卫，不修改主 TXT、主账或其它 agent 文件。按父任务的 6.1-sol high 设置要求工作；当前工具不公开实际模型 ID，不能独立认证运行模型。本轮完成三轮小型有理表示守卫，不重新优化尾部 η。

**最小可操作新接口：有限多个 C 层可以由同一份原来源及同一 Poisson 跳跃场共同耦合。** 它同时包含每层的 logistic 正链，不需要每层重新领 W。把原历史通过 endpoint 条件分布保留后，得到同源总接受数的准确占用身份；尚未取得该总接受数的平方根预算。正链存在没有支付实际 R_angle。

沿用已读 [E01](/Users/zhengzhihao/.codex/skills/math-e01-occupation-measure-flow-lp/SKILL.md) 的流守恒/来源清算工作流、[M04](/Users/zhengzhihao/.codex/skills/math-m04-measurable-supremum-selection/SKILL.md) 的先冻结 selector 约定。这里仅证明有限正测度构造、RN 接受和 Tonelli 身份，不引用未核连续占用 LP 强对偶。新 Lévy 正序依据已独审的 [diagonal_levy_analytic](diagonal_levy_analytic_20261007.md) §3–5；不把其空间 probe 的有限结果当全参数证明。

## 1. 两种 P、原参数和正向 fullfuture

为避免原 fullsoft 与 square root 同名，本稿写

\[
 K^1_{r,L}=\bigotimes_{i=1}^n[(1-r)I+rG_{1,L}^{(i)}],\qquad
 S_{r,L}=(K^1_{r,L})^{1/2},\qquad
 F_{r,L}=h_L^{(n)}*K^1_{r,L}.
 \tag{1}
\]

原 actual soft 响应是 F_{σ(x),L_s(x)}μ(x)=q(x)，不是 Sμ。记

\[
 C_x=\frac{σ(x)}{(1-σ(x))L_s(x)^2},\quad
 r_C(L)=\frac{CL^2}{1+CL^2}.
 \tag{2}
\]

在 0<σ(x)<1 时，r_{C_x}(L_s(x))=σ(x)。对 L≥L_s(x)，r_{C_x}(L)≥σ(x)，且 L∈[a,b]，所以原 fullfuture **确实**约束

\[
 F_{r_{C_x}(L),L}μ(x)\le q(x),\qquad L_s(x)\le L\le b.
 \tag{3}
\]

若需要从较小 L<L_s 接入同一条曲线，此时 r_C(L)<σ，(3) 没有资格。actual hard winner R_h 可位于 L_s 两侧，不能把它自动当作正向 future 时刻。σ=0 是 identity 退化端点；σ=1 对应 C=∞，不在有限总跳率构造中当成有限 C；下节末给其无限活动场。

(3) 是响应帽。它没有给 F_{r,L}(x-y) 的点态单调序，也没有把 actual 外选的 (σ,L_s) 变成辅助链的停时；选点仍在原接收点上。

## 2. 原 Lévy 密度实际上对 C、L 分别单调

在 (2) 曲线上，单轴 S 的 Lévy density 为

\[
 \lambda_{C,L}(z)=\int_{1/L}^{\infty}
 \frac1{2π}\left[\fracπ2-
 \arctan\frac{v^2/C-\log(L^2v^2-1)}π\right]e^{-v|z|}\,dv,
 \quad \|\lambda_{C,L}\|_1=\tfrac12\log(1+CL^2).
 \tag{4}
\]

这个公式包含原密度及原 physical L。增大 C 使 v²/C 减小，phase 增大；增大 L 使 logarithm 增大、下端 support 扩张。因此

\[
 C_1\le C_2,\ L_1\le L_2
 \quad\Longrightarrow\quad
 \lambda_{C_1,L_1}(z)\le\lambda_{C_2,L_2}(z).
 \tag{5}
\]

是空间正测度序。每固定 C 沿 t=log L 的跳率为 r_C(L)，n 轴总率 nr_C(L)。固定 r 跨 L 的旧不可能商仍成立；(5) 使用的是 softness 同时增大的不同路径。

## 3. 一份来源及共同 C 场的严格有限构造

取任意有限原参数集合 θ_j=(C_j,L_j)，均有 0≤C_j<∞、a≤L_j≤b；C=0 解释为空 Lévy 测度。令 C_*=max C_j、L_*=max L_j。取一份轴向 marked Poisson 点集，强度

\[
 \sum_{i=1}^n\mathbf1_{0<u\le\lambda_{C_*,L_*}(z)}\,dz\,du,
 \qquad M_*=\frac n2\log(1+C_*L_*^2)<∞.
 \tag{6}
\]

无需使用无穷强度定理：先抽 Poisson(M_*) 个点，再按 (6)/M_* 独立抽 (i,z,u)；M_*=0 时直接取空集。从同一个原 y 出发，定义

\[
 X_j=y+\sum_{(i,z,u)}z e_i\mathbf1_{u\le\lambda_{C_j,L_j}(z)}.
 \tag{7}
\]

每个 X_j−y 的 law 正是 S_{r_{C_j}(L_j),L_j}。若两参数在 (5) 的序上，新增点落在不交的 marked 区域，故是独立、对称、保质量的真实增量。特别每固定 C 子序列是原 logistic Markov 链。所有 C 层均从同一 y 和同一份点集产生，来源没有复制。

不需要断言矩形混合差 λ_{C2,L2}−λ_{C1,L2}−λ_{C2,L1}+λ_{C1,L1}≥0；(7) 也不提供这种结论。不可比的参数点没有一条已证明的有序一参数链。该共同场是耦合，不是一个可以任意排序后套 Doob 的鞅。

从原 μ/W 抽 y 一次后构造 (7)，总概率质量仍是 1；不是按 C_j 分别发放 W。这个改进解决共同耦合的存在性，不是跨层次数的费用。

任意有限核族本就可以给定 y 后独立抽各 endpoint，从而共用一份来源；所以“共源 coupling 存在”本身不是新上界。这里新增的结构是 (5) 偏序上的**原 S 正增量相关**，尤其固定 C 的原 logistic 链关系；必须实际使用该相关才能期待改进预算。

全参数场也可包括 σ=1。analytic 稿 §3 证明

\[
 \lambda_{\infty,b}(z)\le\frac{e^{-|z|/b}}{2|z|},\qquad
 \int |z|\lambda_{\infty,b}(z)dz\le b.
 \tag{7a}
\]

对每轴取强度 λ_{∞,b}(z)dz du、0<u<1 的 Poisson 点集，参数 θ=(C,L) 选择 u≤λ_{C,L}(z)/λ_{∞,b}(z)；C=0 选空集。可先按不交 dyadic 空间环构造独立有限 Poisson 点集。总绝对跳和期望≤nb，故几乎处处有限，**所有参数子和同时绝对收敛**。定义每点对应的同一个 (7) 即可；任意嵌套选择集的增量独立。有限 r 的边际仍为原 S；C→∞ 的对称特征指数由 |e^{iξz}−1|≤|ξz| 和 (7a) 控制，得到原 S_{1,L}。没有把无限活动过程称作有限总率过程。比值和这些绝对收敛子和有共同可测版本，因而全参数 coupling 存在。

另有更贴近中心几何的 coupling：取两份独立的同类 marked fields Z_θ、Z'_θ（各为原 S 噪声），再取一份共同 U∼Uniform([-1/2,1/2]^n)，定义

\[
 Y_θ=y+Z_θ+Z'_θ+LU.
 \tag{7b}
\]

每个 θ 的边际正是原 fullsoft p_θ。对固定有限 C，沿 t=log L，扩张态 (Y_t,U) 是真正非齐次 Markov 过程；对光滑测试函数其 generator 为

\[
 \mathcal A_t\varphi(Y,U)=LU\cdot\nabla_Y\varphi(Y,U)
 +2\sum_i\int[\varphi(Y+ze_i,U)-\varphi(Y,U)]
       \partial_t\lambda_{C,L}(z)dz.
 \tag{7c}
\]

跳率为 2nr_C(L)，U 保持冻结。两个独立原 S 增量相加给第二项；同一 LU 的变化给第一项。保留 U 是 Markov 性的条件，不把 Y 单独的后验当成同一个 convolution 增量。C=0 时只是 y+LU 的 pure-hard 几何；这已经说明 augmented coupling 存在本身不能给目标弱界。共同 U 不是原 LCA 的共享 grid seed θ，也不能替代原硬来源 z。其 sample paths 给分段仿射径向几何，但原门仍须第 5 节的 endpoint RN。

## 4. 原 early first jump 和 c_t 只能条件接入

原物理尺度 L 下，首标签 t、c_t=1−t、p_t=(σ−t)/(1−t) 的 continuation square root 满足

\[
 S^{c_t}_{p_t,L}=\frac{S^1_{σ,L}}{S^1_{t,L}},\qquad
 \Lambda^{c_t}_{p_t,L}=\Lambda^1_{σ,L}-\Lambda^1_{t,L}.
 \tag{8}
\]

同一 L 时这是原正有序增量。移 L 会改变分母；两个正链之差不推出条件 continuation 的正链。其 effective odds 参数

\[
 \frac{p_t}{(1-p_t)L_s^2}
 =\frac{σ-t}{(1-σ)L_s^2}
 =C_x-\frac{t}{(1-σ)L_s^2}
 \tag{9}
\]

还依原首标签 t。不得把 c_t 重新设为 1、把 t 抹掉，或者重启首跳。新 analytic 稿 §5 给出的 spectral 导数有负高频系数，严格说明 c1 的逐 coefficient 正证明不能直接搬来；它尚未证明条件空间 Lévy 序失败。

合法接入的做法是：**原完整输入和 selector 先冻结，保留原路径的 endpoint 条件分布。** 原 (1−t)^{n−1}G_{1−t,L}(ρ)、早时钟、首次落点 w=y+ρe_i、退出门，以及完整 h_L*K_{σ,t,L} 均继续属于原物理路径。新场 (7) 的跳跃时间或坐标不能命名为这些原标签。

## 5. 有限 actual rows 的同源 RN 接受：保留所有门

这一节限于参数恰属于有限集合的 actual rows，或一个已经独立合法化的实际有限网合同。令不交 E_j 是原接收 rows 的参数分区；不是重新计算一个自由 family 的赢家。原核 k_j(x,y;μ) 在 E_j 外置零，在 E_j 内保留原 current 高交通 k_hi。它仍包含 FIRST/fullfuture、η_C、出生/own/GOOD/score/clock、source cuts、原硬 z、远对、共享 seed、唯一 distinct-child LCA、失败权、strict CP/GP、重捕获、early/history/continuation、双 cutoff/core 和当前角门补集。

原正首跳子核、Σ_Cν_C≤μ、唯一 LCA 及硬端积分≤u，给既有上界（q≥3λ/4、u≤4λ）

\[
 0\le k_j(x,y;μ)\le C_h\kappa_j(x-y),\quad
 C_h=16/3,\quad \kappa_j=h_{L_j}^{(n)}*S_j*S_j.
 \tag{10}
\]

这里只用它定义接受，不把丢门后的上界重新登记为费用。取原 μ(dy)dx joint measure 上的 RN 代表

\[
 a_j(y,x)=\frac{k_j(x,y;μ)}{C_h\kappa_j(x-y)}\in[0,1],
 \tag{11}
\]

零 prior 处取 0。冻结 finite selector 后定义它，避免任意逐参数代表沿不可数 graph 的问题。若需输出原首标签、硬来源 z 或其它历史，接受后再按**原 accepted joint measure 给定 (y,x,j)** 的条件分布抽取；这保留联合测度、t/c_t continuation 与所有原门。它不声称辅助场与原 ON/FIRST 路径逐跳相同。

准确范围是每个固定 index 内的原联合边际；不继承原跨 index 的历史 copula 或 filtration，也不证明原 FIRST 是 (7) 的停时。连续 selector 迁移还需一致可测的参数核/RN 版本及积分论证，不能沿未定义的 argmax graph 替值。

共享 θ、唯一 LCA 和 CP/GP 必须按原 y,z 及原森林判断；不能按 X_j 或新 receiver Y_j 判断。原 gate 若还含完整 μ 的自适应量，在本次抽样中保持冻结。

从 μ/W 抽 y，再按 (7) 得 X_j；给每个 j 独立抽第二个 S_j 增量及 centered h_{L_j} sensor，得

\[
 X_j\xrightarrow{S_j} Z_j\xrightarrow{h_{L_j}}Y_j.
 \tag{12}
\]

边际 Y_j|y 的密度正是 κ_j(x−y)。sensor 即使随 j 变化也合法；它仍需中心立方体几何来控总接受次数。原另一个 hard winner R_h 和来源 z 已保留在 (11) 的原条件历史中，(12) 的 h_{L_j} 是 fullsoft 自身的 outer sensor，二者不能混称。

## 6. 共同 C 占用的准确预算格式

记 T(x)=Σ_j∫k_j(x,y)μ(dy)，限定上述 rows。对其真正高阈值 E_τ={T>τ}，取

\[
 g(x)=\frac{τ}{T(x)}\mathbf1_{E_τ}(x),\quad
 B_j=\mathbf1_{U_j\le g(Y_j)a_j(y,Y_j)},\quad
 N=\sum_j B_j.
 \tag{13}
\]

独立 receiver/history/coin 的构造保证，给定共同场和 y 后各 B_j 条件独立；不宣称无条件独立。每个 conditional acceptance 依 y，不能冒充原 Lebesgue reverse martingale 的投影。正 Fubini 给

\[
 \boxed{\mathbb E_\mu N=\frac{τ|E_τ|}{C_hW}.}
 \tag{14}
\]

不作弱归一化、只把 g 换为原 E 的指示时，同样有 C_hW E N=∫_E T；Palm receiver 此时边际是 T(x)dx/∫_E T，**不是** uniform 原 E。

因此 actual T 未给逐行 lower threshold 时，只能用该强交通校准和原主账；不能把 (14) 的 τ|E_τ| 写成原 λ|E|。本稿来源为完整原 μ/W；旧 fixed 障碍 ν_b 的 S u=ν_b−μ_b、μ_b≤κ 不自动适用于 μ。若用 signed 障碍势，必须另给 actual 来源映射，不能合并两种来源的 W。

对 (13) 则 accepted-occurrence Palm receiver 真为 dx/|E_τ|；定义 P#(dω,j)=B_jP_μ(dω)/E_μN，有

\[
 \mathbb E^\#\frac1N=\frac{\Pr_\mu(N>0)}{\mathbb E_\mu N},\qquad
 τ|E_τ|=C_hW\Pr(N>0)+C_hW\mathbb E(N-1)_+.
 \tag{15}
\]

令 c_m(y)=Pr(N≥m|原来源=y)，则全 C 层共用的精确来源占用为

\[
 C_hW\mathbb E(N-1)_+
 =C_h\int\sum_{m=2}^J c_m(y)\,μ(dy).
 \tag{16}
\]

first-accept 来源 μ(dy)Pr(N>0|y)≤μ(dy)；总质量≤W。其它接受不能另领 W。可检验的充分合同例如

\[
 \Pr^\#(N\le K_n)\ge1-δ
 \quad\Longrightarrow\quad
 τ|E_τ|\le\frac{C_hK_n}{1-δ}W.
 \tag{17}
\]

希望 K_n=√n n^{o(1)} 仍是**待证**。它必须针对 (11) 的原实际门以及 (7),(12) 的 centered sensors 证明；一般 bistochastic sensor 的旧反射反例排除了免费小 count。直接平均每个 C 的 conditional weak 界不能推出 (17)。

也可用 (7b) 替代 (12)：给定两份全 fields、共同 U、原 y 后，直接独立抽 coin，成功概率为 g(Y_{θ_j})a_j(y,Y_{θ_j})。同样的 marginal/Tonelli 证明仍给 (14)–(17)，但 N 的 joint law 和 tail 会改变。旧 q_j=P_jH_jg_j 递推及其 filtration/future-union 不能直接搬到共享 U 的新 coupling；更不能沿 C bands 套一条不存在的全序停止过程。两种 coupling 的 mean 相同，不代表尾预算相同。小守卫未把共享 U 当原空间 sample。

若把 N=Σ_bN_b 按 C bands 分组，则 (16) 保留所有跨 band 重复。只按 band 分别停止，所得源 μ_b^stop=μ Pr(N_b>0|y) 一般满足 Σ_b μ_b^stop>μ。给它们先分配权重 w_b(y)、Σ_b w_b≤1 虽可付一次源，但尚无证据这些权重支配所有 original receiver rows；这正是缺失的共同 C 分层预算。

## 7. 连续 actual 参数仍需具体覆盖或占用引理

(7a) 已给全参数共同场，故连续 C 及 C=∞ 本身不是概率耦合的存在性障碍。但是不能对不可数 receiver 参数逐个计一个 Bernoulli 接受，再把 N 当可积 count。需要先证明有效 finite covering/停止结构，或建立可积占用测度及同源尾合同。

对原 μ=fdy、f≥0、f∈L1 及有限 |E|，存在一个**不重判 gate 的边际交通逼近**。先在 selected joint measure 上定义原 a_selected(y,x)，使

\[
 k_selected(x,y)=C_h a_selected(y,x)p_{σ(x),L_s(x)}(x-y),\quad0\le a_selected\le1.
 \tag{18}
\]

取 measurable finite-valued simple maps (σ_m,L_m)→(σ,L_s)，L_m∈[a,b]，σ_m<1，σ=1 用 σ_m↑1；定义 k_m=C_h a_selected p_{σ_m,L_m}。同一 a_selected 中仍编码原 FIRST 和全部原资格，不沿每个新参数重选 RN。

这里原 fullsoft 密度准确为

\[
 p_{σ,L}(d)=\prod_i\{(1-σ)h_L(d_i)+σ(h_L*w_L)(d_i)\},
 \quad0\le p_{σ,L}\le a^{-n}.
 \tag{19}
\]

active 因子为 L^{-1}(h*w)(d_i/L)，连续依参数；inactive 因子只可能在 |d_i|=L/2 失去连续。对每固定 x，selected L_s(x) 的这些 y 超平面具有 μ 零测度。由 Fubini，它们在 dx|E×μ 下零测；σ 系数为 polynomial，两个端点也连续。因此 dominated convergence 严格给

\[
 \int_E\!\int |k_m-k_selected|\,μ(dy)dx\longrightarrow0,
 \tag{20}
\]

外包为 2C_h a^{-n}、总底测度 |E|W<∞。每个 m 的参数集合有限、最大 C 有限，可用 (6)；其 pointwise source 接受仍用同一 a_selected。固定 index 内原 conditional history可作为标记保留，但新 k_m 的 endpoint marginal暂是近似，只有 (20) 恢复原 selected accepted joint traffic。该论证需要原 L1 绝对连续来源，不能推广给 arbitrary atomic μ。

(20) 没有证明 count/Palm 尾对网格一致。σ_m 或 L_m 的逼近也未被宣布为新的实际 FIRST；原 fullfuture 资格留在冻结标记中，不能未经证明认为每个逼近参数继承 exact cap。若以后给 finite fields 的统一费用，须使该引理适用于这些冻结门的近似 rows，或单列帽误差迁移。σ=1 的 marginal 逼近不依有限 intensity 一致界。

真正下一引理可取以下形式：对当前原实际中央 rows，构造有限/可积候选集合，保留原 conditional accepted joint measure，证明其总 accepted-occurrence Palm 下 (17)，并且常数不依 C bands 数、网格密度或 θ。必要测试要同时保存 μ、original selector、q/u、L_s/R_h、C_x、原 t/c_t、CP/GP、LCA/seed、完整 history 接受权，及跨 C 的同源 count。未构造这些实际数据前，不新增 toy 或旧尾数值。

fullfuture 在这里新增的真实限制只有 (3)：沿原 C_x 的向前曲线，全部完整响应受 q 帽。若欲由它推 (17)，还须证明该响应帽及实际门如何抑制一次接受后的 centered 未来/跨 C 接受。现有 (3) 没有这个条件概率结论。

## 8. 本轮消除的缺口与保留的欠项

1. 已消除：c1 的 fixed-r 移 L 无正商障碍，可改用原 logistic 曲线；全参数 receiver-varying C（含 C=∞）也有同源共同正场，不需每层复制 W。共享 uniform sensor 还给原 fullsoft 的扩张态 Markov 过程 (7c)。
2. 已明确接入：actual FIRST 点是其 C_x 曲线的起点；较大 L 的链参数落在原 fullfuture。原 first jump、c_t、hard z、CP/GP/LCA/allhistory 可以通过冻结 selector 后的 accepted endpoint 条件分布精确保留。
3. 仍未消除：条件首跳后移 L 的直接正链、原外选参数的适应性、continuous actual graph 的有效覆盖及 count 网格一致性、centered moving hard sensor 与跨 C 的 total count 尾预算。共同场/一次源身份及 (20) 的边际逼近均不能替代这些预算。

本稿没有给一般 √n 结论，也没有新增主账 paid 费用。实际 R_angle 及原一般 O(n log n) 基线保持原状态。

## 9. 三轮新表示守卫及证据范围

应用并读取 [L03](/Users/zhengzhihao/.codex/skills/math-l03-rational-and-interval-certificates/SKILL.md) 及 method/provenance/cube-interface，只用 Fraction 精确有限代数。先写 registration，再生成 guard，完成静态核对后只执行一次；脚本拒绝覆盖已有终态。没有重跑旧 2512 项，没有抽样或空间积分。

三轮分别使用 4/5/6 个独立 Bernoulli cells，p_i=1/(i+2)。四参数 subset 为 A00={0}、A01={0,1} 加偶数 extra cells、A10={0,2} 加奇数 extra cells、A11=全部；这只是 marked-field 的有限概率代数 fixture，**不是原 Lévy 点集或原空间来源**。每轮精确枚举 16/32/64 field patterns。核每个 marginal PGF、嵌套旧/新增 counts 独立、不可比参数的共享 intersection 及 covariance。

来源 masses 为 2,3（W=5），四 receiver states 及 acceptance 按登记的有理表给定。核强交通 E N=∫T/(C_hW)、弱归一化 E N=τ|E_τ|/(C_hW)、Palm receiver uniform、source-once stop、unpaid deficit、tail sum 和 Palm 倒数 count；Ch=16/3。没有中心立方体、actual FIRST、CP/GP/LCA 或原 history 样本；不能从这些 fixture 的 count 数值估中央维数阶。

终态为 **479/879/1667 项，合计 3025 项 PASS**；弱均值三轮均为 57/350（这是登记表的校准结果，不是渐近阶）。文件为同前缀 `_registration.json`、`_guard.py`、`_results.json`、`_receipt.json`。保存 SHA-256：registration `f4fa9bc402c15a9d80161f21baa1e68a88098a5ea31e1222fdd7710b6f7c2367`；script `5e2d23f2472bf3aead645be19a6cd7050a5c1521077fde5dd4908bad9aa0fe7f`；results `ab74ff919e446e9f0741c3704057b9051692089e366b2c14aac780f386a37254`。只认证表示代数，paid_fee=false。
