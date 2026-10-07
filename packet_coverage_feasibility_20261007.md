# 预固定 source packets 的常数覆盖要求过强

2026-10-07。本稿独立审计 [无背景 packet square 引理](general_source_packet_square_20261007.md) 后面的覆盖要求。采用 L02 的真实输入重建流程；不借用其待验收文献结论。模型沿用父任务指定的 6.1-sol high；运行时未另暴露可独立核验的模型标识。

**结论：固定 η>0、θ=n^(-1/2) 时，不存在适用于所有完整正 L¹ 输入、所有真实 winner 的维数无关正平均覆盖常数。** 下文的反证量词覆盖所有预固定、不重复来源的 fractional packets，允许分配依赖整份输入和阈值。第一输入近均匀，原连续尺度窗口上核心 winner 唯一。第二输入的全部正来源都在高度 τsqrt(n) 以上，原有限尺度 {1,2} 上核心 winner 仍唯一；同样否定覆盖。两个输入的弱比均趋于1，故都不是弱型端点反例。

与既有 [shift-averaged 小网格稿](shift_averaged_source_dominance_20261007.md) 不同，本稿不预设 packets 是网格小盒、固定几何直径或硬分区；覆盖上包来自每份 packet 的捕获集合自身。也不按 receiver 重组来源。

## 1. 原对象与被审计的合同

记完整边长方体 Q(x,R)=x+[-R/2,R/2]^n。原响应为 μ(Q(x,R))/R^n，μ=f dy，W=μ(R^n)。尺度窗口包含 a，且所有尺度属于 [a,b]。取原真赢家 R(x)，E={Mμ>τ}，m(x)=μ(Q(x,R(x)))，I₁=τ|E|。

任意预固定有限或可数 fractional packets 满足

\[
 d\mu_i=a_i\,d\mu,\quad a_i\ge0,\quad \sum_i a_i\le1,
 \qquad M_i=\mu_i(\mathbb R^n),\quad\sum_i M_i\le W. \tag{1}
\]

零质量 labels 可删。令 m_i(x)=μ_i(Q(x,R(x)))，并取

\[
 E_i=\{x\in E:m_i(x)\ge\theta m(x),\ m_i(x)\ge\eta M_i\},
 \quad
 q_{\rm dom}(x)=\sum_i\mathbf1_{E_i}(x)\frac{m_i(x)}{m(x)}\le1. \tag{2}
\]

已有 source-square 引理给该支 ∫S_dom²dμ≤4I₁/(ηθ)。若额外有 ∫E q_dom≥κ|E|，Cauchy–Schwarz 才给 I₁≤4W/(ηθκ²)。本稿审计的正是这项额外覆盖，未改变 (1)–(2)。一个 label 满足 dominance 仅给 q_dom≥θ，并不自动给 κ 常数。

## 2. 保留 density cap 的捕获集合 packing

先证明一个一般局部事实。设 μ≤H_cap dy，F⊂E 上原 winner 固定为 a，且 m(x)≥M₀a^n，M₀>0。对 F_i=F∩E_i，非空就有

\[
 M_i\ge\theta M_0a^n. \tag{3}
\]

对固定 i，把 μ_i/M_i 写成概率 ν_i。每个 x∈F_i 都满足 ν_i(Q(x,a))≥η。若 k 个代表中心的两两交集满足 ν_i(Q_l∩Q_j)<η²/2，令 T=Σ_l1_Q_l，则

\[
 k^2\eta^2\le(\mathbb E_{\nu_i}T)^2
 \le\mathbb E_{\nu_i}T^2
 \le k+\frac{k(k-1)\eta^2}{2}. \tag{4}
\]

特别是代表数不能超过 K_η=ceil(2/η²)。贪心增加代表，只在找到对所有已选代表交集都小于 η²/2 的中心时增加；有限次后必停。因此每个 x∈F_i 与至少一个代表 x_l 满足

\[
 \mu_i(Q(x,a)\cap Q(x_l,a))\ge\eta^2M_i/2.
\]

这里选代表是证明已经固定的 F_i 的覆盖，不重新定义 μ_i。由 μ_i≤H_cap dy、(3) 及立方体交集体积公式，

\[
 \prod_{j=1}^n(1-|x_j-x_{l,j}|/a)_+
 \ge\frac{\eta^2\theta M_0}{2H_{\rm cap}}.
\]

用 −log(1−u)≥u，得到

\[
 \|x-x_l\|_1\le aH_*,\qquad
 H_*=\log\frac{2H_{\rm cap}}{\eta^2\theta M_0}.
 \tag{5}
\]

每个 ℓ¹ 球体积是 (2aH_*)^n/n!，故

\[
 |F_i|\le K_\eta\frac{(2aH_*)^n}{n!},\qquad
 \int_Fq_{\rm dom}\le
 \frac{W K_\eta}{\theta M_0}\frac{(2H_*)^n}{n!}. \tag{6}
\]

第二步用所有非空 label 的数目≤W/(θM₀a^n)，并用 m_i/m≤1。对可数分配，(3) 本身使这些活跃 labels 只有有限个；Tonelli 无需 label 截断。

这保留了 density cap 的精确作用。仅用原 m>τa^n 时，可取 M₀=τ；半径就是 log(2H_cap/(η²θτ))。没有宣称 arbitrary L¹ source 有统一 H_cap，也没有把高峰分支当 bounded-density 输入。若 H_cap/τ 是 n 的固定幂，而 W/(τ|F|) 保持有界，则 (6) 中的半径仍 O(log n)，factorial 仍阻止常数覆盖。若 H_cap/τ 指数大，本文公式可能不再有效地排除覆盖。

## 3. 近均匀完整输入，连续窗口核心真赢家

取 a=1,b=2，n≥2，D=2n³，ε=n^(-2)，

\[
 \mathcal D=[-D/2,D/2]^n,
 \quad f(y)=\mathbf1_{\mathcal D}(y)
       \left(1-\frac{\varepsilon\|y\|_2^2}{nD^2}\right),
 \quad c=1-\frac1{4n^2},\quad\tau=1-\frac1n. \tag{7}
\]

这是完整非负 L¹ 输入，support 内 c≤f≤1，W=(1−ε/12)D^n，且 c>τ。定义核心 K=[−(D−b)/2,(D−b)/2]^n。对 x∈K，所有 L∈[1,2] 的方体都在来源盒内，准确平均是

\[
 h_L*f(x)=1-\frac{\varepsilon\|x\|_2^2}{nD^2}
               -\frac{\varepsilon L^2}{12D^2}. \tag{8}
\]

它对 L 严格递减。因此在整段连续窗口或任意包含1的有限子集，核心唯一真赢家为1；(8)≥c>τ，所以 K⊂E。全 E 则位于边长 D+b 的 support 膨胀盒。核心 m(x)≥c，故 (6) 可取 H_cap=1,M₀=c。

## 4. 全部来源均为高峰的有限尺度真赢家

该覆盖障碍不限于 near-uniform 原密度。令 h=1/2，并定义周期函数

\[
 p_n(t)=n\mathbf1_{\{\operatorname{dist}(t,h\mathbb Z)\le h/(2n)\}}.
\]

每个周期均值恰为1。用整份来源

\[
 f^{\rm peak}(y)=p_n(y_1)\mathbf1_{\mathcal D}(y)
       \left(1-\frac{\varepsilon\|y\|_2^2}{nD^2}\right). \tag{9}
\]

support 的正密度范围是 [nc,n]，而 nc>τsqrt(n)。因此按原阈值 τ 分出的 f1_{f>τsqrt(n)} 等于全部 f^{peak}（零密度点无质量）；没有加背景质量或删去低源。

取原有限尺度族 {1,2}。设 q=p_n−1，取 mean-zero 周期原函数 Q'=q，再取 mean-zero 周期原函数 R'=Q。两者存在，因为 q、Q 每周期积分零。为免与原 query 混淆，下式 Q、R 仅是这一维原函数。对 L=1,2，两端 x₁±L/2 均与 x₁ 同 phase；积分分部给

\[
 \frac1L\int_{x_1-L/2}^{x_1+L/2}t^2p_n(t)\,dt
 =x_1^2+\frac{L^2}{12}+2x_1Q(x_1)-2R(x_1). \tag{10}
\]

具体地，∫t²q=[t²Q]−2∫tQ；端点项是2x₁LQ(x₁)，而∫tQ=[tR]−∫R=LR(x₁)，因为端点 phase 相同且 ∫R=0。其他坐标是普通均匀平均，所以核心上

\[
 h_L*f^{\rm peak}(x)
 =1-\frac{\varepsilon}{nD^2}
  \left(\|x\|_2^2+\frac{nL^2}{12}
                     +2x_1Q(x_1)-2R(x_1)\right). \tag{11}
\]

除 L² 项外都与这两个候选尺度无关，故核心唯一真赢家仍为1。又因为曲率因子≥c、p_n 的方体平均=1，核心响应≥c>τ。完整来源质量满足 cD^n≤W_peak≤D^n（D 是周期的整数倍），并且 μ_peak≤n dy。因此 (6) 对这份全高源输入取 H_cap=n,M₀=c。

没有声称 (11) 对所有非整周期尺度成立，也没有把 finite winner 偷换为连续 winner。反驳“任意真实 finite winner 的普适覆盖合同”已有 {1,2} 足够；连续窗口反证由 (7) 单独提供。

## 5. 对全部固定 packets 的共同上包

两个输入都满足 |E|≥(D−b)^n，E⊂边长 D+b 的盒；q_dom≤1。令 H_cap=1 或 n，并设

\[
 H_n=\log\frac{2H_{\rm cap}}{\eta^2\theta c},\qquad
 U_n(H_{\rm cap})=
 \left(\frac D{D-b}\right)^n
 \frac{K_\eta(2H_n)^n}{\theta c\,n!}
 +\left(\frac{D+b}{D-b}\right)^n-1. \tag{12}
\]

式 (6) 与 boundary 体积给，对任意满足 (1) 的分配、其所有实际资格 (2)，

\[
 \frac1{|E|}\int_Eq_{\rm dom}(x)\,dx
 \le U_n(H_{\rm cap}). \tag{13}
\]

固定 η∈(0,1]、θ=n^(-1/2) 后，H_cap=1 的 H_n=(1/2)log n+O_η(1)，H_cap=n 的 H_n=(3/2)log n+O_η(1)。由 n!≥(n/e)^n，两种第一项均趋零；第二项是 O(n^(-2))，因为 D=bn³。故 (13) 的右端趋零。η 常数如何选择只改变常数，不改变失败结论。让 η 随 n 大幅变小会改变合同与 square fee，不能当成原常数 η 的解决方案。

两个输入的弱比仍有共同夹逼

\[
 \tau(1-b/D)^n
 \le\frac{\tau|E|}{W}
 \le\frac\tau c(1+b/D)^n\longrightarrow1. \tag{14}
\]

这也解释了为什么失败不否定互补路线。近均匀输入满足 f≤1，因而可以调用本地已核的 dimension-free cube L² 常数 C_B：

\[
 \tau|\{M_\square f>\tau\}|
 \le\frac{C_B^2}{\tau}\|f\|_2^2
 \le C_B^2 C W\quad\text{若 }f\le C\tau. \tag{15}
\]

准确本地常数口径见 [critical soft L² proof §1](../critical_soft_l2/proof.md) 的 C_B，以及 [whole-joint 独审](../full_soft_l2_extension/noise_parameter_review.md) §5；这里未另报未核数值常数，也未增加一般 L² 文献引用。条纹输入的 density cap/τ约 n，直接 (15) 只给 O(n)W；但 (14) 已解析证明该具体输入弱比很小。它说明仅剔除 f≤τsqrt(n) 的已付支，不能保证剩余高源出现 constant-η、sqrt(n)-dominance packets。高源也可以在 query 尺度上完全协同。

## 6. 可证的弱替代：fullness 可以有，份额大小仍未付

对任意原 f∈L¹、合法可测 R(x)∈[a,b]，取预固定嵌套 dyadic 来源网格，μ_i^{(j)}=μ 限制在第 j 层 cell。把完全位于 Q(x,R(x)) 的 cells 作为 fullness labels，定义其总捕获份额 q_full,j(x)。对每个固定 x，来源 μ 不给 query 面质量；每个严格内部来源 y 所在的 cell 最终包含于 query。嵌套性与 dominated convergence 给

\[
 q_{\rm full,j}(x)\uparrow1,
 \qquad\frac1{|E|}\int_Eq_{\rm full,j}\,dx\longrightarrow1. \tag{16}
\]

因此对有限 |E|，一份足够细的、全局预选的网格可有 η=1、平均 fullness≥1−δ；不是 receiver 自己划分 packets。该层可依输入和阈值选，但不因 x 而换。对于一般原子 source 和连续 adaptive radius，selected 面可能承载正来源质量，不能照搬此 L¹ 证明；finite radius family 的面例外则可由 Tonelli 排除。

式 (16) 没有 θ 下限。上述两输入正说明被完整捕获的来源可以分成大量份额远小于 n^(-1/2) 的合法 packets。质量树可把细层向上合并，但一份 packet 达到必要质量后，其完整捕获中心就受 (4)–(6) 限制；任何一次非重复来源分配都会遇到同一 obstruction。允许重用多层来源之后，须另付跨层/同-source 交叉，不能再应用 Σμ_i≤μ 的一次性费用。

这留下一个具体可证而未付的替代接口：保留 η-fullness，按真实 p_i=m_i/m 的份额 band 使用逐 band square，再证明同-source 跨 band 的共同预算；或者直接支付 cooperative source 的原弱交通。不得由 (16) 自动删掉所有小份额 incidence，亦不得把小份额叫几何 diffuse。

## 7. 数值核验与原 geom 距离

本 prefix 先登记后执行一次：[registration](packet_coverage_feasibility_20261007_registration.md)、[guard](packet_coverage_feasibility_20261007_guard.py)、[results](packet_coverage_feasibility_20261007_results.json)、[receipt](packet_coverage_feasibility_20261007_receipt.json)。三轮 n=64,256,1024，η=1/4；49 个严格 Fraction 比较 PASS。用 exp 的有理 Taylor 下界核对 H_n 上包，没有浮点比较作决定。两输入 U 的保守有理上包均小于4/n²；近似如下：

|n|近均匀上包|全高源条纹上包|
|---:|---:|---:|
|64|0.000488401|0.000860685|
|256|0.0000305181|0.0000305181|
|1024|0.00000190735|0.00000190735|

后三位相同来自第一项极小，显示值不表示它为零。守卫仅核新的解析参数、packing 和 boundary 上包；未做空间积分、packet 优化、actual gate 样本或 general endpoint 拟合。不存在新 paid fee 入主账。

本文严格否定一个足以闭合弱输出的附加覆盖合同，并给出有限 winner 的高源版本及 fullness-only 弱替代。原 geom 的 FIRST、fullfuture、CPGP、唯一 LCA、all-history 仍未在这一 hard source 证明中使用；要转回当前原 R_angle/R_capture,heavy 等余项，仍须把这些真实门给出的交通映射到合法 source packets 或另一已付 cooperative 支。没有从本文核表示或 (14) 推出一般 sqrt(n)n^{o(1)} 上界。
