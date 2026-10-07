# 原空间 projection overlap：预注册探路与周期／Rⁿ范围

2026-10-07，执行前落盘。只新增本前缀；不重复旧 weak 或熵测试，不改主账。当前目标是有限原参数网、原 c=1 symbol、完整饱和来源下的 weak-normalized Dθ，另记未截峰的 D 与 pair 上界。结果和审计将在同稿补入。

完整原 n 维周期输入由两个相位 θ₁=Σ_{i≤n/2}x_i、θ₂=Σ_{i>n/2}x_i（mod 2π）表示。每坐标 Gi 仍是原 g₁(ξ²)=log(1+ξ²)/ξ²，原 P_r/K_r 在此子类闭合为二维 Fourier 算子；没有换原跳长或使用 discrete Laplacian。两相位的 Haar 推前均匀，全部接收积分为 normalized Lebesgue phase average。完整 n-torus 体积因子 (2π)ⁿ 同时乘 W、D 和输出，normalized 比率相同。

三轮 n=8/32/128、phase N=128/256/512；每个来源固定 J=8/16/32 的 nested r_j=(j/J)² 网，阈值 α=1/2,1,2（κ=1）全部保留。每个 n 用单 packet 控制、正权四层同心 packet、偏移双 packet 各四层。下界总表实际读取 B41（先混合完整来源的共享分辨率）与 B62（整个来源共享连续同心尺度）；本轮只进入有限正权、多分辨率、嵌套及共同相位机制，不称原 B41 lattice 或连续 B62 完整下界资格。

用固定非负 U 定义原饱和输入 u=A U、Ω={U>0}、ν_b=(Su+κ)1Ω、μ_b=κ1Ω−Su1Ωc。包络 A 在看到输出前由解析原核界预定，不调 λ。原 w₁(t)=∫₁∞e^{−|t|v}v^{−2}dv≤e^{−|t|}，其 2π-periodization 的 supremum ≤(1+e^{−2π})/(1−e^{−2π})<2。若 packet 的一维 support 全宽为 δ、峰高 ≤1，每次原坐标卷积 GiU≤2Σ_lw_lδ_l。故 A=κ/(2nD_eff)、D_eff=Σ_lw_lδ_l，给全连续 periodic ΣGi u≤κ，进而 ν_b≥0、μ_b≤κ，μ_b 在 Ω 内饱和 κ。没有分层重新归一化或独立领取 W。

选 A={Ωc:M_grid>ακ}，M_grid=max_jK_{r_j}ν_b，A_j 用最小 argmax tie rule。主权 g_j=(ακ/M_grid)1_Aj，π_j=P_{r_j}g_j(X_j)。使用真实原 P 增量的 backward survival 与 future occupancy 递推，计算

\[
\alpha\kappa|A|=T_{\rm stop,\theta}+D_\theta,\quad
D_\theta\ge0,\quad
D_\theta\le\sum_{i<j}\mathbb E[\nu_b(X_0)\pi_i\pi_j].
\]

未截峰 a_j=1_Aj 仅是更强对照，不能以它的失败否定主 weak 合同。pair 与精确 D 分开，避免把 pair 放大丢掉的高阶取消误认作 D。每个空间／参数网都直接选择全部 Lebesgue receivers，绝不按 ν_b 抽来源再当 receiver。

FFT 循环网使用原 Fourier symbol 的有限采样；逆 FFT 的 Markov positivity、质量和 probability 屏幕全部保存，负值不截断。FFT 及非线性乘积 aliasing 没有 interval 认证。固定 n32 nested 输入 α=1,J32 另做 N128/256/512 稳定检查；同连续输入与阈值、无验证调参。各 r-grid 最大值不是连续原 winner，top-two response 和 threshold near-margin 单独计数。

周期到 Rⁿ：取 D 为 2π 的整数倍，令 u_D=u_per1_{[-D/2,D/2]^n}、Ω_D={u_D>0}，再构造 ν_D=(Su_D+κ)1ΩD、μ_D=κ1ΩD−Su_D1ΩDc。因为 0≤u_D≤u_per、ΣGi u_per≤κ，这份完整正 L¹ 原来源仍满足 μ_D≤κ、ν_D≥0 与 saturation，而不是直接声称 cutoff ν_per 仍是同一饱和源。原 G₁ 的一维一阶绝对矩为 2/3，故

\[
0\le\nu_D-\nu_{\rm per}1_{[-D/2,D/2]^n},\qquad
\|\nu_D-\nu_{\rm per}1_{[-D/2,D/2]^n}\|_1
\le\tfrac23n\|u_{\rm per}\|_\infty D^{n-1}.
\]

这里不把这条来源 L¹ 误差冒充 actual selector/D 的完整迁移误差：K/P kernel 的边界传播、输出阈值 margin、最大响应 tie margin及 Ω 外新增 receiver 均须另控。原 P 累积每坐标二阶矩 r/2、增量二阶矩 (s−r)/2，提供有限网边界 transport 控制的起点；本轮不会据此认证所有 actual Rⁿ winner 重构。只报告周期有限网的原 symbol 数值证据和已明示的 Rⁿ 来源误差。

## 1. 精确 D 与 pair 的实施公式

设 q_j=P_jg_j，V=V_{j+1,j} 为真实原 P 增量，j=1,…,J（r₀=0 的外域权为零）。使用正 future-union 递推

\[
c_J=q_J,\quad d_J=0,\qquad
c_j=q_j+(1-q_j)V c_{j+1},\quad
d_j=V d_{j+1}+q_jV c_{j+1}.
\]

Cstop=P₁c₁、Cdef=P₁d₁。另递推 h_J=q_J、h_j=q_j+Vh_{j+1}，得到 Ctotal=P₁h₁=Cstop+Cdef。pair 用 e_J=0、e_j=Ve_{j+1}+q_jVh_{j+1}，Cpair=P₁e₁。所有卷积使用 (2) 的原 c1 symbol；θ 主权和未截峰 a_j 同步计算。这保留 future-union 截断，精确 D 不由 pair 线性相减得到。

另核

\[
D_\theta=\sum_{j<J}\int(P_j\nu_b)q_j(Vc_{j+1})dx
=\int\mu_b C_{\rm def}+\int u\,S C_{\rm def}.
\]

结果分别保存正停止源、精确 D、pair 与 signed 势项，不将 ∫uSCdef 改成正能量强合同。每条记录有 source projection、stop+D、direct future-union、pair recursion、obstacle S 作用的残差。源/核/概率场不作数值 clipping；packet 正部是预定义输入本身，而不是输出修正。

## 2. 第一注册批：83 条全部空，不能支持合同

第一注册批 81 条主记录与 2 条固定输入空间稳定记录一次完成，34.09 秒。预设 α=.5、1、2 的非空超阈记录数**均为 0**。尤其 α=1、2 没有任何原空间核心压力。

|n|single packet 的 Ω 外 grid Mmax|nested 四层|shifted nested|
|---:|---:|---:|---:|
|8|0.0735030|0.165221|0.169822|
|32|0.0753761|0.158284|0.160919|
|128|0.0735918|0.156818|0.157186|

表列本轮主空间网 J32 的真实 grid 响应，不是 continuous supremum。83 条中的 Dθ、pairθ、Tstopθ、μCdef、uSCdef 和未截峰对照全为浮点 0，因为接受集合为空；**不能**据此说一般 D 预算小、signed 势项已付、恒等式有非零样本支持，或推出解析空集。

FFT 来源障碍身份最大残差 ≤2.23e−16，ν/μ 总质量误差 ≤6.94e−18。μ 最小值约 −3.33e−17，如实保留；其未修正。核 positivity 是浮点屏幕，非有向区间认证。n32 固定 nested α1/J32 的 N128→256→512 共同节点 M 差分别约 .00536、.00277，来源平均质量差约 .00125、.000616。D 的差为零只是空集合；winner 在空 union 上的零 disagreement 也不构成稳定性证据。

## 3. 空批后的合法幅度强化：独立登记，仍全部空

根据第一批无压力诊断，root 核准了一次解析输入强化。另保存 `projection_overlap_spatial_core_strengthened_registration_20261007.json`、独立脚本/结果/数组，不覆盖第一批。原 λ/κ、α、r 网、空间网、中心、宽度、正权全部保持；仅幅度按更紧的**连续原饱和约束**替换。这是已记录的自适应输入修正，不是独立验证，也不是事后降低 α。

常数证明：原 w≤e^{−|t|}，periodization supremum ≤cothπ。由 π>3、e>8/3，

\[
e^{2\pi}>e^6>(8/3)^6=262144/729>201,
\quad\coth\pi<\frac{202}{200}=C_0=101/100.
\]

packet 单轴 b(t)=(1−(2t/δ)²)³_+ 的积分恰为 16δ/35，故 Gb≤L=C₀16δ/35<1（全部 δ≤1）。对一个 product packet U=b₁b₂，

\[
S U=n b_1b_2-\tfrac n2[(Gb_1)b_2+b_1(Gb_2)],
\]

而 x,y∈[0,1] 时 L(x+y)/2−xy 的双线性角点最大值为 L/2。因此 −SU≤nL/2；完整正混合由线性求和得到

\[
-S U\le\frac{8C_0n}{35}D_{\rm eff},\qquad
A_{\rm strong}=\frac{35\kappa}{8C_0nD_{\rm eff}}.
\]

这同时保证周期 ν_b≥0 和 μ_b≤κ，比原保守幅度强 3500/404≈8.663 倍，不依赖已选 receiver 输出。

本轮多分辨率首先在同一非负 potential U 中正混合，再一次共同构造 ν_b，未把 signed S 分量当作独立正源领取 W。若需检查来源正分解，可给每个 packet-layer 权 w_l 分配 κ_l=κw_lδ_l/D_eff，并以共同 Ω 定义 ν_l=(Aw_lSU_l+κ_l)1Ω。上述逐层下界保证 ν_l≥0，且 Σν_l=ν_b、Σκ_l=κ。它不是各层自行重选 Ω、输出或阈值的 B41 完整原族，只准确记录本轮进入的同源正混合机制。

强化批仍一次完成三轮及固定空间稳定组，共 83 条，34.22 秒。预设 α=.5、1、2 的非空记录数仍**全部为 0**。不继续微调直到过阈。

|n|强化 single 的 Ω 外 grid Mmax|强化 nested|强化 shifted nested|
|---:|---:|---:|---:|
|8|0.308431|0.251558|0.210906|
|32|0.299242|0.249513|0.205315|
|128|0.303498|0.250174|0.204347|

因仍无接受集合，所有主 D、pair、停止源及 signed 势项都是空集的数值零。它没有反驳一般合同，也没有提供支持。单 packet 的 μ 在 Ω 外最大值约 .643–.663κ，nested 约 .304–.328κ，shifted nested 约 .152–.164κ；全局 max μ=κ 来自 Ω 内饱和。即使合法强化，所选来源的外域交通仍不足预设阈值。

强化 fixed n32 nested 的 N128→256→512 M 差约 .01562、.00270；来源平均质量差约 .00115、.000576。FFT/source 离散误差尚未认证，网差不能变成 continuum 或 Rⁿ winner 误差上界。本批正 D 递推的非空分支没有被这些空记录实际压力到，相关代数另有有限滤过守卫，但本批不是非零原模型证书。

## 4. 强化周期来源的合法 Rⁿ 迁移范围

强化幅度只证 Su≥−κ、Ω 外 ΣGi u≤κ，未证全空间 ΣGi u≤κ。所以直接硬截 u 的第一批证明不能原样搬到强化输入：截盒外的新增 Ωᶜ 接收点可能违反 μ≤κ。

已独立审阅 root 新稿 [periodic_saturated_slow_cutoff](periodic_saturated_slow_cutoff_20261007.md)，其缓变截断修复成立。设 Q_R=[−R,R]ⁿ、w_R=√R，χ_R 为 Q_R 上 1、Q_(R+w_R) 外 0 的 coordinate-Lipschitz cutoff。原 G₁ 一阶绝对矩实际为 2/3（使用粗 m₁≤1 也足够）。有

\[
\|[S,\chi_R]u\|_\infty
\le n m_1\|u\|_\infty/\sqrt R=\epsilon_R,
\quad S(\chi_Ru)\ge-\kappa-\epsilon_R.
\]

取 κ_R=κ+ε_R 并**重新构造** ν_R、μ_R，得到完整正 L¹ 饱和来源、μ_R≤κ_R；来源平均质量趋周期 W，κ_R→κ。其 commutator L¹ 界、边界环 o(Rⁿ) 以及 ν_R−χ_Rν 的身份方向均正确。内盒的来源差 ≤2ε_R，结合统一来源峰界及 K 的轴一阶矩 union-tail，得有限网响应的一致趋同；严格 weak 水平集通过阈值余量/Fatou 迁移。

该引理使有界连续周期饱和模型成为一般 Rⁿ L¹ weak 命题的合法必要探针，固定 n 先取 R 极限；没有声称 n/R 的统一复杂度。它不认证当前 FFT 输入/积分，不自动迁移 argmax 并列标签或精确 Dθ；后者仍需阈值/并列集和 source-weighted 积分误差控制。本轮没有据这个解析迁移把空数值记录写成完整 Rⁿ 原模型结论。

## 5. 结论和复现

两批 166 条记录全部保留。结果 JSON 含注册/script hash 与每个 NPZ hash；源数组、原 Gi/P/V symbol、r 网、全部 receiver 标签及 θ 权、Cstop/Cdef/Cpair/survival、signed S Cdef 列都保存，足以重构本次浮点递推。无 source Monte Carlo 或置信区间声明；原 source 和 node variants 全程耦合于同一输入。

本轮实际证据是：这些注册的饱和来源在所测原符号周期离散网下没有非空 weak projection core。一般 source-weighted Dθ、future-union 交通和 signed 势项预算仍未解决；没有新增已付项、没有一般弱界、没有 asymptotic 拟合或 actual geom 反例。
