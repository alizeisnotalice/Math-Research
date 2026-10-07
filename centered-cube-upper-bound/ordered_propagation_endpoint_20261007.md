# 原有序传播的完整入口端点：正吸收源简化与强范数范围

2026-10-07。只写本轮新文件，不改总稿或主账。读取 checkpoint 最新收口、上一轮 ordered_noise_frozen_transfer、原 tex 4878–5000 与 EI/FJL 来源段；使用此前已读 J01/J05 的真实生成元和一次来源合同工作流。技能不是端点定理。本文物理尺度固定为 ell；原 early 来源入口 t≤1/2，而 future receiver 时间可到 1。

本轮结论是准确缩小转移接口：**同一 killed 演化的实际正吸收源只需原 N_s 的固定初始密度弱端点。** 一般时空源可用 root 新 hazard 锚点引理有偿转移，费用为 log n。固定冻结 H 的 O(log n) 端点不能代入原完整 Bernoulli family。本文没有证明该 family 的 polylog 弱端点，没有新增 cube 费用。

## 1. 先确定输出和来源合同

令 B(v)=v/log(1+v)−1，G_c 的 Fourier 乘子为 1/(1+cB(ell²ξ²))。原转移为
\[
 K_{s,t}=\bigotimes_{i=1}^n[(1-r)I+rG_{1-t}^{(i)}],
 \qquad r=\frac{s-t}{1-t},\quad 0\le t\le s\le1.
\]
保持原有序过程、holding 和全部坐标面，不替换成 reset 链。

对一份固定正 L1 时空源 F，支持 t∈[0,1/2]，定义
\[
 {\cal T}F(x)=\sup_{0\le s\le1}
  \int_0^{\min(s,1/2)}K_{s,t}F_t(x)\,dt,\qquad
 M=\int_0^{1/2}\|F_t\|_1\,dt.                    \tag{1}
\]
可允许 L1 密度时间脉冲，指定端点归属；本文反例源没有时间脉冲。sup 先取有理时间和端点，再用原正有限掩码展开及 L1 路径的共同版本。一般 signed forcing 以 |F| 支配，M 是总变差；该操作不保留它的符号取消或来源父组语义。

原实际来源合同更强：固定同一 ν、D、ell，以及原生成元的 killed live 演化 η，令
\[
 \rho_t=V_t\eta_t+\mathsf E_{t,\ell}\eta_t\ge0,\qquad
 \partial_t\eta_t+\mathsf A_{t,\ell}\eta_t=-\rho_t.
\]
原 NB--M 已证明
\[
 N_{s,\ell}\nu=\eta_s+\int_0^sK_{s,t}\rho_t\,dt,\qquad
 \|\eta_s\|_1+\int_0^s\|\rho_t\|_1dt=\|\nu\|_1.   \tag{2}
\]
因此任何正子源 0≤F_t≤rho_t，包括 early/首跳/历史接受的正限制，都有
\[
 {\cal T}F\le \sup_{0\le s\le1}N_{s,\ell}\nu.     \tag{3}
\]
这一步不要求 receiver 与 source 独立，不把 rho 当普通输入另领一份 W，也不对条件弱常数积分。它是 **旧 Duhamel 正恒等式的应用**，不是新增空间定理。不同 ell 的 killed 演化会改变 rho 和 ν；不能把这些来源并成一份 M。

所以该实际合同所需的最小噪声端点是
\[
 \left\|\sup_{0\le s\le1}N_{s,\ell}\nu\right\|_{1,\infty}
 \le A_0(n)\|\nu\|_1                             \tag{4}
\]
对同一任意非负 L1 密度成立。这里入口固定 t=0、c=1。一般时空源的锚点版本则需要每个固定入口 c∈[1/2,1] 的完整 future
\[
 \|\sup_{0\le r\le1}[(1-r)I+rG_c]^{\otimes n}f\|_{1,\infty}
 \le A_{\rm ord}(n)\|f\|_1.                     \tag{5}
\]
不能用上一轮只覆盖 r≤1/(2ceil(sqrt n)) 的正支配来证明 (4) 或 (5)。

## 2. Root 的连续源转移已消除非法条件弱平均

只读审核 hazard_anchor_source_transfer_20261007.md §§1–4。其有限非负弱求和引理为：若 ||f_j||_(1,infinity)≤a_j、A=sum a_j，则
\[
 \|\sum_{j=1}^Nf_j\|_{1,\infty}
 \le[3+2\log(2N)]A.                             \tag{6}
\]
证明先付 {max f_j>lambda}，在补集按 lambda a_j/(2A) 截断；低部分总和≤lambda/2，高部分的截断层饼积分给质量熵费用。没有无截断 L1 积分。

对正 cocycle U，若每个来源格左端 a_j 满足 U(t,a_j)≥theta I，则
\[
 {\cal T}_U F\le\theta^{-1}
  \sum_j\sup_{s\ge a_j}U(s,a_j)\nu_j,\qquad
 \nu_j=\int_{I_j}F_tdt,\quad\sum_j\|\nu_j\|_1=M. \tag{7}
\]
正确方向来自 U(s,a_j)=U(s,t)U(t,a_j)≥theta U(s,t)。对 t>s 尚未输入的部分只在右边加入正上界。分的是一次 F，不是重复存活的 eta。

原 early [0,1/2] 取 2n 格、宽 1/(4n)，完整 holding
[((1-t)/(1-a_j))]^n≥(1−1/(2n))^n≥1/2。结合 (5) 得
\[
 \|{\cal T}F\|_{1,\infty}
 \le2[3+2\log(4n)]A_{\rm ord}(n)M.              \tag{8}
\]
所有 future 终点都保留，没有免费 sup_c；每个 anchor 的输入质量相加只为 M。原实际正 rho 可用更好的 (3)，不用 (8) 的 log n。

Root 同一证明对单个固定冻结 H 给 O(log²(n+2))M；共同来源再结合物理尺度网给 O(sqrt n log²(n+2))M。若 M 本身为 O(sqrt n)W，后者是 O(n log² n)W。它不是原时变传播的费用。分段冻结 U 亦不是一个固定 H；原 signed forcing 先前仅有时空总变差 O(sqrt n)W，不能直接宣布输出误差已付。

## 3. 完整入口的 L2 守卫及投影路线限度

这是输入层面的解析事实，不是弱 L1 证明。设 T_r=[(1-r)I+rG_c]^{tensor n}，Fourier 中 theta_i=1−Ghat_c(ξ_i)∈[0,1]，
\[
 m_r(\xi)=\prod_i(1-r\theta_i),\qquad
 -r\partial_rm_r=\sum_i r\theta_i\prod_{j\ne i}(1-r\theta_j)\le1. \tag{9}
\]
最后一式是参数 r theta_i 的独立 Bernoulli 先验中“恰一成功”的概率恒等式，仅用于固定乘子代数，**不是 actual 后验独立性**。m_r 单调非负，所以
\[
 \int_0^1r|\partial_rm_r|^2dr
 \le\int_0^1(-\partial_rm_r)dr\le1,\qquad
 \int_0^1r\|\partial_rT_rf\|_2^2dr\le\|f\|_2^2. \tag{10}
\]
对 n≥2，r≤1/n 的完整正掩码包络质量≤(1+1/n)^n≤e。对 r≥1/n 用 Cauchy–Schwarz 的 log n 积分。因此
\[
 \|\sup_{0\le r\le1}|T_rf|\|_2
 \le(e+1+\sqrt{\log n})\|f\|_2.                \tag{11}
\]
这里有限多项式展开保证所需逐点全 r 版本。能量守卫没有支付 bad L1 来源。尤其 proper faces 的核对一般空间平移是互相奇异的，不能直接将 (11) 接标准全维 L1 kernel translation Hörmander 条件。

另一个准确但未闭合的结构：Q_c 的乘子 (1+cB)^(-1/2) 是 Gamma(shape 1/2) 从属 Brownian 概率核，G_c=Q_c²。对称 Markov 二步 dilation 及 product extension 给两个 conditional projections E0、E_i，E0 E_i E0=G_c^(i)。令
\[
 S_\tau=\prod_i[e^{-\tau}I+(1-e^{-\tau})E_i],
 \quad r=1-e^{-\tau};
 \qquad T_r=E_0S_\tau E_0.                     \tag{12}
\]
中间确是 commuting-projection semigroup，但只能得到 T_*f≤E0(S_*E0f)。**外层 conditional expectation 不保 weak L1**，不能把中间的弱端点直接移出。一般有限概率 fiber 上取 Z=2^k，概率 2^(-k)，k=1..m，其余概率取零，则 Z 的弱范数≤2、E Z=m；这只是指出该证明步骤不合法，不是原 kernel 或 actual 输入反例。本文没有把 (12) 登记为新端点。

## 4. 原 first-event 退出源：强范数不能替代弱端点

以下只检验新泛函的强版本是否过强；不否定 (4)、(5)、(8)，不声称全 geom 历史合格。固定 ell=1、W=1、n≥512，
\[
 C=[-1/(2n),1/(2n)]^n,\quad
 \nu_\epsilon=\text{uniform density on }[-\epsilon,\epsilon]^n,\quad
 \epsilon=1/(100n^2).
\]
该 C 为原出生小格的 1/n 尺度。T_flat=(12+floor(log2(n)/2))/n；其≤原早首跳截止 T_n=(12+log(n+2))/n，因为 log 2>1/2，且 n≥512 时 T_flat<1/2。

使用 **原过程首次 event 即退出 C** 的正子源，此前没有内部跳：
\[
 d\kappa_0(t,w)=
 {\bf1}_{0<t<T_{\rm flat}}\sum_i(1-t)^{n-1}
 \int\nu_\epsilon(dy)\,
 {\bf1}_{w_i\notin C_i}G_{1-t}(w_i-y_i)\,dw_i\,
 \delta_{y_{-i}}(dw_{-i})\,dt.                 \tag{13}
\]
对该 L1 初始密度它也是 L1 时空密度；不是一个人为 reset 源。它是 D=C、无额外杀势的原 exit rho 的正子源。完整 first-event 源 kappa_all 在 t∈(0,1) 的总质量为一，且精确分解
\[
 N_s\nu_\epsilon=(1-s)^n\nu_\epsilon+
             \int_0^sK_{s,t}\,d\kappa_{\rm all}(t). \tag{14}
\]

原 G_c 的正 Exp 热钟表示给 G_c≤c^(-1)G_1=c^(-1)w；t≤1/2 时 G_(1-t)≤2w。原 w(0)=1。因此 first event 的内部落点概率≤2/n；T_flat 后尚未 first-event 的质量为 (1−T_flat)^n。被 (13) 遗漏的完整 first-event 源质量至多
\[
 b_n=(1-T_{\rm flat})^n+2/n.                   \tag{15}
\]
这是 source 质量，不是把完整捕获质量替换成 rare capture。

旧 EI/FJL 使用的完整系数峰值包络为
\[
 D_n=\sum_{j=0}^n{n\choose j}(j/n)^j(1-j/n)^{n-j},
 \qquad0^0=1.
\]
对任意固定正时空源（甚至 t<1 全区间），正 Tonelli 给
\[
 \int{\cal T}F\le D_n\int\|F_t\|_1dt.          \tag{16}
\]
每个 t 的 G_(1-t) 形状有质量一，系数峰值总质量 D_n；这是 **旧正 strong 包络的合同应用**，不计为本轮新成果。可据此把遗漏源的 maximal 积分支付为 D_n b_n。

构造互不相交 receiver channels，非空 A 的 unvisited 坐标满足 |x_i|≤epsilon，visited 坐标满足 1/(2n)+epsilon<|x_i|≤n+epsilon。在该 channel 选 s=|A|/n。原 w 的真实位移满足
\[
 \Pr\{1/(2n)+2\epsilon<|Y|\le n\}
 \ge1-1/n-4\epsilon-e^{-n}\ge1-2/n.            \tag{17}
\]
这里 w 的近零概率≤区间长度，因为 w≤1；尾概率为 2∫_1^infinity e^(-nu)u^(-3)du≤e^(-n)<2^(-n)。同一初始小盒的全部 y 都满足上述 channel 余量。各 mask 的正贡献积分至少为其 peak weight 乘 (1−2/n)^|A|，可统一降至 (1−2/n)^n。空 mask 在 C 内，故不参与外域。

在这些不交 channels 上用 (14)，holding 初始项为零；再扣掉遗漏源的 (16) 得
\[
 \|{\cal T}\kappa_0\|_1
 \ge(1-2/n)^n(D_n-1)-D_n[(1-T_{\rm flat})^n+2/n]
 \ge(D_n-1)/32\ge\sqrt n/128.                 \tag{18}
\]
常数的全 n 解析核对：n≥4 时 (1−2/n)^n≥1/16（对 n 求导其 log 单调）；n≥512 时 T_flat≥16/n，late≤e^-16<2^-16，D_n/(D_n−1)≤2，且 2b_n≤2^-15+4/n<1/32。至于 D_n−1≥sqrt n/4：j∈{1,..,n−1} 时 Bin(n,j/n) 以 j 为 mode，variance≤n/4；Chebyshev 给半径 sqrt n 区间质量≥3/4、整数数≤2sqrt n+1，故该 mode 概率≥1/(4sqrt n)。求和并加 j=n 的一项即得。无需 Stirling 或维数拟合。

因此即便来源来自原 first-event exit 合同，也不能要求 (1) 的 strong L1 费为 polylog。它不构成 weak L1 下界；channels 的密度高度没有被调成共同阈值。也没有核验原 hardwinner/hardband、soft FIRST、CPGP、SC-F/R、nonconc 或全部历史门，不能叫 actual geom 样本。

## 5. 本轮守卫及终态

先保存同前缀 registration，再运行新 guard 一次。n=512/1024/4096，epsilon、cutoff 与 channel 固定如上。直接 exact integer/Fraction 计算 D_n、late 和 (18)；巨大分子以 binary SHA256 和 bit-length 保存，浮点仅供显示。不是重跑旧 D_n 收据，不采样 G，不以有限实验代替 (13)–(18) 的解析核证明。

36 项 PASS_EXACT，三轮 D_n 显示为 29.03048656、40.77595410、80.88039616；(18) 右边未再粗化的 lower 显示为 3.66528291、5.29293074、10.76586341。源质量≤1，所以这些是原 first-event 子源的强范数范围守卫；不是 weak norm 测量。脚本和 registration hash 在结果 JSON 中。

最终未付接口为 (4) 或 (5) 的 polylog 弱界，以及从它到 moving h_L/真实 receiver 所选尺度的共同来源转移。Root 的 anchor 引理已合法消除一般连续来源的条件弱平均障碍；actual 正 rho 用 (3) 更省。本文的 L2 与投影结构未提供 bad-source 空间支付，强范数限制也没有否定弱路线。原 cube 主账、R_dagger 与一般目标状态均不变。
