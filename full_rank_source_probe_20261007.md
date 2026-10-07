# Full-rank 原饱和来源与 moving-scale dilation 的最小接口

日期：2026-10-07。模型：gpt-6.1-sol high。使用 L02，并阅读其 method、provenance、cube-interface。本文不改旧两相位空批、不调其阈值、不重复 fixed K 的已知预算。无空间 Monte Carlo 运行；新三轮只核原符号的有理区间和完整张量源的代数可行性。

## 1. 本轮能交付什么

完整 n 个坐标的正嵌套张量 potential 可以建立合法的原 c=1 饱和源，幅度不需要两相位模型中的 1/n 缩小。它有可复算的任意 receiver 响应和完整来源 rejection 表示，见 §2–3。

固定原 soft 核、改变 hard receiver sensor 时，accepted-count 恒等式仍成立，见 §4。但原 soft 物理尺度也随 receiver 变化时，固定 r 的原 P 核在任意两个不同尺度之间都**没有概率卷积商**，见 §5。这准确排除直接沿尺度串联原 P 商的路线；不排除斜向 r,L 链、多参数 dilation 或别的原共同参考。

fixed K_r 已有 O(sqrt(n)) 强 L1 正 Bernstein 包络，因此只对这个族压 D_theta 不会产生新的 sqrt(n) 一般付款，见 §6。根据当前任务边界，仅注册可执行工具，停止重型 MC。这里没有新的 actual geom 上下界。

## 2. 完整 full-rank 正来源：无 n 因子的饱和幅度

使用固定周期 2pi 的每坐标原 G=G_1。其原一维 Fourier 符号为

\[
 B(v)=\frac{v}{\log(1+v)}-1,
 \qquad \widehat G(\xi)=\frac1{1+B(\xi^2)}
 =\frac{\log(1+\xi^2)}{\xi^2},
 \quad\widehat G(0)=1.
\]

原 R 密度

\[
 w(x)=\int_1^\infty e^{-s|x|}s^{-2}\,ds\le e^{-|x|}
\]

的质量为 1；Fourier 积分直接给上式符号。周期化密度 w_per 至多 coth(pi)<101/100=:C0：在一个周期内指数周期和的最大值在 0，等于 1+2sum_{k>=1}e^{-2pi k}；pi>3、e>8/3 给 e^{2pi}>(8/3)^6>201，故 coth(pi)<202/200。这不是 FFT 核峰值估计。

取 0<delta_l<=delta_0<=1，

\[
 b_l(x)=\bigl(1-(2x/\delta_l)^2\bigr)_+^3,
 \quad I_l=\int b_l=16\delta_l/35,
 \quad 0\le G b_l\le L_l:=C0 I_l\le1.
\]

设 U_l(x)=prod_{i=1}^n b_l(x_i)，w_l>0、sum w_l=1，u=A sum_l w_l U_l，Omega={u>0}=(-delta_0/2,delta_0/2)^n（含正权 outer 层；边界是零集）。S=sum_i(I-G_i)。则

\[
 -S U_l=\sum_i(Gb_l)(x_i)\prod_{j\ne i}b_l(x_j)-n\prod_i b_l(x_i)
 \le L_l\sum_i\prod_{j\ne i}z_j-n\prod_i z_i\le L_l. \tag{1}
\]

最后一个多重仿射函数在 [0,1]^n 的顶点取最大值。全 1 顶点为 n(L_l-1)<=0；恰一个 0 为 L_l；至少两个 0 为 0。因此对完整 n 坐标而非两相位 reduction，

\[
 A=\frac{\kappa}{\sum_l w_l L_l},\quad
 S u\ge-\kappa,\quad
 \nu_b=(\kappa+S u)1_\Omega\ge0,\quad
 \mu_b=\kappa1_\Omega-Su1_{\Omega^c}\in[0,\kappa],
 \qquad Su=\nu_b-\mu_b. \tag{2}
\]

Omega 外 u=0，Su=-sum_i G_i u<=0，故 mu_b 上界也由 (1) 保证。来源为整份 nu_b；不把 signed tensor 项作为独立正来源，也不对各层重新归一化。原任务提及 B41/B62 的正多分辨率混合提供构造动机；此输入只实现正嵌套层机制，不宣称完整下界资格。

登记的 widths=(1/2,1/4,1/8)、weights=(1,2,4)/7、kappa=1。A=6125/606。该幅度是预先解析确定，不由超阈事件反推。维数大的 smooth tensor potential 的积分为 I_l^n，而 plateau 质量为 delta_0^n；因此固定 widths 下 plateau 可能支配质量。这个可行性问题也必须记录，不能看到空事件后隐去该源。

## 3. 原 K_r 的任意 receiver 精确公式与采源

写 c0=1_{(-delta0/2,delta0/2)}、d_l=c0(Gb_l)。因为 b_l 支撑在 outer interval，完整来源恰为

\[
 \nu_b=\kappa c0^{\otimes n}
 +A\sum_l w_l\left[n b_l^{\otimes n}
 -\sum_i d_l(x_i)\prod_{j\ne i}b_l(x_j)\right]. \tag{3}
\]

定义 T_r=(1-r)I+rG，a0,r=T_r c0、al,r=T_r b_l、dl,r=T_r d_l。K_r=T_r^{tensor n}，故任意 receiver x 的响应为

\[
 K_r\nu_b(x)=\kappa\prod_i a0,r(x_i)
 +A\sum_l w_l\left[n\prod_i al,r(x_i)
 -\sum_i dl,r(x_i)\prod_{j\ne i}al,r(x_j)\right]. \tag{4}
\]

这不是高维相位 reduction。预计算原一维卷积表后，每个 receiver 每个 r 的运算是 O(n * layer count)。脚本保留 prefix/suffix product 实现，避免零分母；它是浮点公式工具，未运行 continuum 表或 FFT 空间 pressure。正负项消去必须另记误差，不能 clipping 然后继续声称原 Su 或精确阈值。

记 J_l=integral c0 Gb_l，则周期 unnormalized Lebesgue 来源质量

\[
 W=\kappa\delta_0^n+A n\sum_l w_l(I_l-J_l)I_l^{n-1}>0. \tag{5}
\]

合法 rejection proposal 是 q=\kappa1_Omega+n u>=nu_b，质量 Q=\kappa delta0^n+n A sum_l w_l I_l^n。选择 plateau 或正 tensor component 再按 nu_b/q 接受，得到 X0~nu_b/W。来源拒绝率 W/Q 可能差，必须计入成本；不存在免费截去来源的步骤。

周期整体来源与 Rn 原饱和源之间，可用已独审的 [periodic_saturated_slow_cutoff](periodic_saturated_slow_cutoff_20261007.md) 的宽度 sqrt(R) 慢 cutoff，固定 n 下来源质量/体积和有限网严格 weak 超水平集可迁移。该引理不认证 FFT、grid winner 稳定性、accepted-count D 或其高阶 cancellation 的迁移；本轮也未进行这些迁移数值。

## 4. accepted-count 表示，以及仅改变 hard sensor 的合法扩展

固定原物理 soft G、有限递增 r_j 网，P_j=K_rj^{1/2}，V_{j+1,j}=P_{j+1}/P_j 是原 Markov 增量，见 [original_ordered_martingale_dilation](original_ordered_martingale_dilation_20261007.md)。取完整來源 X0~nu_b/W，第一步使用 P_1，再以这些 V 增量生成 X_j。令 v_j=K_rj nu_b，M=max_j v_j，E={Omega^c:M>tau}，Aj 是 E 上按最小索引打破 ties 的赢家分区。g_j=(tau/M)1_Aj，sum_j g_j<=1。

给定 path，独立取 Y_j~P_j(X_j,dy)，再以 g_j(Y_j) 接受，N 为接受次数。条件概率 q_j=P_jg_j(X_j) 在 [0,1]；自伴性和 P_j^2=K_j 给

\[
 E N=\frac1W\sum_j\int g_j K_j\nu_b
 =\frac{\tau|E|}{W},\qquad
 P(N>0)=\frac{T_{stop}}W,
 \quad E(N-1)_+=\frac{D_\theta}W,
 \quad E\binom N2=\frac{pair_\theta}W. \tag{6}
\]

它通过自伴核积分换测度，所以来源概率并未冒充 receiver Lebesgue 体积。N 与连续的 q_j 条件独立 coin 模型有同一 law。完整未来 union 给 exact D，而 pair 只是正二阶上界；不能以 pair 膨胀否定 D 预算。

**仅固定 soft、移动 hard 是可表示的。** 另取正、Lebesgue 保质量 receiver sensor H_j，并统一约定 H_j 作用于测试函数：H_j g(x) 是从 x 经该 sensor 采 receiver 后的 g 条件平均，其伴随 H_j^* 作用于来源密度。从 X_j 再经 P_j 以及 H_j 独立采 receiver，q_j=P_j H_j g_j(X_j)，对应 receiver 响应 F_j=H_j^*K_jnu_b。以 M=max_j F_j 定义赢家分区和 g_j，式 (6) 中的 K_jnu_b 换成 F_j 后仍成立。原 even cube sensor h_Lj^{tensor n} 自伴，H_j^*=H_j。相同 r 的不同 H 可以使用 identity path 增量；不必对 hard kernel 开平方。这一表示本身不提供 D 的 source-once 上界。

该 fixed-soft/variable-hard 家族没有原 receiver-dependent G_{c,L_s(x)}，也没有原完整 FIRST/history、source pair、CP/GP 与共同 fullfuture 资格。实际公式的 physical L 同时改变首跳 G、hard seed、continuation，见 [hard_average_transfer_endpoint](hard_average_transfer_endpoint_20261007.md) §1、§5。因此式 (6) 的 sensor 扩展不能直接替换 actual moving hard/soft 交通。其通用 accepted-count 补充由另一个 agent 的 acceptance-count 稿负责；本文只登记适用边界。

## 5. 同 r 跨物理 soft 尺度的原 Markov 商不可能

固定 c>0、0<r<1，每坐标物理尺度 ell 的原 P 符号为

\[
 p_r^{ell}(\xi)=\sqrt{\frac{1+(1-r)cB(ell^2\xi^2)}{1+cB(ell^2\xi^2)}}. \tag{7}
\]

B(v)=integral_0^1[(1+v)^t-1]dt 严格递增且无界。对 ell2>ell1 和 xi!=0，有 p2<p1，故非恒定的原商 R(xi)=p2/p1 严格介于 0 与 1。同时 |xi|->infty 时 p1,p2 都趋 sqrt(1-r)，所以 R(xi)->1。

若 R 是某概率卷积 rho 的特征函数，则

\[
 \frac1{2T}\int_{-T}^T R(\xi)d\xi
 =\int\frac{\sin(Tx)}{Tx}\,d\rho(x)\longrightarrow\rho\{0\}. \tag{8}
\]

Fubini 和绝对值至多 1 的 dominated convergence 证明此式；左侧由 R(xi)->1 趋 1。因此 rho{0}=1，rho=delta0，R 必恒等 1，与严格不等式矛盾。反方向商 p1/p2 在非零 xi 处大于 1，直接违反概率特征函数模长至多 1。故两个方向都不能是 Markov convolution。n 维只需沿单坐标的频率轴限制，结论相同。

这里排除的是保持 r 不变、直接用**原 P 商**跨 L 建共同链的方案。结论不涉及饱和源特选，也无需空间数值；不是 actual geom 或 weak 反例。r=0 的 identity 及 r=1 的零高频极限不由上述证明处理，本文不擅自扩展端点。

周期 Fourier 也不能解除这个商障碍：若周期商的整数频率系数趋 1，Cesaro/Dirichlet average 对任何非零 torus 点趋 0，迫使 rho 在 identity 的质量为 1，同样矛盾。有限 FFT 可能掩盖高频 atom 问题，不可据有限核的近似正性声称原无限频谱 Markov 性。

## 6. 为什么本轮不继续昂贵 fixed K pressure

固定原 Gi 时，K_rnu=sum_{A subset[n]}r^{|A|}(1-r)^{n-|A|}G_A nu。正性直接给

\[
 \int\sup_r K_r\nu\le D_n W,\quad
 D_n=\sum_{k=0}^n\binom nk(k/n)^k(1-k/n)^{n-k}=O(\sqrt n). \tag{9}
\]

端点按 0^0=1；中间项标准 factorial 比较为 O(sqrt(n/(k(n-k))))，求和为 O(sqrt n)。这是既有 fixed-family 正 Bernstein 预算，见 [late_multiface_source_budget](late_multiface_source_budget_20261007.md) 末节；本文不登记为新发现。并且 D_theta<=tau|E|<=integral sup K_rnu，所以新 full-rank 原空间 count 即使非空，也不能突破 fixed K 合同已有的阶数。

登记的可选三轮 n=4/8/16、J=8/16/32、S=256/512/1024、三个正层，receiver winner 全重算的成本约 O(S J^2 n * 3)，每阈值分别约 0.20/3.15/50.33 百万一维表操作，另加来源 rejection。仅在获得新的 actual moving-scale 接口或明确需要 fixed-soft/variable-hard 必要诊断后才激活，当前没有启动。固定阈值 alpha=.5/1/2、tau=alpha*kappa；未来任何空结果保留，调参必须另登记。

若未来用离散 original-symbol FFT 采 P、V，每个核的负质量/守恒与 aliasing 需记录，不允许 clipping；Hoeffding 只覆盖冻结的有限周期模型的 count，不覆盖 continuum 或 Rn。grid max 不等于连续 r 精确 winner。当前 deterministic guard 无随机种子、无 count 置信区间或事件记录。

## 7. 后续精确欠项，而非新 easy 分支

原 P 的每坐标有限 Levy measure 可直接由 (7) 写成

\[
 \Lambda_{r,ell}=\frac12\int_{1-r}^1G_{ca,ell}\frac{da}{a},
 \quad |\Lambda_{r,ell}|=-\frac12\log(1-r),
 \quad P_{r,ell}=\exp\{\Lambda_{r,ell}-|\Lambda_{r,ell}|I\}. \tag{10}
\]

保持 r 的缩放有相同总 Levy 质量而不同 jump measure；§5 给出概率商层面的严格失败。允许 r 随 ell 上升时，总 Levy 质量增加，可能补偿某区域的 density 亏损。一个**待验的充分条件**是每一步原 measure 差 Lambda_{s,ell2}-Lambda_{r,ell1}>=0，这才产生原 compound-Poisson Markov 商；Fourier 点态序不能代替正 measure 序，也不声称此充分条件必要。斜向链是否覆盖 actual 同一 fullfuture 的全部 receiver 参数、或能否构造多参数 dilation，仍未证。本轮没有把 r 自动 retune 成覆盖资格，更没有支付完整 FIRST 后的投影欠额。

## 8. 三轮小守卫与可复算状态

registration 在执行前落盘。新脚本只用 Fraction：log(y) 先精确二进制 range-reduce 到 [1,2)，再用 2 sum z^{2k+1}/(2k+1)、z=(y-1)/(y+1)，余项上界 2z^{2K+1}/((2K+1)(1-z^2))。所有值为正，log 上下界经 B=v/log(1+v)-1、p²=1-rB/(1+B) 时按单调方向反转，再将两个正 p² 相除。未调用浮点 log/sqrt；显示小数不参与断言。

三轮分别 n4/16/64，r、frequency 与 40/56/72 项登记固定；沿一个坐标轴是这些 n 维原张量的真实 marginal，不能把它解释成三轮高维 source pressure。已执行结果为 18 个原符号有理区间和 261 个多仿射顶点分类全部通过，运行约 0.0111 秒。前者严格核 squared quotient<1、reverse>1；后者对 3 个 L_delta 的全部顶点按零坐标数分类，核源可行性。全频率无商由 §5 的解析证明给出，有限高频点不认证极限。全源 nu_b 未产生空间样本，实际 D pressure 仍标 NOT_EXECUTED。

专属文件：full_rank_source_probe_20261007.py、full_rank_source_probe_20261007_registration.json、full_rank_source_probe_20261007_guard_results.json；执行与 hash 收据另存。旧两相位结果全部封存，本稿没有覆盖其中任何记录。
