# 来源抽样的无偏 acceptance count：真实 winner 的 Palm 边际与未付尾部

2026-10-07。只新增本稿及同前缀证书，不改原主账。读取 original_ordered_martingale_dilation、weak_normalized_projection 与当前 actual \(R_{\rm angle}\) 合同。使用 I02 的滤过/补偿子区分、E04 的有限停止与未停收益约定、L02 的真实模型/代理模型边界，已读 SKILL 及相关 method/provenance。统计部分另按 H05 的独立范围合同给 CI；只用自含有限条件期望、正 Fubini 与有界独立试验，不引用未核连续 Snell 或自由熵预算。

本轮固定同一个原 soft family \(P_r=K_r^{1/2}\) 及其真实 Markov 增量。一般 receiver sensor 允许变时，但它不改变 soft 增量。新表示可以不铺设 n 维接收网格；没有因此获得 weak 常数或实际 moving-soft 费用。

## 1. 固定 common soft chain 与一般 receiver sensor

取有限网 \(0=r_0\le r_1\le\cdots\le r_J\le1\)，原 fixed \(c,L\) 的
\[
 P_j=P_{r_j},\quad V_j=V_{r_j,r_{j-1}},\quad
 P_j=V_j\cdots V_1,\quad P_j^2=K_j.
\]
允许重复 r，用 identity increment；这对同一 soft 参数带多个 receiver sensor 有用。真实 P/V 的正性及 cocycle 已在原 dilation 稿直接证明，不能只由 Fourier symbol 非负断言正 Markov。

令 \(H_j\) **作用于测试函数**，为 bistochastic Markov kernel：
\[
 H_j1=1,\qquad \int H_j\phi\,dx=\int\phi\,dx.
\]
不要求 \(H_j\) self-adjoint 或与 \(P_j\) commute。原 even 卷积 sensor \(h_{L_j}^{(n)}\) 满足这些条件且 \(H_j^*=H_j\)。正确响应方向是
\[
 F_j=H_j^*K_j\nu_b,\quad
 M=\max_jF_j,\quad
 E=\{x\in\Omega^c:M(x)>\tau\},
\quad
 g_j(x)=\frac{\tau}{M(x)}\mathbf1_{A_j}(x),              \tag{1}
\]
其中 \(A_j\) 是原最大响应最小索引 tie rule 的不交赢家集合，\(\tau>0\)、\(\nu_b\ge0\)、\(W=\|\nu_b\|_1>0\)、\(\nu_b\) 支撑原 \(\Omega\)。\(g_j=0\) 在 E 外，零 M 无需相除。于是 \(0\le g_j\le1,\ \sum g_j\le1\)，并
\[
 g_jF_j=\tau\mathbf1_{A_j}.
\]
\(\int F_j=W\) 保证 \(|E|\le JW/\tau<\infty\)。这里取真正 \(F_j\) 的 maxweight；不能用旧未加 sensor 响应当分母。

## 2. 从原来源概率开始，条件独立接受

先抽 \(X_0\sim\nu_b(dx)/W\)，沿原 \(V_j\) 抽 \(X_1,\ldots,X_J\)。给定整条路径，独立抽每个新的 receiver：
\[
 X_j\xrightarrow{P_j} Z_j\xrightarrow{H_j}Y_j.
\]
再给每个 j 一枚独立 Uniform coin \(U_j\)，令
\[
 B_j=\mathbf1_{\{U_j\le g_j(Y_j)\}},\quad
 N=\sum_{j=1}^JB_j,\qquad
 \pi_j=q_j(X_j),\quad q_j=P_jH_jg_j.                    \tag{2}
\]
给定整条路径，B_j 是独立 Bernoulli(\(\pi_j\))；**不条件化时它们一般不独立**。Y_j 是新增的接收样本，不是原来源 X0 或原路径停止端点。

Fubini/对称性准确给
\[
 \begin{split}
 \mathbb EN
 &=W^{-1}\sum_j\int(P_j\nu_b)(x)(P_jH_jg_j)(x)dx\\
 &=W^{-1}\sum_j\int g_j\,H_j^*K_j\nu_b
 =\boxed{\tau|E|/W}.                                  \tag{3}
 \end{split}
\]
这不是假定 source distribution 等于 Lebesgue。

更正式地，Lebesgue path measure
\[
 d\mathfrak m=dx_0\prod_jV_j(x_{j-1},dx_j)
\]
是 σ-finite；\(\mathbb P_\nu=W^{-1}\nu_b(X_0)\mathfrak m\) 是概率测度。加入 receiver kernels/coins 后仍为概率。所有 (3) 的非负积分由 Tonelli 合法；有限 J 下 N≤J，后续 signed 有限和也可积。Lebesgue stationarity 只用于 kernel 对偶和 \(\int q_j=\int g_j\)，没有“均匀抽一个 Rn 点”的操作。

源换测度后原 Lebesgue reverse martingale 不再以相同公式成立。例如在 \(\mathbb P_\nu\) 下，
\[
 \mathbb E_\nu[a(X_0)\mid\mathcal F_j]
 =\frac{P_j(a\nu_b)(X_j)}{P_j\nu_b(X_j)}
\]
在正分母处成立。不能把 \(P_jg_j(X_j)\) 当作 \(\mathbb P_\nu\) 下原 \(g_j(X_0)\) 的条件期望；后者因源支撑/外域可为零。式 (2) 的 \(\pi_j\) 是**新 receiver sensor coin**的接受概率，正是这一点使无偏表示合法。

## 3. Stop、deficit、pair 与尾和，全部保留同一来源

给定路径的 PGF 为 \(\prod_j(1-\pi_j+\pi_j z)\)。因此
\[
 \begin{split}
 W\Pr(N>0)&=W\,\mathbb E_\nu[1-\prod_j(1-\pi_j)]
              =T_{\rm stop}\le W,\\
 W\mathbb E(N-1)_+
 &=W\,\mathbb E_\nu[\sum_j\pi_j-1+\prod_j(1-\pi_j)]
              =D_\tau,\\
 W\mathbb E\binom N2
 &=W\sum_{i<j}\mathbb E_\nu[\pi_i\pi_j]
              ={\rm Pair}_\tau.                         \tag{4}
 \end{split}
\]
这里 Tstop 与 D 是 weak_normalized_projection 中同一接受 union/overlap 的数值，不因把 source-probability 表示出来而改变。first-coin 的观察方向可变，union 事件不变；不要将其称为原 actual FIRST。

记 \(p_m=\Pr(N\ge m)\)，准确有
\[
 \mathbb EN=\sum_{m=1}^Jp_m,\quad
 D_\tau/W=\sum_{m=2}^Jp_m,\quad
 {\rm Pair}_\tau/W=\sum_{m=2}^J(m-1)p_m.                 \tag{5}
\]
所以尾预算可比强 pair 预算弱；没有因此自动获得任何 p_m 的维数界。

## 4. 补偿与首次接受后的真实未来 count

按 forward 顺序观察。设 \(\mathcal G_{j-1}\) 已含 \(X_0,\ldots,X_{j-1}\) 和此前 receiver/coins；观察 X_j 后、抽 receiver/coin 前另设 \(\mathcal G_{j-1/2}\)。正确两层补偿为
\[
 \mathbb E[B_j\mid\mathcal G_{j-1/2}]=q_j(X_j),\qquad
 \mathbb E[B_j\mid\mathcal G_{j-1}]
       =(V_jq_j)(X_{j-1}).                              \tag{6}
\]
\(N_j-\sum_{\ell\le j}(V_\ell q_\ell)(X_{\ell-1})\) 是有限 bounded-integrable martingale。另在增量揭示与 coin 揭示分成两步的 filtration 中，coin-centered 和 \(\sum(B_j-q_j(X_j))\) 也是 martingale；不能省掉所用 filtration。

令 \(T=\min\{j:B_j=1\}\)，允许 T=∞、未停收益0。这是真实辅助 forward 停时。有限可选抽样只给
\(\mathbb E\sum_{j\le T}(V_jq_j)(X_{j-1})=\Pr(N>0)\le1\)；它不支付停止后的接受次数。

把存活 probability source 按
\[
 a_0=\nu_b/W,\quad
 \ell_j=V_ja_{j-1},\quad
 \zeta_j=q_j\ell_j,\quad a_j=(1-q_j)\ell_j               \tag{7}
\]
递推，\(\zeta_j\) 正是 \(\Pr(T=j,X_j\in dx)\)，\(\sum\|\zeta_j\|\le1\)。令
\[
 R_j(x)=\sum_{k>j}V_{r_k,r_j}q_k(x),\quad
 R_J=0,\quad R_j=V_{j+1}(q_{j+1}+R_{j+1}).
\]
则保持所有过去 first-accept 排除后，严格
\[
 \boxed{D_\tau/W=\sum_j\int R_j\,d\zeta_j.}             \tag{8}
\]
这不是强 pair：过去次数只计一次 first-accept measure，总质量≤1。未知的是这份停止来源上的真实未来 count R_j；将它套一个未知原极大常数仍会循环，不能凭 \(\sum\|\zeta_j\|\le1\) 宣称 \(\int R\,d\zeta\) 小。

有限尾事件的 Doob process 也可直接构造：给定阶段接受数 k、状态 x，后续至少 m-k 次的概率按 V/q 两项递推，得到 bounded martingale \(\Pr(N\ge m\mid\mathcal G_j)\)。它的初值是未知 p_m，而不是来自 W 的统一小数。E04 的有限 Snell 对至多一次成功的0/1收益给≤1；若奖励改为累计 N，最优停止/终端值就是未付的 count，不能免费使用相同停止费。

## 5. Palm 接收边际真为 uniform Lebesgue，倒数 count 是弱合同

若 \(\mathbb EN>0\)，定义 accepted-occurrence Palm law
\[
 \mathbb P^\#(d\omega,j)
       =\frac{B_j(\omega)}{\mathbb EN}\,\mathbb P_\nu(d\omega).
\]
它下的 N 是 size-biased count，且选中 receiver Y_j 有准确边际
\[
 \Pr^\#(Y_j\in dy,\ {\rm index}=j)
 =\frac{g_j(y)F_j(y)}{W\mathbb EN}dy
 =\frac{\mathbf1_{A_j}(y)}{|E|}dy.                     \tag{9}
\]
因此 Palm receiver 为 Lebesgue|E 的 uniform **概率**；这是源采样经接受后验得到的结论，不能倒过来把初始 source law 当 Lebesgue。

准确 identities 为
\[
 \mathbb E^\#\frac1N
 =\frac{\Pr(N>0)}{\mathbb EN},\quad
 \mathbb E^\#\frac{N-1}N=\frac{D_\tau}{\tau|E|},\quad
 \mathbb E^\#(N-1)=\frac{2\,{\rm Pair}_\tau}{\tau|E|}.
                                                               \tag{10}
\]
新的较弱充分预算格式是：**只对原 truewinner/原 source** 证明
\[
 \mathbb E^\#(1/N)\ge1/C_n
 \quad\Longrightarrow\quad \tau|E|\le C_nW.            \tag{11}
\]
或只需 \(\Pr^\#(N\le K)\ge1-\delta\)，就有 \(C_n=K/(1-\delta)\)。它不要求 E N² 的强预算。也可写成 source sampling 可测的截尾条件
\(\mathbb E[N\mathbf1_{N>K}]\le\delta\mathbb EN\)。

这只是更贴近弱 volume 的真实新随机合同，**不是付款证明**；倒数身份本身没有减少复杂度。Palm source 也被 N 重新偏置，不再是 \(\nu_b/W\)。

## 6. 饱和目前只给 signed 尾势，任意 sensor 的免费小尾确实失败

令 \(c_m(y)=\Pr(N\ge m\mid X_0=y)\)。在 Lebesgue path measure 下由 \(m\mathbf1_{N\ge m}\le N\) 与 bistochastic sensor 得
\[
 0\le c_m\le1,\qquad
 \int c_m\,dy\le\frac1m\sum_j\int g_j
       =\frac1m\int_E\frac{\tau}{M(x)}dx.               \tag{12}
\]
对原饱和 \(Su=\nu_b-\mu_b,\ 0\le\mu_b\le\kappa\)，在原稿有界 S、u∈L1 或已合法的 duality 域中，
\[
 Wp_m=\int\mu_bc_m+\int uSc_m
 \le\frac{\kappa}{m}\int_E\frac{\tau}{M}
                     +\int uSc_m.                    \tag{13}
\]
未付的是 signed 尾势 \(\int uSc_m\)。逐 m 加 (13) 反而引入 harmonic 和，不能代替已有更好 \(\int C_{\rm def}\le\int_E\tau/M\) 的整体 cap 身份。饱和没有把 (11) 给出的 Palm 小 count 自动变成真结论。

以下严格排除 **任意 sensor 类** 的免费小尾，未排除 centered h_L。有限循环群 \(Z_{2J+1}\) 取 source \(\nu=\delta_0,W=1,\Omega=\{0\}\)、\(P_j=I\)。给不同非零 a_j，取 self-adjoint reflection sensor
\[
 (H_j\phi)(x)=\phi(a_j-x).
\]
它 bistochastic、可逆、无算子方向歧义；\(F_j=\delta_{a_j}\)，严格阈 \(\tau=1/2\) 下 M=1、A_j={a_j}、\(g_j=\tfrac12\mathbf1_{a_j}\)。source path 固定于0，而每个 receiver 永远为 a_j，故
\[
 N\sim{\rm Bin}(J,1/2),\quad \mathbb EN=J/2,\quad
 \mathbb E^\#(1/N)=\frac{2(1-2^{-J})}{J}.
                                                               \tag{14}
\]
这个反向压力保持真正响应 winner 与严格阈，但并非原 centered sensor、非非退化原 G、非指定饱和坏来源或实际 FIRST。它证明 (11) 的常数不可能只依 Markov/bistochastic 与 winner 归一化；必须使用 centered 几何或其它原结构。不能用它否定原 fixed G 的 weak 目标。

## 7. 原来源抽样的可证明 CI 设计

若原来源、P/V/sensor sampler 及原 truewinner g_j evaluation 精确，则 M 个**独立完整联合试验**给 \(N^{(a)}\in[0,J]\)，直接估计
\[
 A=\mathbb EN,\quad s=\Pr(N>0),\quad
 d=\mathbb E(N-1)_+,\quad p=\mathbb E\binom N2.
\]
对应目标为 \(WA=\tau|E|,\ Ws=T_{\rm stop},Wd=D_\tau,Wp={\rm Pair}\)。不同 trial 独立；同一 trial 中的接受币无条件不能假独立。已知范围分别 \(J,1,J-1,J(J-1)/2\)。Hoeffding 的有限有界变量证明由单变量中心化 mgf≤exp(t²B²/8)与 trial 独立乘积得到：
\[
 \Pr\{|\widehat a-\mathbb Ea|>
 B\sqrt{\log(8/\delta)/(2M)}\}\le\delta/4.
\]
四个量 union bound 后同时覆盖至少 \(1-\delta\)。无需估来源 entropy；也没有免费消除 J 的采样范围。

更适合 count-tail 的纯有理 CI：每个 m 记录 \(S_m=\#\{a:N^{(a)}\ge m\}\)。它在独立 trials 间是 Bin(M,p_m)，不同 m 的 S_m 相互依赖无妨。取 \(\sum_m\delta_m\le\delta\)，用精确 binomial-tail inversion 得同时区间 \([\ell_m,u_m]\)。可选择预先有理 endpoint grid，核
\[
 \Pr_{{\rm Bin}(M,\ell_m)}\{S\ge S_m\}\le\delta_m/2,\quad
 \Pr_{{\rm Bin}(M,u_m)}\{S\le S_m\}\le\delta_m/2
\]
对应的最大合法 lower/最小合法 upper；S_m=0 下 lower0，S_m=M 下 upper1。单调 tail inversion 直接给 coverage，选 grid 只使区间更宽。用 (5) 相加得到 A/d/p，s=p1；可用 p_m 单调关系作 intersection，但不能当各 m 独立。

approximation 与 sampling error 分账：若整个 joint trial law 与目标 TV≤ε_sampler，则任一 [0,B] statistic 的均值误差≤Bε_sampler。uniform kernel/source/coin误差可通过逐步耦合相加；没有证明这种 TV 时，CI 仅覆盖近似 sampler 的均值。原 max/strict-threshold/tie oracle 的误差尤其不能凭 response 浮点误差自动换成 coin uniform 误差；应认证 g_j上下围，模糊 receiver 行保守成区间并用同 coin \(N^-\le N\le N^+\)。CI 与 kernel truncation/response oracle误差分别保存，不能把蒙特卡洛CI当原 G 认证。

本稿只给上述可证明设计，**不执行原来源 Monte Carlo**，也没有声称现有 arbitrary L1 来源与原核已具有效率足够的精确 sampler/oracle。

## 8. 真实价值与 actual moving 接口范围

fixed c/L 的 \(K_r\) 正 Bernstein 全参数 strong L1 \(O(\sqrt n)W\) 已在旧稿保留；本文不把重新得到 fixed \(D_\tau=O(\sqrt n)W\) 当突破。新的用途是 count-tail/Palm 直接针对 receiver weak volume，能将任意变时 outer **hard sensor** 同源纳入表示，且 source-stop 仍≤W。

但原 actual soft 物理尺度 \(L_s(x)\) 同时改变 seed \(h_{L_s}\) 与 \(G_{c_t,L_s}\)。若只外面 \(h_{L_j}\) 变化而 soft G/P 链共同，则 (1)–(14) 合法；若 soft \(G_{L_j}\) 也变，尚无同一 source path 的正 cocycle，不能逐 j 重领 \(\nu_b/W\) 后称 union-stop≤W。source时间 c_t、原 ON/FIRST history、两来源/LCA/seed/fullfuture门也不能由 sensor 的名字自动继承。

固定共同 base kernel 若另有原 joint RN接受密度≤1，可用额外接受保留门，但 source/history依赖使 \(\pi_j\) 变成 augmented-state 条件概率，不再是裸 \(P_jH_jg_j\)；必须重新核响应与 (9) 的真实分母。没有该 canonical common law 时，不能按自由 sensor 抽样假装通过 actual门。

同 r 跨物理 soft 尺度直接使用原 P 商建共同链，还已有具体解析障碍：[full_rank_source_probe_20261007.md](full_rank_source_probe_20261007.md) §5 在 \(0<r<1\) 证明两个方向都不是 Markov convolution。小尺度至大尺度商在非零频率严格小于1、高频却趋1，概率特征函数的 Cesàro 平均会迫使全部质量在0，从而矛盾；反方向商大于1。这里只引用该证明，不重复其数值，也不扩展到 r=0/1、斜向 r/L 或多参数 dilation。该稿把 H 定义为 density 算子，故其 \(q=P H^*g,F=HK\nu\) 与本文 test-operator 约定 \(q=PHg,F=H^*K\nu\) 一致。

因此 (11) 对原 centered sensor/common source 的小 count 合同，或其与 actual moving soft/完整门可兼容版本，仍是未付关键。当前原主账不变。

## 9. 新三轮精确守卫终态

先保存同前缀 registration，再只执行新 exact_guard 一次。未读取旧240 fixture；没有随机种子。有限空间为 \(\mathbb Z_2^d\)、counting measure，\(G_i=(2I+\mathrm{flip}_i)/3\)。令 P 的每轴 Walsh eigenvalue \(\rho_j=1-j/(3J)\)，\(r_j=\tfrac32(1-\rho_j^2)\)，则精确 \(P_j^2=K_{r_j}\)，增量每轴 eigenvalue \(\rho_j/\rho_{j-1}\) 为正 Markov。新 u 在0、1两点取1、2，按外域 \(-Su\) 最大值选 κ，再构造饱和 \(\nu,\mu\)。传感器取 identity 与变时对称 flip mixture；四阈值 κ/16、κ/8、κ/4、κ，响应/winner 全部从新 fixture 重算。

24个 fixture 用 forward state/count recursion 与独立 backward PGF 给完整 N law，并另核 pair 双时积分、union/overlap recursion、first-accept future count (8)、source-conditional signed tail、Palm uniform intensity、条件整路径 coin 的 PGF。exact CI guard 用8个独立 trials 的理论 binomial law、1/16有理 endpoint grid、每尾 \(\delta_m=1/(20J)\)，对登记真参数枚举全部观察计数核 coverage；这不是执行来源 Monte Carlo。3个空 strict event 原样保留，无除零或调参。

|轮|d|J|全部精确 checks|空 event fixture|reflection \(\mathbb E^\#1/N\)|
|---|---:|---:|---:|---:|---:|
|1|3|3|579|1|7/12|
|2|4|6|835|1|21/64|
|3|5|9|1097|1|511/2304|

总计 **2512/2512 Fraction PASS**，含注册存在检查1项；运行约0.53秒，全部终态，无 live session。reflection 每轮另在 \(\mathbb Z_{2J+1}\) 精确核 self-adjoint/bistochastic 反射、完整 Bin(J,1/2) law 和 \(N^\#=1+\mathrm{Bin}(J-1,1/2)\)。它认证 (14) 的任意 sensor 类边界，不能外推成 centered h 或实际障碍反例。

专属文件为 [registration](acceptance_count_projection_20261007_registration.json)、[exact guard](acceptance_count_projection_20261007_exact_guard.py)、[results](acceptance_count_projection_20261007_results.json)。执行 script SHA256 为 95d6e45c797b38ae763be77863e0fa51a5481fc6068d1c75e5fb962d17bcc698；registration SHA256 为 9588c3b8daec3597ca77be088dac8d5b534ae546771015a1c064ebbbd0c6b1ed。这些是有限条件期望/coin/CI 的独立有理认证，不是原 \(\mathbb R^n\) G 或 centered-h_L 样本，不认证原 fullhistory、LCA、FIRST、moving-soft 或目标阶数。固定 family 的旧 \(O(\sqrt n)\) 正包络也未重跑或新领费用。

本轮最终剩余合同仍为原 common-source/truewinner 下的 (11)，或 (8) 中 first-accept 来源上 future count 的有偿估计；实际 moving soft 还须先给兼容共同路径。Palm 恒等、停止来源质量≤1和置信区间均不自动完成这两项。
