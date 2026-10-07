# Aggregate mass-band：固定参考支付、高源剩余与真实整体例子

2026-10-07；父任务指定 6.1-sol high。结论分三项：固定参考响应支有严格 source-once 交通费；参考不足门能保留一个明确的高源剩余核，并只要求其单 cutoff 尾而非逐 pair 衰减；一个完整高源 slab 在这个门下占 Θ(n) 个 mass bands，但整体 source profile≤1/4。尚未证明任意输入的剩余尾，也未付一般 geom。

采用已读 L02 的真实输入/资格重建流程、A01/E01 的同一完整来源守恒。查重 [direct mass-band](direct_capture_mass_band_square_20261007.md)、[gated occupancy](gated_occupancy_bridge.md)、[source-tail](source_occupation_tail_geometry_20261007.md) 和 [source-optimal coupling](source_optimal_acceptance_coupling_20261007.md)，检索本课题目录的 fixed-scale/anchor baseline。未找到本文的固定参考好支与参考核正剩余完整组合；单 cutoff 尾的充分方向则早已在 source-tail §3、source-optimal §4 证明，不能重新记为新结论。逐 pair 质量档差 Gaussian 衰减已被 [真实赢家独审](direct_mass_band_cross_audit_20261007.md) 否定；本文不假设它。

## 1. 预固定尺度：支付完整响应，而不只是体积

冻结完整 μ=f dy≥0、W=||f||₁，原 true winner R(x)∈[a,b]、M(x)=h_R(x)*μ(x)、m(x)=M(x)R(x)^n，以及原 E⊂{M>τ}。不重选 winner 或任何原 FIRST 门。先取 U(x)=A_af(x)=h_a*f(x)，K>1，并切

\[
 E_A=\{x\in E:M(x)\le K U(x)\},\qquad
 E_D=E\setminus E_A.
\]

因为 U 是一份固定 probability kernel 的完整源响应，

\[
 \boxed{\tau|E_A|\le\int_{E_A}M(x)dx
                    \le K\int U(x)dx=KW.} \tag{1}
\]

若实际硬交通是 ∫E_A∫h_R(x−y)g(x,y)dν_h(y)dx，且原已经积分后的真实 g∈[0,1]、ν_h≤μ，则它也≤KW。这是真交通支付，与任意 source-square 假设无关。

多个预固定尺度合法的总账方式是：选同一概率权重 w_ℓ≥0、Σw_ℓ=1，令 k₀=Σw_ℓh_Lℓ、U=k₀*f。也可用一份 probability measure 对 L 积分。尺度/权重可由整份 input 与阈值预选，但不能因 receiver x 改变。仍有 ∫U=W，(1) 仍仅为 KW。

若改用 U=max_ℓA_Lℓf，只有 ∫U≤J W，一般好支费是 KJ W；不能让每个 anchor 另取得 W 而仍报 KW。对有限 J 个原尺度，等权混合使 M≤J U，故 K≥J 时整支都已付；下面的例子不反驳这种有限尺度基本支付。

## 2. 参考不足门给出的真实量化结构

单 anchor a 的 E_D 上有

\[
 \frac{\mu(Q(x,a))}{m(x)}
 <\frac1K\left(\frac a{R(x)}\right)^n
 =\frac{e^{-\beta(x)}}K,\quad\beta=n\log(R/a). \tag{2}
\]

这是 complete source 的实际捕获比例，不是 source 密度或角向独立性。原 R=a 的 winner 全被好支支付。剩余来源绝大部分在 query 的 Qa 外；但 (2) 并未限制实际质量 m(x) 的 band 个数。

可以对任意上述共同 k₀ 写同一个准确的正核分解：

\[
 h_R=\min(h_R,k_0)+(h_R-k_0)_+.
\]

在 E_D 上的原 weak-normalized source profile 分为

\[
 C(y)=\tau\int_{E_D}\frac{\min(h_R(x-y),k_0(x-y))}{M(x)}dx,
\quad
 V(y)=\tau\int_{E_D}\frac{(h_R(x-y)-k_0(x-y))_+}{M(x)}dx. \tag{3}
\]

由于 M>τ、∫k₀=1，C(y)≤1。完整来源 first moment 则满足

\[
 \int C\,d\mu\le\tau\int_{E_D}\frac{U(x)}{M(x)}dx
                         <\frac1K\tau|E_D|. \tag{4}
\]

这是实际的可吸收 reference 部分，并同时有 ∫C≤W；不是把任意共源 reference 当成正交投影。单 anchor 且 R≥a 时 min(h_R,h_a)=R^(-n)1_Qa，(3) 正是 Qa-inner 与 Qa 外 annulus 的分解。

不作 weak-normalization 的实际 g 门交通中，common-kernel 部分也有严格费≤W（min≤k₀，ν_h≤μ，∫k₀=1）。若原硬幅度为 τ<M≤2τ，则在 E_D 上还≤2τ|E_D|/K。该 core 与旧 gated cone 的核心 source-once 费属于相关同一机制，不能重领旧 paid 分支。

取 k₀ 为窗口内 kernels 的 probability mixture时，k₀≤k*=sup_Lh_L，剩余强包络是 k*−k₀，质量 Z−1，Z=1+nlog(b/a)。因而 kernel 分解本身只把原 O(n) strong envelope 减去1，未产生平方根。

## 3. 合用已核低密度支付后，一个真正更弱的 aggregate 合同

下文给纯 hard weak route 的完整校准；不宣称新 tail 已证。固定 ρ>0、ν∈(0,1)，在原 source 上一次切

\[
 f_{\rm lo}=f\mathbf1_{f\le\rho\tau},\quad
 f_{\rm hi}=f-f_{\rm lo},\qquad W_{\rm lo}+W_{\rm hi}=W.
\]

把 E_D 中 M_□f_lo>ντ 的 receiver 作为 E_L；由本地已核 dimension-free hard L² 常数 C_B，

\[
 \tau|E_L|\le C_B^2\rho\nu^{-2}W_{\rm lo}. \tag{5}
\]

本地口径见 [critical soft proof §1](../critical_soft_l2/proof.md)，以及 [whole-joint 独审](../full_soft_l2_extension/noise_parameter_review.md) 的硬 C_B 输入。本稿未新增未核常数或外部引用。选 E_R=E_D\E_L 后，原 selected low response≤ντ<νM，因此 high selected response≥(1−ν)M。原 σ、R、FIRST、fullfuture 与 all-history 仍是完整 input 的对象；没有为 f_hi 重新选停时。

对剩余高源 μ_hi=f_hi dy，定义

\[
 V_{\rm hi}(y)=\tau\int_{E_R}
              \frac{(h_R(x-y)-k_0(x-y))_+}{M(x)}dx.
\]

对 μ_hi 的 first moment J_hi，由 min(h_R,k₀)*μ_hi≤U<M/K，准确得

\[
 \boxed{J_{\rm hi}:=\int V_{\rm hi}d\mu_{\rm hi}
                 \ge(1-\nu-K^{-1})I_R,\quad I_R=\tau|E_R|.} \tag{6}
\]

因此一个足够且比强 square/逐 pair 衰减更弱的剩余合同是：只在一个 H>0 截点证明

\[
 \int(V_{\rm hi}-H)_+d\mu_{\rm hi}\le\delta I_R,
              \qquad\delta<1-\nu-K^{-1}. \tag{7}
\]

合同 (7) 若真成立，J_hi≤H W_hi+δI_R，故由 (1)、(5)、(6)，

\[
 \boxed{\tau|E|\le K W+C_B^2\rho\nu^{-2}W_{\rm lo}
              +\frac{H}{1-\nu-K^{-1}-\delta}W_{\rm hi}.} \tag{8}
\]

例如 K≈sqrt(n)、ρ≈sqrt(n)、H≈sqrt(n)polylog(n)、ν和δ固定且总和<1，会给目标级费用。**只有 (1)、(5)、(6) 已证明；任意输入的 (7) 仍未证明。** 这不是用新名字声明已付预算。相对旧 full-profile tail，具体进展是已付 reference 比较支、已付 low-density 输出支，以及一个不改原赢家、保留至少 1−ν−1/K 交通的真实高源正剩余核。

按原完整 m(x) 的 D_k 分割 V_hi，single-band 仍≤原 T_k≤2，全部 diagonal≤2J_hi≤2I_R；(7) 允许这一部分的 mass bands 彼此高度相关，也不要求 ∫V_hi²是 polylog·I_R。

对 actual heavy，(1) 的交通支付直接合法，common-kernel 支也直接合法；但 (6) 的 row lower bound属于完整 hard source，不能自动转给带 FIRST/CPGP/LCA 的 g 子交通。若以 V_hi^g替换，其 first moment一般更小。要向原 geom 主账迁移，须使用原未付交通自己的下界/预算校准，并分配其真实吸收常数；不能把ν=δ=1/4的纯 weak 示例直接塞入已有较小吸收账。

## 4. 真正 anchor-cold 的全高源连续赢家，占线性 bands 却整体小

下面不是 uniform input，也不是仅源矩阵反例。取 n≥2、δ=1/(100n)、D=16n²、τ=B=3/8，完整输入

\[
 f(y)=\delta^{-1}\mathbf1_{[-\delta/2,\delta/2]}(y_1)
                \prod_{j=2}^n\mathbf1_{[-D/2,D/2]}(y_j),
 \qquad W=D^{n-1}. \tag{9}
\]

所有正来源密度均为100n>τsqrt(n)，故 ρ=sqrt(n) 的 f_lo恰为零；未加背景质量。使用整个连续尺度窗口 [1,2] 的真实 winner。

receiver 核心为

\[
 \frac{4/3-\delta}{2}<|x_1|<\frac{2-\delta}{2},
 \qquad |x_j|<(D-2)/2\quad(j\ge2). \tag{10}
\]

这里所有 tangent queries都完整在 support 盒内，而 normal slab 直到 R_*(x)=2|x₁|+δ 才全部捕获。部分捕获阶段的平均为

\[
 \frac1R\frac{R/2-|x_1|+\delta/2}{\delta},
\]

其 R 导数是 (|x₁|−δ/2)/(δR²)>0；全捕获后平均为1/R，导数严格负。因此整个窗口上唯一 true winner就是

\[
 R(x)=2|x_1|+\delta\in(4/3,2),\qquad
 M(x)=1/R(x)\in(1/2,3/4),\qquad m(x)=R(x)^{n-1}. \tag{11}
\]

Aa f(x)=0，因为 |x₁|>1/2+δ/2。故 (10) 对 **任意** K>1 都保留在单 anchor E_D，同时满足原严格阈值及上幅度带 τ<M≤2τ。没有用 near-uniform receiver 去猜这个非均匀门。

核心的完整捕获质量比为 (3/2)^(n−1)，因此实际非空 mass bands为 Θ(n)。但对每个 source y，原 source profile有

\[
 S_{\rm core}(y)=\tau\int_{E_{\rm core}}
                    \frac{\mathbf1_{Q_x}(y)}{R(x)^{n-1}}dx.
\]

固定 x₁ 后，R仅依 x₁，故 tangent center积分≤R^(n−1)，正好抵消分母。双侧 normal receiver总长度是2/3，于是

\[
 \boxed{S_{\rm core}(y)\le\tau(2/3)=1/4,
 \quad\int S_{\rm core}^2d\mu\le I_{\rm core}/4,
 \quad I_{\rm core}\le W/4.} \tag{12}
\]

对 tangent source |y_j|≤(D−4)/2，上面每个 slice都取满，故 S_core=1/4。在核心质量档 k内令 w_k为该档双侧normal receiver长度；该档 I_k=τw_k(D−2)^(n−1)，而这个 inner source上 T_k=τw_k。所有档在同一来源上都正；整体小费用来自 w_k 的总和，而不是跨档正交或指数衰减。

式 (12) 仅支付 E_core 这个原 receiver 子集；其中 I_core=τ|E_core|，不是完整 E={M>τ} 的弱比或全空间费用。

完整位于 (4/3,2) 的 dyadic massband满足

\[
 w_k=r_k(2^{1/(n-1)}-1)
       \ge\frac{(4/3)\log2}{n-1}\ge\frac2{3(n-1)}.
\]

因此相距 Θ(n) 的两个完整档有真实 normalized cross下包

\[
 \frac{\int T_kT_l\,d\mu}{\sqrt{I_kI_l}}
 \ge\tau\sqrt{w_kw_l}
         \left(\frac{D-4}{D-2}\right)^{n-1}
 \ge\frac{\tau\,2}{3(n-1)}
         \left(\frac{D-4}{D-2}\right)^{n-1}. \tag{13}
\]

这只是 O(1/n) 下包，但对线性档距仍大于任何固定 C exp(−cΔ²/n) 的大 n 量级。单 anchor gate不能复活已失败的 pair-decay；同时 (12) 说明这种失败完全兼容常数 aggregate预算。本文不把这个特定 cancellation推广给任意相关输入：一般 R可依全部 coordinates，tangent积分便不能这样抵消。

这是一份合法结构例子，不是一般 weak反例。没有认证它满足原 soft FIRST/fullfuture、CPGP、LCA与全部历史门。它证明 hard gate M>K Aa本身不限制 massband 数，也不使 massband source supports正交；它没有反驳任何额外保留实际门的 aggregate尾合同。

## 5. 新三轮证书与剩余接口

先登记后执行一次：[registration](aggregate_mass_band_budget_20261007_registration.md)、[guard](aggregate_mass_band_budget_20261007_guard.py)、[results](aggregate_mass_band_budget_20261007_results.json)、[receipt](aggregate_mass_band_budget_20261007_receipt.json)。三轮 n=16,64,256，93个Fraction/整数断言 PASS。B=3/8 时，核心实际非空质量档的索引范围分别为[7,16]、[27,64]、[107,256]；这两个端档可以被核心质量范围部分截断。用于 (13) 的完全包含的两极档是 k_lo+1 与 k_hi−1，guard 核验的档距为 k_hi−k_lo−2，分别为7、35、147。I_core/W近似为0.248175、0.249520、0.249878，全都严格小于1/4；这些仍只是 E_core 子集的弱比。

守卫核原source mass/高度、query包含、partial/full导数符号、strict幅度、anchor零、线性档距、inner-source margin和解析aggregate cap。没有积分任意空间 input、没有 MC、没有优化放松矩阵，也没有把这三条input的零截断尾解释成一般阶数数据。不存在旧数值重跑。

当前可操作的下一接口是 (7) 在 **原共同参考不足门、真实 high-source限制及正核剩余** 上的 tail，而不是为所有 massband pairs制造通用衰减。已付分支和量化保留率明确；仍缺利用真实 winner/原全未来几何控制 high-source 总负载的方法。只凭 (2) 的 anchor-cold和原 n档窗口不够；还需要控制整份实际 cross流量或其尾，不能从档位个数、normalized Gram或例子 (12) 推出一般 sqrt(n)n^{o(1)}。
