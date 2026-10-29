# 真实 PB 条件均值的确定误差接口

2026-10-29。先运行 `probe.py` 并保存数值，再选用下列有原始文献支持的常数。本结论是独立非同分布 Bernoulli 和的误差接口，不把 Gaussian 过程当作 cube 来源。

## 1. 原始定理与可用常数

采用 Ilya Tyurin, *New estimates of the convergence rate in the Lyapunov theorem*, arXiv:0912.0726v1 (2009)，§1 式(1)和定理7（PDF第2、6页，证明第15–17页）。已读原文的独立、中心化、有限三阶绝对矩假设及其数值常数区分。原文给一般独立非同分布情形 **C=0.5606**；同分布的另一常数不用于本接口。引用的是文献定理，未独立重做其计算机辅助常数证明。较保守的 C=1 同样可用；不能仅由该文把0.5606向下取成0.56。

直接来源：[arXiv摘要及元数据](https://arxiv.org/abs/0912.0726)、[原始全文](https://arxiv.org/pdf/0912.0726)。检索及实际阅读范围见 `search-log.md`。

设独立 B_i~Bernoulli(p_i)，S=ΣB_i，q_i=1−p_i，

\[
\mu=\sum_i p_i,\quad \sigma^2=\sum_i p_iq_i>0,\quad
\beta_3=\sum_i p_iq_i(p_i^2+q_i^2),\quad
\epsilon=C\beta_3/\sigma^3.
\]

因为 E|B_i−p_i|³=p_iq_i(p_i²+q_i²)，上述独立非同分布定理直接给

\[
\sup_x\left|P(S\le x)-\Phi((x-\mu)/\sigma)\right|\le\epsilon. \tag{1}
\]

这里 x 为所有实数。某些 p_i=0 或1仍允许；若 σ=0，则 S 为确定整数，直接计算所有段概率。

## 2. 所有整数段的统一上下包

对整数段 E=[a,b]∩{0,…,n}，先与支持相交。空段概率精确为0。非空段令

\[
u=\begin{cases}-\infty&a=0,\\(a-1/2-\mu)/\sigma&a>0,\end{cases}\quad
z=\begin{cases}+\infty&b=n,\\(b+1/2-\mu)/\sigma&b<n,\end{cases}
\]

以及 G=Φ(z)−Φ(u)，k=1_{u有限}+1_{z有限}。S是整数，所以在有限半格端点的CDF恰等于所需整数CDF。分别应用(1)，无穷端点使用精确的0或1，得到

\[
\boxed{\max(0,G-k\epsilon)\le P(S\in E)\le\min(1,G+k\epsilon).} \tag{2}
\]

内段 k=2；两端尾段 k=1；全支持 k=0。半格 normal 质量只是选择CDF评估点，本式并未假设某个额外的 continuity-correction 改善定理。(1)已对所有实数统一，故(2)对所有整数段、所有允许参数、任意有限表族同时成立，不需再做概率并集或付统计 α。

若 p_i∈[p0,1−p0]，0<p0≤1/2，则

\[
\epsilon\le {C[1-2p0(1-p0)]\over\sigma}
\le {C[1-2p0(1-p0)]\over\sqrt{np0(1-p0)}}. \tag{3}
\]

证明是 p_i²+q_i²≤1−2p0(1−p0)，然后 Σp_iq_i≥np0(1−p0)。实际矩式通常比(3)更紧，计算费用只需遍历参数或不同参数计数。p0=1/4 时右侧为 5C/(2√(3n))；p0=1/3 时为 5C/(3√(2n))。

## 3. 宏观 cube 的合法 p 范围

固定非负周期一维单位质量密度 ρ、固定源类别 χ∈{0,1}，记 η=∫_0^1ρχ，假定0<η<1。对任意 query中心和 L=q+ℓ≥1，q=floor L≥1，每个窗口包含 q 个完整周期及一段余弧。令余弧两类质量 e_+,e_-，则

\[
0\le e_+\le\eta,\quad0\le e_-\le1-\eta,\qquad
p={q\eta+e_+\over q+e_++e_-}.
\]

分式对 e_+递增、对 e_-递减，故

\[
{q\eta\over q+1-\eta}\le p\le{(q+1)\eta\over q+\eta},\qquad
\boxed{{\eta\over2-\eta}\le p\le{2\eta\over1+\eta}}. \tag{4}
\]

每坐标条件于 query bin 后对 class 核取均值，同样先对 e_±取均值，界仍成立。若不同表的 η_b 都位于 [η0,1−η0]，统一可取

\[
p0={\eta0\over2-\eta0}. \tag{5}
\]

η=1/2 给 [1/3,2/3]。只有已知 |η_b−1/2|≤1/10，才可据(5)统一使用 [1/4,3/4]。不要求严格等于1/2，但 η_b 或其严格包围必须来自实际源质量计算。

重要区别：范围(4)不准许把真实 prior Bin(n,η) 换成 Bin(n,1/2)。哪怕η−1/2很小，标准化偏移约2√n(η−1/2)，也可能很大。若各表全周期类别质量不同，gated raw源的总质量 Z=Σ_b P_{ρ_b^n}(S∈E_b)还可能不等于1；归一化必须单独携带，BE本身不会修复来源质量。

## 4. 正加权真实条件均值

固定 bin信息及一个共同L。令 μ_{i,b,±} 是真实 class 核的条件均值，s_{i,b}=μ_++μ_->0，p_{i,b}=μ_+/s，A_b=∏_i s_{i,b}。正 multilinear 恒等式给真实的未归一化条件均值

\[
m=\sum_b A_b T_b,\quad T_b=P\{\textstyle\sum_i\mathrm{Bern}(p_{i,b})\in E_b\}. \tag{6}
\]

Bernoulli变量仅是多项式系数的精确表达，不是新增来源或 query 的概率假设。若某 s=0，整表项为0。对各表使用(2)得 T_b^-、T_b^+，则

\[
\boxed{m^-:=\sum_b A_bT_b^-\le m\le m^+:=\sum_b A_bT_b^+.} \tag{7}
\]

如用未截断的 normal和 g=ΣA_bG_b，还有确定绝对误差

\[
|m-g|\le\sum_b A_bk_b\epsilon_b\le2C\sum_b A_b\beta_{3,b}/\sigma_b^3. \tag{8}
\]

若同源实际 prior权重为π_b，且 π_b=0 的表项确实为空，和 prior平均 m_prior=Σπ_bA_b 比较，可得

\[
{|m-g|\over m_{\rm prior}}\le
\max_{b:\pi_b>0}{k_b\epsilon_b\over\pi_b}. \tag{9}
\]

当π_b约1/R且p统一远离0、1时，(9)为 O(R/√n)，而非自动与R无关。应优先使用加权实际矩预算(8)或截断包(7)。表权极小、或T_b极罕见时，绝对BE可能完全不能给正下包。

若参数按query bins计数，设 N_j为各bin坐标数，则 μ_b=Σ_jN_jp_{j,b}、σ_b²=Σ_jN_jp_{j,b}(1−p_{j,b})、β_{3,b}=Σ_jN_jp_{j,b}(1−p_{j,b})[p_{j,b}²+(1−p_{j,b})²]。因此这一步与n本身的DP长度无关；表乘积A_b用log表达，可构成大n接口。保存的普通浮点μ、Φ和logA不是区间证书，应用(7)作认证时还须向外支付参数、CDF和指数运算误差。未知真实Z还须以独立严格Z区间除开。

## 5. 与条件方差及共同尺度的连接

另线正多项式引理给 Var(F|bins)≤Vhat_L 后，若 m_L^->λ，便有

\[
P(F_L>\lambda\mid bins)\ge
{(m_L^- -\lambda)^2\over Vhat_L+(m_L^- -\lambda)^2}. \tag{10}
\]

否则取0。有限尺度上取这些下界的最大值合法：选择由bin信息确定的同一个L，再对隐藏位置应用Cantelli。它不假设各尺度独立。若F是raw未归一化源，而目标是F/Z，需先将阈值变为λZ或用Z的可靠上包。

## 6. 实际数值与边界

`results.json` 保存57个明确PB向量，n31/63/255各19个。最大Kolmogorov误差依次0.0948632358、0.0671305973、0.0335837582；最大误差/Lyapunov为0.3986816871。全部456个分箱满足所引用0.5606预算；这只是有限数值核验，不决定最优常数。真实37 box公式产生的本轮p向量位于约[0.4,0.6]；n255例为新有限query扩展，不是旧MC样本。

最稀有的保存段为n255、p_i=1/4的顶段：真实DP值约3.35969×10^-23，而BE绝对预算约0.05067134，下包为0。因而本接口可给真均值包，却不能仅凭绝对BE认证极罕见相对尾或任何弱范数增长。普通float DP也不是对上述小概率的舍入认证。

独立脚本 `check_bounds.py` 已实际运行：81组Fraction验证(4)，另一个3坐标8点正多项式例严格核验均值替换及方差帽；对另线12个保存小例的有理方差和分层平方根区间只读重算全部通过。详见 `checks.json` 和 `AUDIT.md`。
