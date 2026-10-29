# Count posterior通道：无表数损失的固定尺度L1界

2026-10-29。先运行 `probe.py`、`trace_checks.py`、`multitable_checks.py`，保存NumPy与Fraction数值，再写本推导。以下证明自含，不使用Gaussian近似或外部未核验定理。

## 1. Binary input与signal fraction

一个binary-input通道由 C~Bernoulli(η)、query X及其联合律构成，0<η<1。令Q为X的边际律，p(X)=P(C=1|X)，u=η(1−η)。全期望和全方差给

\[
E_Qp=\eta,\qquad f={\operatorname{Var}_Q p\over u}\in[0,1]. \tag{1}
\]

f是由query揭示的binary variance fraction。n个坐标复制该独立通道，S=ΣC_i的先验为π_s=binom(n,s)η^s(1−η)^(n−s)。给定query向量后，C_i条件独立、参数p(X_i)，所以T_s(X)=P(S=s|X)是**精确PB概率**。令L_s=T_s/π_s。

## 2. Count Gram只依赖f

令 Z_i=(C_i−η)/√u。对k=0,…,n定义

\[
e_k(S)={1\over\sqrt{\binom nk}}
\sum_{|D|=k}\prod_{i\in D}Z_i,\quad e_0=1. \tag{2}
\]

表达式对坐标置换不变，故确为S的函数。两个不同子集的乘积内积为0，因其对称差中某坐标只出现一次、EZ_i=0；同一子集范数为1。因此这些n+1个函数在Bin(n,η)上正交归一，又因S有n+1个正质量状态，它们组成全部count函数的正交基。

写 g(X)=(p(X)−η)/√u，E_Qg=0、E_Qg²=f。由坐标条件独立，

\[
v_k(X):=E(e_k(S)|X)
={1\over\sqrt{\binom nk}}\sum_{|D|=k}\prod_{i\in D}g(X_i).
\]

同样的子集正交性给 E_Q^n(v_jv_k)=1_{j=k}f^k。将singleton indicator在基(2)中展开，得到

\[
L_s(X)=\sum_{k=0}^n e_k(s)v_k(X),\quad
\boxed{E_{Q^n}(L_sL_t)=\sum_{k=0}^n e_k(s)e_k(t)f^k.} \tag{3}
\]

所以在**共同η、n固定、每坐标同一通道**下，完整count posterior Gram只依赖f，而非query通道其余细节。这不声称所有query概率律本身相同。

对任意count集合E，记π_E=P(S∈E)>0，a_{E,k}=E[1_E(S)e_k(S)]，L_E=P(S∈E|X)/π_E。则

\[
\boxed{\operatorname{Var}_{Q^n}(L_E)
={1\over\pi_E^2}\sum_{k=1}^n a_{E,k}^2 f^k.} \tag{4}
\]

特别是每一个bin的diagonal variance随f单调。系数平方非负；不能以完整Gram非对角元素的单调性替代这一步（非对角元素未必单调）。

若各坐标通道不同但binary prior均为同一个η，令各自signal为f_i。上述证明将f^k换为 `binom(n,k)^−1 Σ_{|D|=k}∏_{i∈D}f_i`，仍不超过(f*)^k，其中f*=max_i f_i。若坐标的η_i不同，此处count基与结论不能直接照搬。

## 3. 不同表、不同query倾斜的共同diagonal帽

令 E_1,…,E_R是{0,…,n}的任意固定分割，允许空段，π_b=P(S∈E_b)。对每个表b可采用**不同**binary-input通道，但η共同，signal f_b≤f*<1。其query律Q_b可以不同。记T_b为该表的第b个bin后验，L_b=T_b/π_b。空段项完全为0，不除π_b=0。

各a_{E_b,k}仅依赖共同Bin(n,η)及分割。由(4)，

\[
\sum_{b:\pi_b>0}\pi_b\operatorname{Var}_{Q_b^n}(L_b)
\le\sum_{k=1}^n(f^*)^k
\sum_{b:\pi_b>0}{a_{E_b,k}^2\over\pi_b}.
\]

每个bin内Cauchy给 a_{E_b,k}²/π_b≤E[1_{E_b}e_k²]；对分割求和为1。因此

\[
\boxed{\sum_b\pi_b E_{Q_b^n}(L_b-1)^2
\le\sum_{k=1}^n(f^*)^k\le{f^*\over1-f^*}.} \tag{5}
\]

这里未将各表放在同一个Q下计算，也未偷换成跨表相关性假设。它是分别积分后的diagonal和。相同证明允许表内坐标通道不同，只要所有binary priorη共同、各signal≤f*。

单个通道取全部singleton bins时(5)的有限几何和是精确值：

\[
\sum_s\pi_s E_{Q^n}L_s^2=\sum_{k=0}^n f^k. \tag{6}
\]

分段/coarsening只会降低这种加权trace；不降低所有可能的后验最大值。

## 4. 真正gated与prior的固定尺度L1界

假设各表query律Q_b^n相对共同基准测度μ有非负密度A_b，且∫A_b dμ=1。定义

\[
G=\sum_b A_bT_b,\qquad P=\sum_b\pi_bA_b.
\]

G和P各自积分为1：对G用 E_{Q_b^n}T_b=π_b及Σπ_b=1。正系数及两次Cauchy给

\[
\begin{aligned}
\|G-P\|_{L^1(\mu)}
&\le\sum_b\pi_bE_{Q_b^n}|L_b-1|\\
&\le\sum_b\pi_b\sqrt{\operatorname{Var}_{Q_b^n}L_b}\\
&\le\sqrt{\sum_b\pi_b\operatorname{Var}_{Q_b^n}L_b}
\le\boxed{\sqrt{f^*/(1-f^*)}}. \tag{7}
\end{aligned}
\]

这是不同表真实query倾斜下的结论，不要求A_b相同。若只考察一个通道、所有表相同，则G=P=A恰好完全抵消；后验max仍可很大，二者不矛盾。

**Cube实现。** 每表取非负单位周期密度ρ_b、共同源类别χ，假定∫ρ_bχ=η共同。先y~ρ_b，再独立取宽L的均匀位移得到query x（模1），就是一个binary channel，C=χ(y)。其query密度为a_b=z_{b,+}+z_{b,−}，posterior为p_b=z_{b,+}/a_b。n坐标乘积query密度A_b=∏a_b。令真实gated源 `F=Σ_b ρ_b^n 1_{Σχ(y_i)∈E_b}`，prior源 `Fprior=Σ_bπ_bρ_b^n`。固定同一个L时，真实cube平均分别恰为G、P。因此(7)直接适用于这一真实cube接口；不需要Gaussian转移。

## 5. 宏观L≥1的universal corollary

设m=floor L≥1。对任意表、任意query，把窗口分成m完整周期及余弧，其两类余弧质量满足0≤e_+≤η、0≤e_-≤1−η。于是

\[
{m\eta\over m+1-\eta}\le p\le{(m+1)\eta\over m+\eta},\qquad E_Qp=\eta.
\]

若随机变量p∈[a,b]且均值η，从E[(p−a)(b−p)]≥0直接得Var p≤(η−a)(b−η)。代入上述端点、u=η(1−η)，得到

\[
f_b\le {u\over(m+\eta)(m+1-\eta)}
={u\over m(m+1)+u},\qquad
{f^*\over1-f^*}\le{u\over m(m+1)}.
\]

故对这些宏观真实cube来源和任意count分割、任意n、R，固定L有

\[
\boxed{\|A_LF-A_LF_{\rm prior}\|_1
\le\sqrt{\eta(1-\eta)\over m(m+1)}
\le{1\over2\sqrt{m(m+1)}}.} \tag{8}
\]

m=1时上界1/(2√2)≈.353553。整数L时各p=η，真实界其实是0。所有量词都在共同η、product坐标、共同固定L、非负单位质量来源及count分割框架内；不是一般输入C_n定理。

第38轮连续query f与bin-channel f的区别：精确bin均值满足pbar=E_Q(p|query-bin)，故f_bin≤f_cont。若class核在每bin的range≤B/L，且a=z_++z_-≥1/L，则

`p−pbar=[(1−pbar)(z_+−m_+)−pbar(z_-−m_-)]/a`，

所以|p−pbar|≤B，全方差给

\[
0\le f_{\rm cont}-f_{\rm bin}\le B^2/[\eta(1-\eta)]. \tag{9}
\]

38轮参数可用B=29879/44945280、η∈[.49,.51]，得到加量≤22318866025/12620438518754304≈1.76846993×10^-6。η范围也可由 `a=127/255`、`|η−a|≤(1/1020)(245/81)`严格保证。该解析差额不能认证已保存float f_bin本身；约.01586只是数值读数。universal界(8)不依赖这些FP值。

## 6. Gaussian通道的max后验：有熵而无一致Lq矩

本节仅为独立Gaussian基准，**不是cube近似定理**。令Y,Z独立N(0,1)，X=τY+sZ，0<τ<1、s=√(1−τ²)。对X等先验quantile分段，H_R(y)=max_b P(X∈I_b|Y=y)/P(X∈I_b)。条件密度比为

\[
\ell_y(x)={1\over s}\exp\left[-{\tau^2x^2-2\tau xy+\tau^2y^2\over2s^2}\right].
\]

每个bin后验比是ℓ_y在该bin按先验的平均。最大值位于x=y/τ，为s^-1 exp(y²/2)，故H_R≤该值。沿dyadic R，分割细化使H_R不减；含x=y/τ的bin直径趋0，所以

\[
\boxed{\sup_m H_{2^m}(y)=s^{-1}e^{y^2/2}.} \tag{10}
\]

单调收敛和Gaussian积分给，对0<q<1，

\[
\sup_m E H_{2^m}^q=s^{-q}(1-q)^{-1/2};
\]

q≥1则该sup为∞。特别是每个任意小但非零τ都有R增长的平均后验熵来源；这不是“gated必无收益”。τ=0时H_R恒1，τ→0与R→∞不能交换。

与此同时包络的weak-L1仍为s^-1：对t>s^-1，写r=√(2log(st))，则 `tP(s^-1e^{Y²/2}>t)=s^-1 e^{r²/2}2Φbar(r)≤s^-1`；不等式来自φ(r+z)≤e^-r²/2 φ(z)对z≥0积分。t↓s^-1给等号。故非一致Lq矩不等于弱范数无界，也不等于真实gated weak gain。

## 7. 无dimension一致L1(max)的真实product-channel反例

回到任意固定binary channel，设f>0、d=ess sup p−ess inf p>0。令U=Σp(X_i)。给定C=1或0，p的均值分别为η+Var(p)/η、η−Var(p)/(1−η)，差恰为f。因此给定S=s，

`E(U|S=s)=n m0+s f`，`Var(U|S=s)≤n d²/4`。

第二式在给定完整label向量后由独立性及区间方差帽成立；均值只依赖count，故再条件到S不会产生额外方差。选

\[
D=\lfloor2d\sqrt n/f\rfloor+1,
\quad\mathcal S=\{0,D,2D,\ldots,\lfloor n/D\rfloor D\}.
\]

各中心nm0+sf相距严格大于2d√n，其半径d√n的区间互不相交；Chebyshev给每个条件律在本区间内概率至少3/4。L_s是query在条件S=s与Q^n之间的真实RN，所以

\[
\boxed{E_{Q^n}\max_s L_s
\ge {3\over4}(\lfloor n/D\rfloor+1),\quad
\liminf_{n\to\infty}{E\max_sL_s\over\sqrt n}
\ge{3f\over8d}>0.} \tag{11}
\]

任意q≥1的dimension一致Lq控制也随之失败。这里最大化包括极小先验质量的count状态。若要求整数quantile bins，取R>min(η,1−η)^−n即可保证每个count jump都被quantile网格击中，所有非空段成为singleton；η=1/2可直接取R=2^n。因此(11)给明确的n/R非一致矩反例，**没有证明多项式R即可达到√n**，也不将其改写成真实weak阶。

此反例可由第4节的真实cube channel实现，只要固定ρ、χ、L的p非恒定。特别是ρ严格正、χ为一个非空非全圆周弧、L非整数时p不恒定：令h=(χ−η)ρ，若其窗口积分恒0，则绝对连续积分的a.e.导数使h具有非整数平移周期L mod1；这要求h的正号弧在该非零平移下不变，不可能。于是f>0。此处n变化时可固定这一个一维来源，不声称满足原匹配n依赖参数资格。若所有表恰采用同一个来源，尽管max后验满足(11)，G=P仍完全相同。

## 8. 收益边界

(5)–(8)是本轮新的严格固定尺度接口，允许不同表的query倾斜，不付R损失。Gaussian与(11)则说明，若先取best posterior再忽略对应来源权重及prior，确会留下可增长的熵；不能用一个并不存在的dimension/R一致Lq(max)界消掉它。

本稿不控制 `∫sup_L|A_LF−A_LFprior|`，不允许每坐标/每表选不同cube尺度，亦不从逐L L1界推出weak maximal排序。J随Kv增长时跨尺度费用仍需独立支付；真实共同尺度的正log增益目标保持开放。
