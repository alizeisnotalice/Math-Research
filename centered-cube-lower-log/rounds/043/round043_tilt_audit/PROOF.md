# 第43轮独立证明：任意查询律的正确归一化trace与原先验付款障碍

本文件在本机实际运行 `checks.py` 后整理。已执行144个Fraction有限查询律例、24个精确histogram trace例、90个精确倾斜恒等式及90个旧先验Jensen界，另完整重放9000个数值查询律。来源数值探路先于本轮正式目标式；没有把扫描当成全参数证明。

## 1. 固定posterior、任意查询律

令一周期上的非负可积子密度ρ_c满足Σ_cρ_c=ρ、∫ρ=1，原class质量η_c=∫ρ_c>0，Ση_c=1。质量零的class可删除。取整数q≥1、0≤ell≤1，L=q+ell。令h为长度ell的周期指示核，定义

\[
r_c(x)=\int\rho_c(y)h(x-y)dy,\quad r(x)=\sum_cr_c(x),\quad
a(x)=\frac{q+r(x)}L,
\]
\[
p_c(x)=\frac{q\eta_c+r_c(x)}{q+r(x)}.
\]

该核与真实中心窗仅相差一个确定query平移；q奇数时为半周期平移。对于任意查询律ν，包括奇异概率测度，本文均使用这个显式公式定义的p_c版本，而非任意更改过零测集值的后验版本。因为0≤r_c≤η_c、0≤r≤1，公式在每个x上定义良好。

保持p_c不变，改用任意概率ν，定义新class先验

\[
\mu_c=E_\nu p_c>0,\qquad
\tau_\nu=\sum_c\frac{\operatorname{Var}_\nu p_c}{\mu_c}.
\]

这是新联合律P_ν(dx,c)=ν(dx)p_c(x)的class先验。一般μ≠η，不能继续声称该联合律的class边缘是原来源的η。

## 2. 一般标量引理及维数无关帽

若V∈[a,b]且0<a≤b，令m=EV。由(V−a)(V−b)≤0得

\[
\frac{\operatorname{Var}V}{EV}
\le a+b-m-\frac{ab}{m}
\le a+b-2\sqrt{ab}=(\sqrt b-\sqrt a)^2.
\]

第二步来自(m−√ab)²/m≥0。第一步在V只取a,b时取等；选择均值m=√ab即可同时取等。因此这是只知道区间[a,b]时的最优标量常数。

对V_c=p_c/η_c，有

\[
\frac q{q+1}=:a\le V_c\le b:=\frac{q+1}q,\qquad ab=1.
\]

所以

\[
\frac{\operatorname{Var}_\nu p_c}{\mu_c}
=\eta_c\frac{\operatorname{Var}_\nu V_c}{E_\nu V_c}
\le\frac{\eta_c}{q(q+1)}.
\]

求和得到自含的严格结论

\[
\boxed{\tau_\nu\le c_q:=\frac1{q(q+1)}<1.}
\]

本结论不要求ν相对于原Q有密度，也不要求ν是score倾斜；不依赖K、来源形状、class质量或ell。q=1给τ_ν≤1/2。标量引理取等并不证明本channel族能同时实现所有class的取等条件；本文不声称c_q是channel最优常数。第42轮在原Q下较强的1/[4q(q+1)]帽，不能未经证明就替换到任意ν。

## 3. 正确新先验下的histogram trace

定义

\[
H_{cd}=\frac{E_\nu[p_cp_d]}{\sqrt{\mu_c\mu_d}},\quad
u=(\sqrt{\mu_c})_c.
\]

H是Gram矩阵，Hu=u，H−uu^T为归一化协方差矩阵，因而半正定。还可由条件Jensen证明H≤I：对v∈R^K，令Z=v_C/√μ_C，在P_ν下有

\[
v^THv=E_\nu(E[Z|X])^2\le EZ^2=\|v\|_2^2.
\]

所以H有谱1,λ_2,...,λ_K，其中λ_j≥0、Σ_{j≥2}λ_j=τ_ν≤c_q。

取n个独立的同一个联合律P_ν。源端class计数N的正确先验是

\[
\pi_\mu(\alpha)=\frac{n!}{\prod_c\alpha_c!}\prod_c\mu_c^{\alpha_c}.
\]

其后验g_α(x_1,...,x_n)是多项式∏_iΣ_cp_c(x_i)z_c的z^α系数。g_α作为函数只取决于固定p；但归一化先验随ν改变。定义

\[
G^{(n)}_{\alpha\beta}
=\frac{E_{\nu^n}[g_\alpha g_\beta]}{\sqrt{\pi_\mu(\alpha)\pi_\mu(\beta)}}.
\]

给出精确谱式的理由如下。单坐标class函数1_{C=c}/√μ_c构成L²(μ)的正交基，其条件期望Gram是H。独立n坐标给H^{⊗n}。归一化histogram指示1_{N=α}/√π_μ(α)构成源端对称张量子空间的正交基，其Gram恰为G^(n)。把H的正交特征基张量后对称化，按次数β、Σβ=n编号的特征值为∏_jλ_j^{β_j}，包括λ_1=1。于是

\[
\operatorname{tr}G^{(n)}
=h_n(1,\lambda_2,\ldots,\lambda_K)
=\sum_{d=0}^nh_d(\lambda_2,\ldots,\lambda_K)
\le\sum_{d=0}^n\tau_\nu^d
\le\frac1{1-c_q}.
\]

h_d≤(Σλ_j)^d，因为右边展开的每个非负单项式的多项系数至少1。没有极限分布或Gaussian近似。q=1给

\[
\boxed{\operatorname{tr}G^{(n)}\le2\quad\text{对全部 }n,K,\rho,\nu.}
\]

相应histogram χ²=trG−1≤c_q/(1−c_q)。该结论必须使用Mult(n,μ)；原Mult(n,η)不具备这个Gram的常数特征向量与源端正交归一化。

## 4. Score倾斜的精确混合恒等式

令原查询律Q(dx)=a(x)dx。因为a∈[q/L,(q+1)/L]且q≥1，对任意实t都有

\[
m(t)=E_Qa^t\in(0,\infty),\qquad
Q_t(dx)=a(x)^tQ(dx)/m(t).
\]

固定t和任意非负可测F(x_1,...,x_n)，测度乘积直接给

\[
\boxed{E_{Q^n}\left[e^{t\sum_i\log a(X_i)}F(X^n)\right]
=m(t)^nE_{Q_t^n}F(X^n).}
\]

对可积有符号F也成立。该式不是独立性近似；独立性已经包含在Q^n与Q_t^n中。数值Gauss求积以dx为基准时，倾斜权重是w_j a(x_j)^(1+t)，主程序的 `(1+t)*log(a)` 因而正确。

为避免T_E的符号混淆，设μ_t=E_{Q_t}p，π_t=Mult(n,μ_t)：

- 若定义type向量的平方范数 `T_E²=Σ_{α∈E} g_α²/π_t(α)`，则E_{Q_t^n}T_E²≤trG_t≤1/(1−c_q)，从而混合矩≤m(t)^n/(1−c_q)。
- 若定义scalar event posterior `g_E=Σ_{α∈E}g_α`，用π_t(E)归一化，则E(g_E²/π_t(E))≤1，因为0≤g_E≤1且Eg_E=π_t(E)。这不是与上一个type平方和相同的量。
- 若T_E仍用原π_η归一化，倾斜恒等式仍精确，但不能直接套新prior trace帽。必须支付归一化比，下一节给出无法统一删除的严格障碍。

即使归一化正确，m(t)^n仍然存在；不同table的a_b、m_b(t)、Q_{b,t}、μ_{b,t}也通常不同。本结论没有把这些table置于共同查询律，没有控制score乘积与selector/共同尺度max的联合积分。

## 5. 原先验分母的指数付款：一般严格障碍

保持旧π_η为分母，记

\[
\widetilde{\mathcal T}_n
=\sum_\alpha\frac{E_{\nu^n}g_\alpha^2}{\pi_\eta(\alpha)}.
\]

由于Eg_α=π_μ(α)，Jensen与多项式定理给

\[
\widetilde{\mathcal T}_n
\ge\sum_\alpha\frac{\pi_\mu(\alpha)^2}{\pi_\eta(\alpha)}
=\left(\sum_c\frac{\mu_c^2}{\eta_c}\right)^n
=\left(1+\sum_c\frac{(\mu_c-\eta_c)^2}{\eta_c}\right)^n.
\]

只要μ≠η，底数严格大于1；旧分母trace就不可能有n无关上界。这不是仅仅缺少一个证明技巧，而是对“新prior帽直接支付原prior量”的严格反例机制。正确的新prior trace仍可≤2，两者没有矛盾。

简单精确例：q=1、ell=1/2、均匀源，两个class为半周期格。取查询Dirac在恰捕获class1的canonical残余窗中心，p≡(2/3,1/3)，η=(1/2,1/2)，μ=(2/3,1/3)。g_α≡π_μ(α)，所以新trace=1而旧trace=(10/9)^n。本机已对n=1,2,4,8,16,32,64精确展开核验。

## 6. 连续正来源、真实有限score倾斜的精确反例

旧先验障碍也能在实际连续查询倾斜中严格实现，不必仅依赖Dirac例。取两个class子密度

\[
\rho_1=\tfrac1{10}1_{[0,1/2]}+9\,1_{[9/40,11/40]},\qquad
\rho_2=1_{[1/2,1]}.
\]

周期延拓后η=(1/2,1/2)，总ρ≥1/10。使用真实中心窗全长度L=5/4、q=1、ell=1/4，直接定义

\[
z_c(x)=\frac1L\int_{x-L/2}^{x+L/2}\rho_c(y)dy,\quad a=z_1+z_2,\quad p_c=z_c/a.
\]

取t=1，m(1)=∫a²dx，μ_{1,c}=∫a z_c dx/m(1)。本审计独立按全长5/4的窗端点，使用源边界±5/8取模分成9段；每段a,z_c均仿射，积分二次乘积精确给

\[
m(1)=\frac{505697}{500000},\qquad
(\int az_1,\int az_2)=\left(\frac{128591}{250000},\frac{49703}{100000}\right),
\]
\[
\mu_1=\left(\frac{257182}{505697},\frac{248515}{505697}\right),\qquad
\sum_c\frac{\mu_{1,c}^2}{\eta_c}
=\frac{255804572698}{255729455809}>1.
\]

因此这个真实正来源在有限t=1下，原先验分母trace至少为该有理底数的n次方。直接全长积分没有沿用主脚本canonical ell窗的变量，故也独立核对了奇q半周期平移的积分不变量。数值扫描中其他有限t的μ偏移只是浮点探索证据，本节的有理例则严格证明旧归一化障碍。

## 7. 实际核验与范围

`checks.py`实际完成：144个任意正η/非负有限u与ν的Fraction帽；24个独立query枚举hist trace与矩阵幂trace的完全齐次公式；90个Fraction score倾斜恒等式；90个旧先验Jensen因子；7个指数toy展开；直接物理全长5/4例；9000个原law的只读重放及8个实际score-tilt归一化例。中心方差公式与主 `ΣEv²/μ−1` 的最大差约2.10e−15，解释了浮点负零量级。

已严格闭合的是任意查询律下**正确新先验归一化**的固定channel/histogram trace帽、混合倾斜恒等式及旧先验指数付款障碍。尚未闭合的是原来源/selector归一化下的实际score与多table共同积分；m(t)^n没有消失。以上均不提供C_n的新上界、下界或增长阶改善。

## 8. 原先验指数付款的精确谱定位

不仅能证明原先验分母量必然增长，还能把其全部n依赖定位为单个leading eigenvalue的n次方。令

\[
H_\eta=E_\nu[(p/\sqrt\eta)(p/\sqrt\eta)^T],\quad
u=\sqrt\eta,\quad d=\mu/\sqrt\eta,\quad
\chi=\sum_c\frac{(\mu_c-\eta_c)^2}{\eta_c},\quad \kappa=c_q.
\]

这里u是单位向量，H_ηu=d、u^TH_ηu=1、||d||²=1+χ。令γ=λ_max(H_η)。沿用V_c=p_c/η_c和ab=1，标量不等式V_c²≤(a+b)V_c−1给

\[
\operatorname{tr}H_\eta=\sum_c\eta_cE_\nu V_c^2
\le(a+b)\sum_c\eta_cE_\nu V_c-1
=a+b-1=1+\kappa.
\]

H_η≥0且H_η²≤γH_η，故

\[
1+\chi=\|H_\eta u\|^2\le\gamma u^TH_\eta u=\gamma.
\]

因此γ≥1，余下非负特征值之和≤κ。原先验归一化的histogram Gram仍有精确对称张量公式：把线性映射e_c↦p_c/√η_c张量化，normalized sym基的像正是g_α/√π_η(α)。这只是多项式/Gram恒等式，不要求η是ν下真实class边缘，也不声称该映射是Markov收缩。所以

\[
\widetilde{\mathcal T}_n=h_n(\operatorname{spec}H_\eta),
\quad
\boxed{\gamma^n\le\widetilde{\mathcal T}_n
\le\frac{\gamma^n}{1-\kappa}}.
\]

上界来自h_n(γ,λ_rest)=Σ_{j=0}^nγ^{n−j}h_j(λ_rest)≤γ^nΣ_j(κ/γ)^j≤γ^n/(1−κ)。q1给γ^n到2γ^n之间的夹逼。

进一步取正交分解R^K=span(u)⊕u^⊥，写

\[
H_\eta=\begin{pmatrix}1&\delta^T\\\delta&C\end{pmatrix}.
\]

C≥0且tr C=trH_η−1≤κ，||δ||²=||d−u||²=χ。因为γ≥1>κ≥||C||，leading eigenvector在u方向分量非零；由特征方程消元得到精确Schur恒等式

\[
\gamma-1=\delta^T(\gamma I-C)^{-1}\delta
\le\frac\chi{\gamma-\kappa}\le\frac\chi{1-\kappa}.
\]

χ=0时γ=1，同式仍成立。因此

\[
\boxed{1+\chi\le\gamma\le1+\frac\chi{1-\kappa}.}
\]

这定量定位了原先验付费，严格强于只报告“先验改变了”。对于score倾斜Q_t，第四节恒等式进一步给

\[
[m(t)\gamma_{\eta,t}]^n
\le E_{Q^n}\left[e^{t\sum\log a}\sum_\alpha\frac{g_\alpha^2}{\pi_\eta(\alpha)}\right]
\le\frac{[m(t)\gamma_{\eta,t}]^n}{1-\kappa}.
\]

这个原先验leading price不能被新先验trace≤2删除；对其与score的实际付款仍需要新估计。

第六节物理例的原先验Gram，直接全长积分z_cz_d/(η m)严格给

\[
H_\eta=\frac1{1517091}
\begin{pmatrix}822002&721090\\721090&770000\end{pmatrix}.
\]

其数值谱约(.04906916024,1.000308892886)。本审计独立精确重算矩阵，并对其二次特征多项式作Fraction符号核验，认证γ位于[1+χ,1+2χ]，χ=75116889/255729455809。144个有限law另检查了trace帽的Fraction方向及数值谱夹逼；谱小数本身只作数值说明。主 `price_probe.py` 的2×2完全齐次公式和log计算方向正确；百万维log价格是浮点描述，不充当C_n证据。
