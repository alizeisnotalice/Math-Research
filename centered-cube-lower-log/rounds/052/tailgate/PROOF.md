# 单侧计数尾门：matched-amplitude 混合集中到一个合法产品

本轮先运行真实连续周期核的数值，再独立证明下述结构。结论覆盖合法仿射softsplit的任意单侧count-tail（包括极稀尾）；一般非负countweights先给出精确混合恒等式，只有满足计数方差O(n)的子类进入集中结论。没有使用自由Gaussian通道。

## 1. 模型和所有计数门的精确混合

取质量1的一维密度ρ、π∈(0,1)，并假设

\[
\rho_1=\pi\rho+c(\rho-1)\ge0,\qquad
\rho_0=(1-\pi)\rho-c(\rho-1)\ge0.
\]

其质量分别为π、1−π。n个独立source-label对的标签总数记K，先验为Bin(n,π)。F_k是条件K=k后积分标签的质量1输入密度，G_{k,L}是它在同一L下的张量平均。若w(k)≥0、0<Z_w=Ew(K)<∞固定于source端，则

\[
F_w=\sum_{k=0}^n\nu_kF_k,\quad
\nu_k=\frac{w(k)\binom nk\pi^k(1-\pi)^{n-k}}{Z_w},\qquad
G_{w,L}=\sum_k\nu_kG_{k,L}.
\tag{1}
\]

所有尺度共用同一ν；禁止按查询或L重新选权。单侧尾是w=1_{K≤k₀}或1_{K≥k₀}，其Z_w可以极小。

记Y_L=A_Lρ、d_L=Y_L−1，则严格有

\[
G_{w,L}(x)=\prod_iY_L(x_i)\,
\frac{E_{p_i}[w(K)]}{E_\pi[w(K)]},\quad
p_i=\pi+c\frac{d_L(x_i)}{Y_L(x_i)}.
\tag{2}
\]

数值里的PB递推及tail CDF就是(2)，分母始终保留。

## 2. 对全部k统一的matched产品domination

置τ=c/[π(1−π)]、z=k/n、α_z=1+τ(z−π)。归一化指数source倾斜给合法密度

\[
\rho_z=1+\alpha_z(\rho-1),\quad
Y_{z,L}=1+\alpha_zd_L,\quad T_{z,L}=\prod_iY_{z,L}(x_i).
\]

它来自ρ₀与ρ₁的正组合；z0/1为端点密度。假设来源序列在指定尺度族满足|d_L|≤h_n→0、|τ_n|≤T<∞。于是|α_z|≤A=1+T、m_n=1−Ah_n>0最终成立。

51轮的统一单计数比较是

\[
G_{k,L}\le \gamma_nT_{k/n,L}\quad\text{对全部k,L,x},\qquad \gamma_n\to1.
\tag{3}
\]

为使本稿自含，重述论证。指数倾斜后

\[
G_{k,L}/T_{z,L}=P_{q_i}(K=k)/\operatorname{Bin}(n,z)\{k\},\quad
q_i=z+z(1-z)e_i,\quad |e_i|\le b_n:=Th_n/m_n\to0.
\]

令η_n=b_n+b_n²/4，则Σq_i(1−q_i)≥(1−η_n)nz(1−z)。PB的Fourier反演以
$I(V)=(2\pi)^{-1}\int_{-\pi}^{\pi}e^{-2V\sin^2(t/2)}dt$
上界点质量。变量替换u=√Vt与支配收敛给I(V)√(2πV)→1；Stirling公式对min(k,n−k)→∞给中心Bin点质量乘√(2πnz(1−z))→1，且误差由最小阶乘参数一致控制。故在min(k,n−k)大的区域，比值≤(1−η_n)^−1/2(1+o(1))。

小计数端直接逐标签向量比较，点质量比≤exp(2b_n min(k,n−k))。若b_n=0则所有后验和匹配产品已相同；否则取M_n=min(√n,b_n^−1/2)，则M_n→∞、b_nM_n→0；在min(k,n−k)≤M_n用后界，余下区域用前界。端点k0/n严格相等。这证明(3)，不使用tail质量的倒数帽。

由(1)与非负性立刻得到

\[
G_{w,L}\le \gamma_n\underbrace{\sum_k\nu_kT_{k/n,L}}_{T_{\nu,L}}.
\tag{4}
\]

这是逐点凸组合，不是weak-L¹函数和三角不等式；此时还没有降到单个产品。

## 3. 单侧截断Binomial的方差≤n/4：自含证明

先证Binomial CDF序列F_j=P(Bin(m,π)≤j) logconcave。其pmf f_j的相邻比递减。对内部j，令r=f_{j+1}/f_j，则

\[
rF_{j-1}=\sum_{i<j}rf_i\le\sum_{i<j}f_{i+1}\le F_j.
\]

所以F_j f_j≥F_{j−1}f_{j+1}，即F_j²≥F_{j−1}F_{j+1}；范围外以CDF=0或1延拓，边界同样成立。

现条件化n个iid Bernoulli标签于K≤k₀。n≥2时，取剩余n−2标签的CDF F，令D=P(K≤k₀)>0。任意两标签的条件2×2概率为

\[
p_{11}=\pi^2F_{k_0-2}/D,\quad
p_{10}=p_{01}=\pi(1-\pi)F_{k_0-1}/D,\quad
p_{00}=(1-\pi)^2F_{k_0}/D.
\]

于是

\[
\operatorname{Cov}(C_1,C_2\mid K\le k_0)
=p_{11}p_{00}-p_{10}p_{01}
=\frac{\pi^2(1-\pi)^2}{D^2}
[F_{k_0-2}F_{k_0}-F_{k_0-1}^2]\le0.
\]

条件标签仍交换。若其单个均值为q，则

\[
\boxed{\operatorname{Var}(K\mid K\le k_0)
=nq(1-q)+n(n-1)\operatorname{Cov}(C_1,C_2\mid E)\le n/4.}
\tag{5}
\]

上尾K≥k₀对补标签n−K应用同一证明；n1直接成立。结论对任意π、阈值及极稀正质量尾统一，不要求高斯或中偏差近似。

## 4. matched幅度混合的L¹集中

设ν为尾门先验，令\bar z=E_νK/n，β=α_{\bar z}=E_να_{K/n}，并取同一个合法产品T_{β,L}。β介于端点幅度之间，因此ρ_β合法，且不依赖L或查询。

固定L，记d=d_L、Y_α=1+αd、T_α=Y_α^{⊗n}。对整段合法幅度，Y_α≥m_n。利用∫d=0、乘积律与Cauchy，

\[
\int|\partial_\alpha T_\alpha|
=E_{Y_\alpha^{\otimes n}}\left|\sum_i\frac{d(X_i)}{Y_\alpha(X_i)}\right|
\le\sqrt{n\int d^2/Y_\alpha}\le\sqrt n\,h_n/\sqrt{m_n}.
\]

每项score的均值是∫d=0，故没有额外n倍漂移。对α积分，再对ν平均，得一般方差版本

\[
\|T_{\nu,L}-T_{\beta,L}\|_1
\le\frac{|\tau|h_n}{\sqrt{m_n}}\sqrt{\operatorname{Var}_\nu(K)/n}.
\tag{6}
\]

特别由(5)，单侧tail满足

\[
\boxed{\|T_{\nu,L}-T_{\beta,L}\|_1
\le |\tau|h_n/(2\sqrt{m_n}).}
\tag{7}
\]

所以同样结论适用于任何固定非负countweights，只要其归一化计数先验满足Var_ν(K)=O(n)。任意权并不自动满足：双端点等权的方差n²/4，uniform-count prior为n(n+2)/12，均被本次数值保留。

## 5. 真正支付阈值split：有限增长网格结论

取同一有限尺度族Λ、J=|Λ|，令

\[
Z=\max_{L\in\Lambda}|T_{\nu,L}-T_{\beta,L}|,
\quad\|Z\|_1\le J|\tau|h_n/(2\sqrt{m_n}).
\]

由(4)，$\max_LG_{w,L}\le \gamma_n(\max_LT_{\beta,L}+Z)$。这里是**函数和**，不能沿用第50轮max/max的同阈值union。对任意θ∈(0,1)必须分阈值，严格得到

\[
\boxed{
\|\max_LG_{w,L}\|_{1,\infty}
\le\frac{\gamma_n}{1-\theta}\|\max_LT_{\beta,L}\|_{1,\infty}
+\frac{\gamma_nJ|\tau|h_n}{2\theta\sqrt{m_n}}.}
\tag{8}
\]

在接受源族局部尺度网格，h_n=O(1/K)、J=O(Kv)、|τ|有统一界、v→∞、K/v^s→∞。取θ=1/v，第二项O(v²)=o(√n)，第一系数$\gamma_n/(1−1/v)=1+o(1)$。因此

\[
\boxed{\|\max_\Lambda G_{w,L}\|_{1,\infty}
\le(1+o(1))\|\max_\Lambda T_{\beta,L}\|_{1,\infty}+o(\sqrt n).}
\tag{9}
\]

它把这一单侧count-tail机制降到**一个合法且尺度无关的matched产品**，不提供该产品本身的弱型上界。没有把未知速率的$(\gamma_n−1)$乘节点数J后偷当o(√n)，故本稿不擅自把(9)改写为系数恰1的加性比较。

本集中路线不覆盖Var_ν(K)为n²量级的任意非负权；它们仍只有(4)的幅度混合结构，需另用root的一般prior比较。也不把有限Λ的L¹求和扩成连续所有尺度，不处理任意空间macroclass或非仿射softsplit。

## 6. 本轮真实数值的范围

连续周期来源ρ=1+.9cos、π.5、c.05的softsplit最小值为ρ₁≥.005、ρ₀≥.095，合法。实际区间平均严格是Y_L=1+.9 sinc(L)cos(2πx)，因此n维立方体平均可用真正张量核计算；不是Gaussian代理。该固定cos来源的h不随n→0，所以不能仅凭n16/64/128实验调用(9)的接受源族渐近前提。

固定query seed、同3尺度、先验不随L或query改变；训练4096查询选阈值及单产品幅度，另4096独立查询核对。PB DP保留n128、k0的原质量2^−128，约2.94e−39（log质量−88.723），并用完整归一化后输出。h-grid步长.125，无法分辨极小profile差；经验最佳值不作真实weaknorm证书。

`verify_exact.py`另用Fraction检查1824个单双侧尾先验，最大Var(K|tail)/n恰1/4；没有删除极稀尾。最小交付两份源码、冻结计划、两份结果及本证明；NOTE为本地说明，不需重复归档。没有C_n新log下界，也没有一般C_n上界。
