# 第52轮证明：Fenchel–Ward恒等式与soft支持内的双符号障碍

日期：2026-11-07。N=2b，b≥4，c>0，A>1。
沿用第50—51轮的有界实际Φ_N预算及其KKT正则性条件。
以下恒等式和符号构造是严格结果；余项的小量估计尚未证明。

## 1. 正象限上的边界闭合

令R_b(t)=exp(t)f_{b,c}(t)，R(t₁,t₂)=R_b(t₁)R_b(t₂)，
dξ=exp(−t₁−t₂)dt₁dt₂，dν=R dξ。R_b在(0,∞)严格递增，由0到∞，
且exp(−t)R_b是Gamma–多项式密度。
设H=h(λR)为高度A内的Fenchel最优响应，P=∫Φ_N(H)dξ。
记Φ_{N,A}为Φ_N在[0,A]上的闭凸限制，a=λR，并置
\[
 F=\Phi_N(H)-\lambda RH=-\Phi_{N,A}^*(a).
\]
标量共轭函数是A-Lipschitz，故F局部绝对连续；几乎处处
∂_iF=−λH∂_iR。
a=1处H虽跳跃，F两侧均为0，没有额外delta；soft、平台与cap交界也连续。
因此弱意义下
\[
 \partial_i(e^{-t_1-t_2}F)=e^{-t_1-t_2}[-F-\lambda H\partial_iR].\tag{1}
\]

边界不能省略。这里−λAR≤F≤0，且R_b(0)=0。
边界截面满足
\[
 \int_0^\infty e^{-u-t_j}|F(u,t_j)|dt_j
 \le\lambda A f_{b,c}(u)\longrightarrow0
 \quad(u\downarrow0\text{ 或 }u\to\infty).
\]
又∂_iR≥0，∫∂_iR dξ=1，所以式(1)的导数绝对可积。
积分可得无需移动输入域的严格Ward恒等式
\[
 \boxed{\lambda\int H(1-k_i)\,d\nu=P,\qquad k_i=\partial_i\log R.}\tag{2}
\]
该证明处理了正象限的轴边界，不依赖假定整个域平移不变。

在D=[0,h]²上记
Γ_i=∫_{另一坐标∈[0,h]} exp(−h−t_j)F(h,t_j)dt_j≤0。
则准确关系是
\[
 P_D=\lambda\int_DH(1-k_i)d\nu-\Gamma_i,
 \qquad |\Gamma_i|\le\lambda A f_{b,c}(h).\tag{3}
\]
不能把有限方窗当成式(2)没有边界的全域。

## 2. 完整曲率的均值平移分解

记m(c)=1+c/(1+exp(c))，s_b=∂_c logR_b，
\[
 e_{b,c}(t)=\frac{\partial_c^2R_b(t)}{R_b(t)}-b m''(c)(1-k_b(t)),
\]
\[
 \mathcal E_N=\mathbb E_\nu H\{e_{b,c}(T_1)+e_{b,c}(T_2)
                       -2s_b(T_1)s_b(T_2)\}.
\]
沿(c+θ,c−θ)的对比路径，乘积求导和式(2)给
\[
 T_{--}=N m''(c)P/\lambda+\mathcal E_N.
\]
结合第50轮完整预算Hessian和交换对称性，严格得到
\[
 \boxed{(\log P)''_-=N[-m''(c)]-
       \frac\lambda P(\mathcal E_N+B_{--}).}\tag{4}
\]
B必须含下端coarea与全部soft响应，交叉score项也不能靠交换对称性删除。
在方窗上另有−b m''(c)(Γ₁+Γ₂)/P_D，大小至多
b|m''(c)|λA[f_b(h)+f_b(h)]/P_D。

μ=1+exp(−μ)，p_μ=(μ−1)/μ，有m″(μ)=−p_μ。
在c=bar c_N=μ−√(μ/N)z_N/2处，−m″(c)=p_μ+O(√(logN/N))。
于是
\[
 \frac\lambda P\bigl(|\mathcal E_N|+B_{--}\bigr)=o(N)\tag{5}
\]
是(logP)″_−∼Np_μ的明确充分条件。式(5)未闭合；
已有值渐近、下端局部展开或有限维数值均不能替代它。

## 3. 一次分部积分的精确接口

写p=(1+exp(c))⁻¹，f_i=f_b(t_i)，g_i=f_{b−1}(t_i−c)，
w_i=f_{b−2}(t_i−2c)，r_i=g_i/f_i，u_i=w_i/f_i，S=bp(r₂−r₁)。
第51轮参数导数给出
\[
 \frac{\partial_\theta^2(f_1f_2)}{f_1f_2}
 =bp\left[A_0+\frac{g'_1}{f_1}+\frac{g'_2}{f_2}\right],
\]
其中A₀=−2(1−p)+r₁+r₂+(b−1)p(u₁+u₂)−2bp r₁r₂。
设dh按分布/coarea意义包含H的下端跳跃和普通soft导数。
将g′项分部积分，得到
\[
 T_D+B_D=bp\int_D HA_0d\nu+
 \int_D a\,dh(a)\,\mathcal C\,d\nu+bp(E_{1,D}+E_{2,D}),\tag{6}
\]
\[
 \mathcal C=S^2-bp(k_1r_1+k_2r_2).
\]
例如E_{1,D}=∫f₂(t₂)[H(h,t₂)g₁(h)−H(0,t₂)g₁(0)]dt₂。
全正象限中cap保证顶边消失，底边因g_i(0)=0消失；有限窗必须保留E。

## 4. 两种符号均出现在实际soft支持

先由有限多项式得到两个严格渐近式：

- 若0<t₁<c且t₂=T→∞，则r₁=0，
  r₂∼(b−1)/(2pT)，k₂∼(b−1)/T，
  故C∼b(b−1)²(b−2)/(4T²)>0。
- 若0<t₂<c且t₁=c+ε↓c，则r₂=0，
  C∼−b(b−1)²c^(−b)ε^(b−2)<0。

仅在任意点找符号还不够：远尾可能已进入cap，dh为0。
下面对**每个固定λ>0、A>1**把两种符号放入普通soft集合N<a<Φ′_N(A)。

正号：固定任意soft值a*，令t₂=T，唯一解
λR_b(t₁(T))R_b(T)=a*。由R_b的正多项式两端渐近，
t₁(T)∼K/T<c，K>0，故仍有上面的正号，且a=a*始终在soft。

负号：若λR_b(c)²<Φ′_N(A)，取
max{N,λR_b(c)²}<a*<Φ′_N(A)，对角解t₁=t₂=t>c。
于是S=0、C=−2bpkr<0。
若λR_b(c)²≥Φ′_N(A)，固定任意soft值a*，取t₁=c+ε，
再唯一解出t₂<c使λR_b(t₁)R_b(t₂)=a*，即得到近结点负号。
这些点附近各有soft内的开集，故不是cap外或零权重处的形式符号。

因此式(6)中的soft局部系数不具任一固定符号；
仅对该系数作一次分部积分后宣称T+B有单向符号的路线被排除。
这没有决定完整积分Hessian的符号，因为第一项及其加权相消仍须结算。
