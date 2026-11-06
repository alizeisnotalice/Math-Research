# 第53轮证明：完整KKT响应在漂移尺度上为o(N)

日期：2026-11-07。这里的P是第50—52轮的**单接收、两个相等块**实际Orlicz最优预算，
不是一般强极大算子端点常数C_N。保留其唯一内点乘子和局部Hessian正则条件。

## 定理与范围

记f_{b,c}为b个Y_c=E+cD独立和的密度，R_{b,c}(t)=exp(t)f_{b,c}(t)。
令N=2b，L=logN，μ=1+exp(−μ)，固定K<∞，
|c−μ|≤K√(L/b)，高度上限为2，平均目标为1。
统一使用Φ_N(h)=h(1+log⁺h)^(N−1)。
在全正象限用dξ=exp(−t₁−t₂)dt，R=R_{b,c}(t₁)R_{b,c}(t₂)，
dν=R dξ，H=h(λR)，P=∫Φ_N(H)dξ，E_νH=1。
令s_b=∂_c logR_b，I_b(c)=E s_b(T)²，B为对比路径(c+θ,c−θ)的完整响应。
对固定K且N充分大，严格有
\[
 I_b(c)\le C_K(1+b(c-\mu)^2)^2=O_K(L^2),\qquad \lambda/P\le2N,\tag{1}
\]
\[
 B^{\rm soft}_{--}=O_K(L^2/N),\qquad
 B^{\rm low}_{--}=O_K(L^2/\sqrt N),\tag{2}
\]
因而
\[
 \boxed{\frac\lambda P B_{--}=O_K(\sqrt N\,L^2)=o(N).}\tag{3}
\]
常数不依赖具体漂移序列。参数窗口最终进入下文所用固定邻域；
不能把“充分大”误解为本轮所有数值维数都已进入显式阈值。

## 1. 紧参数、固定未标准化位置的Fourier引理

令Y_d=E+dD，E为单位指数变量，D为参数p_d=(1+exp(d))⁻¹的Bernoulli变量。
记m(d)=1+dp_d，v(d)=1+d²p_d(1−p_d)，g_{n,d}(x)为ΣY_d−nm(d)的密度。
对任意D₀⊂(0,∞)紧集及固定M，有一致估计
\[
 g_{n,d}(x)=\frac{1}{\sqrt{2\pi n v(d)}}[1+O_{D_0,M}(n^{-1})],
 \quad d\in D_0,\ |x|\le M.\tag{4}
\]

给出证明以免把第45轮μ邻域展开无依据地扩展到任意紧参数。
中心化特征函数ϕ_d满足
|ϕ_d(z)|≤(1+z²)^(−1/2)，因为Bernoulli因子的模至多1。
紧参数上各阶矩一致有界，v(d)≥1，且在固定小邻域有一致的对数展开
\[
 \log\varphi_d(z)=-v(d)z^2/2-i\kappa_3(d)z^3/6+O_{D_0}(z^4).
\]
Fourier反演令z=u/√n。在|u|≤n^(1/12)内，被积函数为
\[
 e^{-v(d)u^2/2}\left[1-
 \frac{i(xu+\kappa_3(d)u^3/6)}{\sqrt n}
 +O_{D_0,M}\left(\frac{u^2+u^4+u^6}{n}\right)\right].
\]
用稍弱的统一Gaussian包络控制余项；该截断内指数扰动一致趋零。
一阶项为奇的纯虚函数，在对称区间积分恰为0。
从n^(1/12)到小频率固定端点的积分为O(exp(−a n^(1/6)))。
外频利用Gamma因子，以固定可积幂分出指数衰减，故为O(exp(−a'n))。
Gaussian主项积分与余项的多项式Gaussian矩给出(4)，不是O(n^(−1/2))的相对误差。

## 2. 全状态score与漂移Fisher界

第45轮给出c在μ的某固定邻域J₀内的全状态界
\[
 -C[1+b(c-\mu)^2]\le s_{b,c}(t)\le bp_c.\tag{5}
\]
其精确倾斜桥公式为
\[
 s_{b,c}(t)=bp_c\{1-u(d)B_b(d)/u(c)\},\quad
 t/b=c\,m(d)/d,\quad u(d)=d/(1+e^d),
\]
其中B_b(d)=g_{b−1,d}(m(d)−d)/g_{b,d}(0)。
式(4)给任意紧d区间上一致B_b(d)=1+O(1/b)。
u在μ处有二阶严格极大；c固定邻域内u(c)有正下界。
对t/b处在包含所有m(c)的固定紧区间，逆关系的导数
\[
 \frac{d}{dd}\{c m(d)/d\}=-c/d^2-cp_d(1-p_d)
\]
一致远离0，所以|d−c|≤C|t/b−m(c)|。
Taylor估计遂给
\[
 |s_{b,c}(t)|\le C\left[1+b(c-\mu)^2+
                         \frac{(t-bm(c))^2}{b}\right].\tag{6}
\]
紧区间外，最后一项至少为常数乘b，而(5)给|s|≤Cb，故(6)仍成立。
独立和的四阶中心矩为3b²v(c)²+bκ₄(c)=O(b²)，常数在J₀上一致。
平方(6)并积分得(1)中的Fisher界。这里只使用固定K漂移窗口，
没有从c=μ处的Fisher极限直接推断漂移结论。

## 3. 乘子界及显式有限阈值

对c∈J=[6/5,27/20]和b≥40000，证明λ/P≤2N。
此阈值版只需c∈J；漂移窗口最终进入J∩J₀后即可用于渐近定理。

以下常数均可由24项正Exp级数及几何尾严格比较，代码使用精确有理运算：
\[
 e^{27/20}<4,\quad e^{13/10}>10/3,\quad
 (6/5)e^{129/100}>e^{6/5}+1.
\]
因此m(c)>31/25，而m(c)的全局最大为m(μ)=μ<13/10；
第三式给1+log(6/5)−log(1+exp(−6/5))>91/100。

k=0 Gamma项及单调性给
\[
 \log R_b(t)\ge\log R_b(6b/5)
 \ge .91b-\log(6b/5)-1\ge .90b\quad(t\ge6b/5).\tag{7}
\]
这里logΓ(b)=Σ_{j=1}^{b−1}logj≤∫₁ᵇlogx dx=b logb−b+1，
无需未控制的Stirling余项。
在b=40000，log48000<16且.01b−16−1>0；此后差的导数.01−1/b>0。

Var(Y_c)<3/2，由Chebyshev及并界，
ν{T₁,T₂≥6b/5}≥1−1875/b>.95。
若λ exp(.90N)≥Φ′_N(11/10)，该集合上H≥1.1，违反E_νH=1。
故
\[
 \lambda<\Phi'_N(11/10)e^{-.90N}.\tag{8}
\]

由log2>.693、log(11/10)<.096及exp(.43)<1693/1096，严格有
α=log[(1+log2)/(1+log(11/10))]>.43。
若H=2，则(8)给
logR>.90N+.43(N−2)≥1.32N（N≥86）。
一步密度在两段上的最大值均≤1；卷积也满足f_b≤1。
所以logR≤T₁+T₂，cap概率由总和上尾控制。

在0≤η≤.01的指数倾斜下，
Λ″_c(η)≤(100/99)²+(27/20)²/4<3/2。
取η=.01得Λ_c(η)≤1.3η+.75η²，因此
\[
 \nu\{H=2\}\le e^{-N/8000}\le e^{-10}<1/4.\tag{9}
\]
最后，非cap处aH≤NΦ_N(H)，cap处贡献2λν{H=2}，从而
λ(1−2ν{H=2})≤NP。式(9)给所需λ/P≤2N。

## 4. 全部soft响应

令S=s_b(T₁)−s_b(T₂)。对1≤H≤2及N≥8，
\[
 \frac{\Phi'_N(H)}{\Phi''_N(H)}=
 \frac{H(1+\log H)(N+\log H)}{(N-1)(N-1+\log H)}
 \le\frac{4(N+1)}{(N-1)^2}\le\frac6N.
\]
soft上a=Φ′_N(H)，故
\[
 B^{\rm soft}_{--}\le\frac6N E_\nu S^2
                      =\frac{12}{N}I_b(c).
\]
用了块独立及E s_b=0；不需要限制在局部下端带。
乘上(1)的乘子界，得到(λ/P)Bsoft≤24I_b=O_K(L²)。

## 5. 对任意KKT层一致的下端coarea界

U=logR_b(T)严格单调，令q_b为其推前密度：
\[
 q_b(u)=f_b(t(u))/k_b(t(u)),\qquad k_b=\partial_t\log R_b.
\]
正截断多项式逐项给k_b(t)≥(b−1)/t。
记γ_b为单位率Gamma密度，Gamma–二项混合及xγ_b(x)=bγ_{b+1}(x)给
\[
 \sup_t t f_b(t)\le\sup_x x\gamma_b(x)+cb\sup_x\gamma_b(x)=O_J(\sqrt b).
\]
这里supγ_b≤2/√b（b≥8）可直接由Fourier绝对积分得到：
|z|≤1用(1+z²)^(−b/2)≤exp(−bz²/4)，外频用
2^(−b/4)(1+z²)^(−2)。于是
\[
 \|q_b\|_\infty\le C_J/\sqrt b\quad
 (c\in J\text{ 时可取 }10/\sqrt b).\tag{10}
\]

在a=1处δ(a−1)=δ(loga)，没有额外λ因子。
记w=−logλ、A_b(u)=s_b(t(u))。coarea等于卷积层密度：
\[
 B^{\rm low}_{--}=\int[A_b(u)-A_b(w-u)]^2q_b(u)q_b(w-u)du
 \le4\|(q_b)\|_\infty I_b(c).
\]
这是对所有w一致的界，不要求定位最优下端层或控制特定角向区域。
由(1)、(10)得Blow=O_K(L²/√N)，再乘λ/P≤2N即得(3)。

## 6. 尚未闭合的形状项

第52轮Ward分解现变为
\[
 (\log P)''_-=N[-m''(c)]-(\lambda/P)\mathcal E_N+o(N).
\]
本轮只消除了完整非负响应B的量级障碍；
仍须证明(λ/P)|E_N|=o(N)才可断言对比曲率∼Np_μ。
以上定理没有改善一般C_N的对数级上下界。
