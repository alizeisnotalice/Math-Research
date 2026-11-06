# 第44轮：反馈模型的锐对数项及一般输入上界

日期：2026-11-07。令mu=1+e^(−mu)、B=e^mu、L=log N、sigma=sqrt(mu)，以barPhi_G(z_N)=sigma z_N/sqrt(N)定义d_N=sigma sqrt(N)z_N。记

c_Q=(9−mu)/24，g=(mu−1)/8，c_R=c_Q+g=(3+mu)/12。

## 1. 结论与范围

反馈模型R_N定义沿用r011：对非负有界非减h及可预测策略pi，在E_pi h(S_N)>1时取真实预算J(h)=∫₀∞e^(−t)Phi_N(h(t))dt之倒数上确界。Phi_N(t)=t(1+log⁺t)^(N−1)。本轮证明

log R_N=Nmu−d_N+c_R log N+o(log N)。                         (1)

结合已证明的Q_N≤C_N≤R_N与r022的Q_N展开，有

c_Q≤liminf (log C_N−Nmu+d_N)/log N
   ≤limsup (log C_N−Nmu+d_N)/log N≤c_R。                     (2)

不声称C_N的该极限存在，不把R_N等号转成C_N等号。模型相对独立比例机制的分离现为精确对数阶：R_N/Q_N=N^(g+o(1))。

## 2. 匹配参考的均值、方差及隐式根

取任意固定有界C，令s=Nmu−d_N+C L，匹配a=a_N(C)满足ell(a)=2a−m(a)=s/N。记h=a−mu、M=Nm(a)、V=Nnu(a)、K₃=N kappa₃(a)。在mu处

m(mu)=nu(mu)=mu，m'(mu)=0，m''(mu)=−p_mu，
nu'(mu)=mu−1，kappa₃(mu)=mu(3−mu)，p_mu=(mu−1)/mu。

隐函数定理及Taylor给h=−d_N/(2N)+O(L/N)，因此

M=Nmu−gL+o(L)，
V=Nmu−(mu−1)d_N/2+O(L)，K₃=Nmu(3−mu)+O(sqrt(NL))。          (3)

更精确地均值损失是g z_N²+o(L)，而z_N²/L→1。所有余项对固定C紧集一致。

对V附近定义z_V、d(V)=sqrt(V)z_V，barPhi_G(z_V)=sqrt(V)z_V/N。r022的隐式微分恒等式给

d(V)=d_N+[d_N/(2mu N)](V−mu N)+o(L)=d_N−2gL+o(L)。           (4)

其中O(L)的方差余项只造成o(L)的d(V)余项。令y=M−s，则

y=d_N−(C+g)L+o(L)，y−d(V)=(g−C)L+o(L)。                     (5)

## 3. 有限终端参考payoff的统一展开

记f_N、F_N为该固定a的和密度及CDF，x=y/sqrt(V)。r022自包含Fourier证明的一项绝对Edgeworth展开为

F_N(s)=barPhi_G(x)+[K₃/(6V^(3/2))](1−x²)phi_G(x)+O(1/N)。     (6)

其常数在a趋mu的紧邻域统一。局部密度估计同时给

f_N(s)=[phi_G(z_V)/sqrt(V)](1+o(1))∼L/N。                    (7)

说明：(5)使x−z_V=O(L/sqrt(N))，Gaussian密度比趋1；绝对O(1/N)密度误差相对L/N趋零。Gaussian Mills比和隐式根使phi_G(z_V)/sqrt(V)∼z_V²/N∼L/N。

在(6)的Gaussian项围绕z_V展开，用barPhi_G(z_V)=d(V)/N及(5)，得到

barPhi_G(x)−y/N=(C−g)f_N(s)L+o(f_N(s)L)。

根右侧差(d(V)−y)/N=O(L/N)，比f_N L≈L²/N小一阶。Edgeworth项为

−[(3−mu)/6]f_N(s)L+o(f_N(s)L)，

因为K₃/(6V)→(3−mu)/6、x²/L→1。因此

F_N(s)−y/N=[C−g−(3−mu)/6]f_N(s)L+o(f_N(s)L)。                (8)

这里保留了均值损失g和方差补偿2g，不能只从log Q_N的最终系数推测参考payoff。

实际终端phi_N(r)=e^(−r)（r>0）、1−r/N（r≤0）。其参考值精确为

v_N(s)=1+y/N−F_N(s)+R_N^exp(s)+E(s−S_N)_+/N，
R_N^exp(s)=E[e^(S_N−s)1_(S_N<s)]。                            (9)

此处R_N^exp是低侧加权项，与反馈常数R_N不同。Gamma结构给R_N^exp≤(s/N)f_N=O(f_N)。共同负倾斜lambda∼sqrt(L/N)/sigma及全左移密度比给短缺量≤C f_N/lambda²=O(1)。故(9)最后两项分别O(f_N)、O(1/N)，均为o(f_N L)。结合(8)与g+(3−mu)/6=c_Q，得到核心参考式

v_N^a(s)=1−(C−c_Q)f_N(s)L+o(f_N(s)L)。                      (10)

固定h₀的状态平移只改变O(f_N)：v_N'=−R_N^exp−P(S_N>s)/N，第二项必须保留，但由(7)它相对f_N仅为O(1/L)。因此(10)也适用于s−h₀，误差仍o(f_N L)。这一步未假定最优轮廓或最优根被达到。

## 4. 全局势的端点值

使用r043已经证明的同一有限N全空间势与整个屏障平移初始化：对固定alpha>1，存在固定h₀、Lambda及有界源系数s_N，使

mathcalU_(N,N)(s)≤W_(N,N)(s−h₀)≤I_(N,N)(s−h₀)。             (11)

第二个不等式只在当前临界窗口使用：该窗口在正侧强制区之外。不能声称W≤I在所有状态成立。

在精确匹配点s=Nell，Poisson公式的二次部分为0，故

J_N(s)/f_N(s)=(p_a/8)(C_N(s)/N)(H_N−1)。

辅助密度中心Fourier比给C_N/N=Lambda_a''(−lambda)+O(1/N)→mu；pa→p_mu。所以

J_N(s)/f_N(s)=gL+o(L)。                                      (12)

实际上差J_N/f_N−gL趋g(gamma_E−1)，但锐项仅需(12)。固定平移使J变化o(f_N L)。其余代价逐项为

s_N f_N=O(f_N)，Lambda b²=O(L/N)=O(f_N)，固定平移代价O(f_N)。

这里b=m(a)−a≈d_N/(2N)。由(10)有v_N→1，U=alpha J+s_Nf_N=O(f_N L)→0，因此指数化只多出O(U²)=o(f_N L)。由(10)–(12)，

mathcalU_(N,N)(s)
 ≤1−[C−c_Q−alpha g]f_N(s)L+o(f_N(s)L)。                     (13)

给定epsilon>0，先选满足全局构造小参数条件的固定alpha>1，使(alpha−1)g<epsilon/2，再取C=c_Q+g+epsilon。所有屏障常数包括h₀均先于N固定；(13)及(7)给

mathcalU_(N,N)(s)≤1−c_epsilon L²/N，c_epsilon>0，             (14)

对充分大N成立。不能令alpha依赖N后沿用这些常数。

## 5. 实际预算、反馈上界及模型等号

r032已经对全部允许轮廓和全部可预测策略统一证明：若J(h)≤e^(−s)、s=Nmu−d_N+O(L)，则

E_pi h(S_N)≤E_pi phi_N(s−S_N)+O(L/N)
             ≤mathcalU_(N,N)(s)+O(L/N)。                    (15)

其高层尾用一个实际Orlicz积分预算的Fenchel吸收，未重复花费预算。由(14)，(15)右侧严格小于1，排除J≤e^(−s)的所有可行轮廓。因此log R_N≤s。这里只需一个临界预算；若上确界超过e^s，按上确界定义即有预算更小的可行h，不需最优者存在，也不需证明非临界全局根唯一。

epsilon任意，得到limsup(log R_N−Nmu+d_N)/L≤c_Q+g=c_R。r029的真实可预测策略及实际有界软轮廓分母已经给匹配liminf≥c_R，完成(1)。r011的Q_N≤C_N≤R_N给(2)。

真实C_N仍有宽度g的未闭合对数系数区间；(1)说明仅使用这个反馈重排上界模型无法把上界系数降到c_Q，不能据此断言真实矩形极大算子也取得反馈收益。
