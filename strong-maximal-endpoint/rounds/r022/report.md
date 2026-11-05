# 第22轮：独立比例机制的精确对数项

日期：2026-11-06（Asia/Shanghai）。状态：Q_n的严格渐近及一般上界初始化的范围障碍；不声称一般C_n的最优对数系数已确定。

## 1. 结论与机制范围

沿用mu=log B、sigma=sqrt(mu)、L=log n，以barPhi_G(z_n)=sigma z_n/sqrt n定义d_n=sigma sqrt n z_n。则
\[
\boxed{\log Q_n=n\mu-d_n+\frac{9-\mu}{24}\log n+o(\log n).}\tag{1}
\]
Q_n允许任意确定性比例向量(r_1,...,r_n)及任意非负非减标量轮廓；坐标增量独立。本轮不是把比例限制为相等，也没有把轮廓限制为有限阶梯或平台。
第19轮实际输入给(1)的下界。本轮证明匹配上界，关键是保留均值缺口D及阈值下方的精确指数加权项。
一般输入仍只有
\[
 n\mu-d_n+\frac{9-\mu}{24}L+o(L)\le\log C_n\le n\mu-d_n+O(L).
\]
比例依赖过去的反馈过程不满足本轮独立Fourier分解，故不能用(1)替代C_n或mathcalR_n的一般上界。

## 2. 统一结构与输入预算

对每个r_i，写Z_i=Y_i+B_i，Y_i~Exp(1)，B_i非负且为0或c_i=−log r_i；P(B_i=c_i)=r_i/(1+r_i)。所有变量独立；r_i=0时定义B_i=0，避免无穷乘零。
令S=ΣZ_i、M=Σm_i=nmu−D、V=Σv_i、K_3=Σκ3_i。
前轮的显式均值、方差及唯一非退化均值极大点给
\[
 n\le V\le V_{\max}n,\quad
 |V-\mu n|+|K_3-n\mu(3-\mu)|\le C\sqrt{nD}.\tag{2}
\]
第三累积量的稳定性与方差相同：在极大点附近用Taylor，在其余紧集用正均值缺口；端点连续，因为c^j e^-c→0。
各阶固定矩一致有界，且每个特征函数均含Exp因子的模(1+t²)^-1/2。

对固定实数C，测试s=nmu−d_n+C L，假设存在轮廓h及比例向量使
\[
 J=\int_0^\infty e^{-t}\Phi_n(h(t))dt\le e^{-s},\qquad E h(S)>1.
\]
h非减严格给J≥e^-t Phi_n(h(t))。所以t<s时h(t)≤e^(t−s)，t≥s时反解模数。
令b_n=sqrt n L，u=t−s∈[0,b_n]，y=log h(t)>0。由y≤2u/n及
u≥ny−(n−1)y²/2，得到一致二阶界
\[
 (h(t)-1)_+\le u/n+C u^2/n^2.\tag{3}
\]
远尾t>s+b_n沿用密度与正mgf分解，低振幅贡献≤e^-cL²，高振幅贡献≤e^-c'n，均为o(L/n)。所有常数一致于比例和轮廓。
记d=nmu−s、y_0=M−s=d−D。初步有
E(S−s)_+≤(y_0)_++sqrt V≤d+sqrt V，且E(S−s)_+²≤d²+V=O(nL)，后式利用S−s≤d+(S−M)。
因此基准1分解及h≤n^-2于t≤s−2L给粗约束
\[
 F_S(s-2L)=O(\sqrt{L/n}).\tag{4}
\]

## 3. 粗定位先排除大均值损失

若D≥2d，Chebyshev使F_S(s−2L)≥1−O(1/L)，与(4)矛盾。故D≤2d=O(sqrt(nL))，由(2)得V/(mu n)→1。
独立非同分布Berry–Esseen误差一致为O(n^-1/2)，因此
\[
 F_S(s-2L)=\overline\Phi_G((y_0+2L)/\sqrt V)+O(n^{-1/2}).
\]
若该距离小于(1−epsilon)sqrt L，正态尾将比sqrt(L/n)大一个幂次，与(4)矛盾；另一方面y_0≤d~sigma sqrt(nL)。所以
\[
 y_0\sim\sigma\sqrt{nL},\qquad D=o(\sqrt{nL}).\tag{5}
\]
这些结论一致于所有本测试中可行的比例和轮廓。特别地s<M、s>0，且充分大n时s<nmu。

## 4. 倾斜局部密度与精确Gamma低侧项

取唯一lambda>0使负倾斜e^-lambda S后的均值等于s。Exp部分的均值下降至少n lambda/(1+lambda)，故lambda=O(y_0/n)→0；均值Taylor再给
lambda~y_0/V~sqrt(L/n)/sigma。
倾斜后各Exp速率为1+lambda，Bernoulli参数为(1+exp((1+lambda)c_i))^-1。
它们的矩一致有界、总方差asymp n。对倾斜居中和密度g，Exp特征函数的乘积衰减及三阶局部Taylor直接给
\[
 g(0)\asymp n^{-1/2},\quad\|g\|_\infty=O(n^{-1/2}),\quad\|g'\|_\infty=O(n^{-1}).\tag{6}
\]
证明与第16轮相同但逐项累加：小频区log特征函数余项为O(n|t|³)，积分给中心密度绝对误差O(n^-1)；外区抽取Exp因子的一半幂给指数衰减。无需坐标同分布。
负倾斜公式与(6)给
\[
 f_S(s)/F_S(s)=\lambda(1+O(L^{-1/2})),\quad
 f_S(s-u)\le C f_S(s)e^{-\lambda u}\quad(u\ge0).\tag{7}
\]
第一式将∫_{−∞}^0e^(lambda v)g(v)dv与g(0)/lambda比较，相对误差≤C/(lambda sqrt n)；第二式只用全局g上界和中心下界。

低侧还可用更强的精确结构。S=Gamma(n,1)+B，B=ΣB_i≥0独立，故对s>0，
\[
\begin{aligned}
 f_S(s)&=e^{-s}E[e^B(s-B)_+^{n-1}]/\Gamma(n),\\
 R(s):=E[e^{S-s}1_{S<s}]&=e^{-s}E[e^B(s-B)_+^n]/\Gamma(n+1).
\end{aligned}
\]
于是
\[
\boxed{R(s)\le(s/n)f_S(s).}\tag{8}
\]
不存在B>s的额外项，因为两边相应幂均为0。
由Eh>1、(3)和粗高侧预算，F_S(s)≤R(s)+O(sqrt(L/n))。结合(7)—(8)，lambda→0，吸收R≤C lambda F后得到F_S(s)=O(sqrt(L/n))。
再由(7)积分，
\[
 E(s-S)_+=\int_0^\infty u f_S(s-u)du
 \le C f_S(s)/\lambda^2=O(1).
\]
所以E(S−s)_+=y_0+O(1)，R=O(L/n)。将其代回真实逆模数预算(3)，最终保留均值缺口得到
\[
\boxed{F_S(s)\le y_0/n+O(L/n).}\tag{9}
\]
直接丢掉D、或使用sqrt V的粗超额误差，都不足以确定本轮系数。

## 5. 三角阵列的一项绝对Fourier展开

令X_i=Z_i−m_i。对所有比例向量，全实轴上一致有
\[
 F_S(M-\sqrt V x)=\overline\Phi_G(x)
 +\frac{K_3}{6V^{3/2}}(1-x^2)\phi(x)+O(n^{-1}).\tag{10}
\]
这里不要求近最优比例；近最优性仅用于(2)的累积量识别。
完整误差论证如下。每坐标中心四阶矩一致有界，Vasymp n，标准化和的特征函数满足
\[
 \log\chi(t)=-t^2/2+K_3(it)^3/(6V^{3/2})+O(t^4/n)
\]
于|t|≤epsilon sqrt n。由纯虚三阶项及|e^(iv)−1−iv|≤v²/2，特征函数减Gaussian修正项的误差≤C(t⁴+t⁶)e^-ct²/n。
外区使用
|chi(t)|≤(1+t²/V)^(-n/2)≤(1+t²/(Vmax n))^(-n/2)，其除以|t|后的积分指数小。
全线误差除以|t|可积且O(1/n)，CDF Fourier反演给(10)。Gaussian修正项可为有符号密度，但总质量1，反演合法；与第18轮同一证明不依赖同分布。

## 6. 隐式根反演与D=O(L)

对任意V∈[n,Vmax n]，定义z_V>0、d(V)=sqrt V z_V，使
barPhi_G(z_V)=sqrt V z_V/n。则z_V²~L一致，且
\[
 d'(V)=\frac{z_V}{2\sqrt V}
 \left(1+\frac{\sqrt V}{n\phi(z_V)}\right)^{-1}
 =O(\sqrt{L/n}).\tag{11}
\]
由(2)、(5)，a_n:=K_3/(6V)→a:=(3−mu)/6一致。
对任意固定epsilon>0，在y=d(V)−(a+epsilon)L处代入(10)。正态尾位移与三阶修正的差给
\[
 F_S(M-y)-y/n\ge c_\epsilon L^2/n,
\]
因为sqrt V归一化Gaussian密度asymp L/n，a_n z_V²=aL+o(L)，而(10)绝对误差及根右侧位移只有O(1/n)、O(L/n)。
函数F_S(M−y)−y/n严格递减；(9)迫使实际y_0大于该测试点。因此
\[
 s\le M-d(V)+aL+o(L).\tag{12}
\]
(11)与方差稳定性给|d(V)−d_n|≤C sqrt(DL)。将s=nmu−d_n+C L代入(12)，有
D≤C_1L+C_2 sqrt(DL)，故
\[
\boxed{D=O(L).}\tag{13}
\]
这一步先排除较大D，再做锐方差展开，不循环假设近最优参数已经局部化。

现在V−mu n=O(sqrt(nL))。在该区间(11)给
\[
 d(V)=d_n+\tfrac12\lambda_n(V-\mu n)+o(L),\quad
 \lambda_n=d_n/(\mu n).
\tag{14}
\]
理由：d'(V)/(lambda_n/2)=1+O(1/L+|V−mu n|/n)；积分误差为O(1)+o(1)，特别为o(L)。

## 7. 每坐标的全局方差补偿

定义H_lambda(r)=mu−m(r)+lambda(v(r)−mu)/2。lambda→0时，它一致趋向非负均值缺口，只有r*=B^-1处为0；所以全局极小点进入r*的固定邻域。
在c=−log r变量内m''(mu)=−(mu−1)/mu、v'(mu)=mu−1。隐函数与Taylor严格给
\[
 \min_rH_\lambda(r)=-\frac{\mu(\mu-1)}8\lambda^2+O(\lambda^3),\quad
 c_{\min}=\mu-\mu\lambda/2+O(\lambda^2).\tag{15}
\]
远离该邻域的固定正均值损失排除了其他全局极小点；端点已包含在r∈[0,1]的紧性中。
由(15)，
\[
 D+\tfrac12\lambda_n(V-\mu n)
 \ge-\frac{\mu-1}{8}z_n^2-Cn\lambda_n^3
 =-\frac{\mu-1}{8}L+o(L).\tag{16}
\]
把(14)、(16)代入(12)，得到s≤nmu−d_n+[(3−mu)/6+(mu−1)/8]L+o(L)。
任取测试C>(9−mu)/24便矛盾，故log Q_n≤nmu−d_n+(9−mu)L/24+o(L)。
第19轮给同系数下界，完成(1)。全部误差一致于任意比例和可行轮廓；不需要Q_n的上确界被某个轮廓达到。

## 8. 一般上界初始化的机制障碍

对常数移位Gaussian初值W_0(q)=Phi_G((q−b)/(sigma sqrt t_0))，若q<0时必须≤epsilon=N^-p，则连续性给
b≥sigma sqrt t_0 Phi_G^-1(1−epsilon)~sigma sqrt(2p t_0 log N)。
在第21轮使用的全带加权矩控制条件t_0≥c log N下，这种初始化必须付出Omega(log N)的有效位移。非线性初值也必须把其在q=0处的总有效位移计入；仅锐化局部K不能自动确定最终log N系数。
这是当前Gaussian初始化方法的范围障碍，不是所有可能屏障都不可能改进的定理。

正支撑的Gamma(2,1) CDF G_2(q)=1−e^-q(1+q)（q≥0；q<0取0）是有效一步低基准：G_2(q)≤tanh(q/2)，由1+q−2/(1+e^-q)在0为0且导数≥1/2得到。
但它不能直接接普通Gaussian热半群：q<0时E G_2(q−Z_r)=0，而G~Normal(mu,mu)给E G_2(q−G)>0。
本机独立积分在q=−1给正确值0.001895398038066558…。
Gamma参考半群或精确短时Bellman拼接仍需新的比较证明，不能从该正支撑初值单独推出锐系数。

## 9. 本机软件、证据与下一轮

MATLAB R2026a实际优化H_lambda，lambda=0.2…0.001，c搜索区间[0,20]、TolX=1e-10。minH/lambda²趋向−mu(mu−1)/8，lambda=0.001时约−0.04450415636。
对维数2,4,8,16,32,64,128,256，使用两组不同确定比例的精确Gamma–二项双混合，检查CDF与(10)、以及精确R/f≤s/n。测试标准化位置0、sqrt(log n)，不使用Monte Carlo。
例如x=sqrt L时，n乘CDF误差在n=8,32,128,256分别约−0.04290、−0.02134、−0.000744、0.006212；这些采样不是统一误差证书。
Wolfram15独立符号核验(15)的两个系数与Gamma加权低侧恒等式；80位积分复核第8节反例，和MATLAB吻合。
复现：MATLAB `addpath('code'); sm_independent_log_checks;`；Wolfram `wolframscript -file code/independent_log_checks.wl`。
MATLAB Code Analyzer无报告。确定性算法，无随机种子，双精度值未区间认证；严格结论由第2—7节证明。

下一轮两项动作：① 研究一般上界在正支撑参考核下的初始化或短时拼接，明确能否去掉非锐的同阶位移；② 检验能超越独立比例机制的新矩形/相关输入，而不是重复已闭合Q_n族内搜索。
一般C_n的最优log n系数及乘法等价仍未解决。基础来源沿用sources.md及前轮已核验文献；本轮新增估计已直接证明。无新增优先权主张，无冗余跟踪材料需清理。
