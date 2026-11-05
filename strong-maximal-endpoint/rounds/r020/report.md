# 第20轮：一般上界的对数平方余项

日期：2026-11-06（Asia/Shanghai）。状态：一般输入上界。使用有限时域Gaussian离散屏障；本轮定理不依赖鞅中偏差误差。

## 1. 主要结论

令L=log n，mu=log B，sigma=sqrt(mu)，barPhi_G(z_n)=sigma z_n/sqrt n，d_n=sigma sqrt n z_n。Phi_G为标准正态CDF，Phi_n仍为本项目指定的Orlicz模数。
存在固定K和n_0使n≥n_0时
\[
\boxed{\log C_n\le\log\mathcal R_n\le n\mu-d_n+K(\log n)^2.}\tag{1}
\]
第19轮实际输入下界给
\[
 n\mu-d_n+\frac{9-\mu}{24}\log n+o(\log n)
 \le\log C_n\le n\mu-d_n+O((\log n)^2).\tag{2}
\]
因此log C_n=nmu−d_n+O((log n)²)。这将第14轮上侧O(n^(3/8)(log n)^(1/8))大幅缩小，但仍不确定一般log n系数或一般C_n的乘法渐近等价。
以下先对全部可预测策略给全状态CDF下界，再由第11轮已证Q_n≤C_n≤mathcalR_n转为一般输入结论。

## 2. 一步均值、方差及统一矩

一步Z_r≥0的密度满足f_r(z)≤2e^-z，对r∈[0,1]一致。设m_r、v_r为均值及方差，D=mu−m_r。前轮已由显式公式及唯一非退化均值最大点证明
\[
 0\le D\le D_{\max}:=\mu-1,\quad |v_r-\mu|\le C_m\sqrt D.\tag{3}
\]
因密度界，对固定eta_0=1/8，E|Z_r−mu|^j e^(eta_0|Z_r−mu|)在j≤5上一致有界。
参考增量G~Normal(mu,mu)有同样的有限矩性质。这里不假设自适应总和的无条件方差为O(n)。

## 3. Gaussian一步比较引理

记g_tau(x)=Phi_G(x/(sigma sqrt tau))，f_tau=g_tau'。
取固定A=10，L≥1，tau≥C_0 L，且|x|/(sigma sqrt tau)≤2A sqrt L。
当C_0足够大时，存在与C_0的进一步增大无关的固定C_*，使所有r满足
\[
 E g_\tau(x+\mu-Z_r)\ge g_{\tau+1}(x)-C_*\frac L\tau f_\tau(x).\tag{4}
\]
证明如下。以Y=mu−Z_r及Y_0=mu−G作三阶Taylor比较，Gaussian半群给E g_tau(x+Y_0)=g_(tau+1)(x)，故差为
\[
 D f_\tau(x)+\tfrac12(v_r+D^2-\mu)f_\tau'(x)+\mathrm{Rem}.
\]
Gaussian导数及精确密度比给，对z=x/(sigma sqrt tau)，
\[
 \frac{|g_\tau'''(x+u)|}{f_\tau(x)}
 \le\frac C\tau(1+z^2+u^2/\tau)
 e^{|z||u|/(\sigma\sqrt\tau)}.
\]
选C_0使2A/(sigma sqrt C_0)≤eta_0，第2节的加权矩控制Taylor积分余项，得到
|Rem|≤C(1+z²)f_tau(x)/tau≤C' L f_tau(x)/tau。
这一估计允许无界r参数及无界Z；没有截断增量，也没有把全局绝对误差当作相对密度误差。

令R=|f_tau'/f_tau|=|z|/(sigma sqrt tau)。再使D_max R≤1，则(3)使主项除以f_tau至少为
\[
 D-\tfrac12D^2R-\tfrac12 C_m\sqrt D R
 \ge\tfrac12D-\tfrac12 C_m\sqrt D R
 \ge-C_m^2R^2/8.
\]
最后一个不等式由完成平方精确得到。R²≤4A²L/(mu tau)，于是(4)成立。
常数选择次序：先固定A、eta_0和统一矩界，再固定C_*；所有充分大的C_0共享同一个C_*。

## 4. 有限时域平滑与右移密度的下界

固定大整数N，令L=log N、tau_0=C_0 L、epsilon=N^-4、delta=N^-5。
取K≥8C_*+8，然后继续增大固定C_0。定义
\[
 H_0=4\sigma\sqrt{\tau_0 L},\quad
 h_k=\frac{K L}{k+\tau_0},\quad H_k=H_{k-1}+h_k,
\]
\[
 W_k(q)=\Phi_{\rm G}\left(\frac{q-\mu k-H_k}{\sigma\sqrt{k+\tau_0}}\right),
 \qquad V_k(q)=W_k(q)-\epsilon-k\delta.
\]
H_N≤C_H L²，C_H固定；tau_0只为平滑初值，不对应真实额外增量。

在当前中央带|z|≤A sqrt L、z=(q−mu k−H_k)/(sigma sqrt(k+tau_0))，置tau=k−1+tau_0、x=q−mu k−H_(k−1)。则|x|/(sigma sqrt tau)≤2A sqrt L，(4)可用。
Gaussian参考一步给g_(tau+1)(x)=W_k(q+h_k)。必须定量估计右移收益，不能只说密度可比。
对0≤u≤h_k，精确公式给
\[
\log\frac{f_{\tau+1}(x-u)}{f_\tau(x)}
=-\tfrac12\log(1+1/\tau)+\frac{x^2}{2\mu\tau(\tau+1)}
 +\frac{2xu-u^2}{2\mu(\tau+1)}.
\]
其下界为−1/(2tau)−|x|h_k/(mu(tau+1))−h_k²/(2mu(tau+1))。
由于|x|≤2sigma A sqrt(tau L)、h_k=K L/(tau+1)，这三个负项分别不超过
1/(2C_0)、2AK/(sigma C_0^(3/2))、K²/(2mu C_0³)。
继续增大C_0使其和≤log2，就有
\[
 W_k(q+h_k)-W_k(q)
 =\int_{x-h_k}^{x}f_{\tau+1}(v)dv
 \ge\tfrac12h_k f_\tau(x)
 \ge\frac K4\frac L\tau f_\tau(x).\tag{5}
\]
中央带映到较宽的2A带；右移区间由显式密度比处理，无需假设移后仍在原带内。
由(4)—(5)及K的选择，中央带中E_r W_(k−1)(q−Z_r)≥W_k(q)。

## 5. 两个外带与初值

左外带z<−A sqrt L时，W_k(q)≤N^(-A²/2)≤delta，而E W_(k−1)≥0，所以允许误差delta的一步界成立。
右外带z>A sqrt L时，引入独立标准正态U，
\[
 E W_{k-1}(q-Z_r)
 =P\{Z_r+\sigma\sqrt\tau U\le\mu+h_k+\sigma\sqrt{\tau+1}z\}.
\]
其补事件由Z_r>mu+h_k+sigma sqrt(tau+1)z/2，或U>sqrt((tau+1)/tau)z/2覆盖。密度尾界给
\[
 1-E W_{k-1}(q-Z_r)
 \le2e^{-\sigma\sqrt{\tau+1}z/2}+e^{-z^2/8}
 \le2N^{-\sigma A\sqrt{C_0}/2}+N^{-A^2/8}\le\delta
\]
对充分大N成立，增大C_0使第一个指数至少8；A=10使第二个指数12.5。
故全实轴、全部r、1≤k≤N均有
\[
 E_r W_{k-1}(q-Z_r)\ge W_k(q)-\delta.\tag{6}
\]

初值方面，q<0时W_0(q)≤Phi_G(−4sqrt L)≤N^-8/2≤epsilon；q≥0时W_0≤1。
所以V_0(q)≤1_{q≥0}，且(6)给E_r V_(k−1)(q−Z_r)≥V_k(q)。
逐步条件化允许动作依赖完整过去，并不要求策略为固定或Markov策略。
从终端V_0下界反向迭代，得到对每个可预测策略pi及所有实q，
\[
\boxed{P_\pi(S_N\le q)\ge
 \Phi_{\rm G}\left(\frac{q-N\mu-H_N}{\sigma\sqrt{N+\tau_0}}\right)-2N^{-4}.}\tag{7}
\]
V可以为负，这不影响其作为终端指标函数下界及条件期望迭代的合法性。

## 6. 全状态分布界给精确到对数平方的超额期望

令s=Nmu−d，其中d=d_N+O(L²)。取T=(mu+1)N。
对(7)的上尾积分至T，超过T的部分用前轮一致正mgf界
log E exp(eta S_N)≤Nmu eta+2VN eta²，V=1+4/e²，eta=1/(4V)，得到
\[
 E(S_N-s)_+\le E(G_N-s)_++O(N^{-3})+O(e^{-cN}),\tag{8}
\]
其中G_N~Normal(Nmu+H_N,mu(N+tau_0))。这是分布比较的推论，不是对原自适应方差作假设。
设a=d+H_N~d_N，Gaussian负尾积分不超过var(G_N)/a乘exp(−a²/(2var(G_N)))=O(sqrt L)。
因此
\[
\boxed{E_\pi(S_N-s)_+\le d+C_E L^2.}\tag{9}
\]
(7)中的小常数尾误差只积分到T，不能直接把它积分到无穷；(8)已经单独控制剩余尾。

## 7. 近阈值输入账本需要同步锐化

沿用第11轮mathcalR_N定义：非负有界非减h，J=∫e^-t Phi_N(h(t))dt，存在pi使E_pi h(S_N)>1。
输入模数仍为Phi_N(t)=t(1+log^+t)^(N−1)；正态CDF记为Phi_G。
假设J≤e^-s。单调性给Phi_N(h(t))≤e^(t−s)，故h≤N^-2于t≤s−2L，h≤1于t≤s。
取中段长度b_N=sqrt N L²，反解模数精确得到
\[
 (h(t)-1)_+\le(1+O(b_N/N))(t-s)/N\quad(s\le t\le s+b_N).
\]
远尾按h<e^(4t/N)与其补集分开。低振幅部分用theta=c b_N/N的正mgf Chernoff界，贡献≤exp(−c'L^4)；高振幅部分用条件密度联合界，贡献≤J 2^N(N/4)^(N−1)/(N−1)!≤e^-c''N。
这里b_N/d→∞、b_N/N→0，theta最终≤1/2；常数c取足够小即可。
结合严格可行性与(9)，基准1分解给
\[
\boxed{P_\pi(S_N\le s-2L)\le d/N+O(L^{5/2}/N).}\tag{10}
\]
误差来自(1+O(b_N/N))(d+C_E L²)/N，以及可忽略远尾。不能继续使用第14轮的粗负尾代理超额估计：它不足以提供本尺度误差。

## 8. 最后比较与一般输入转移

选择固定K_0>C_H+2，令s=Nmu−d_N+K_0 L²、d=d_N−K_0 L²。
在q=s−2L处，(7)正态左尾的距离至多为d_N−L²，且其方差比mu N更大。
故以beta=d_N/N、f_G=phi(z_N)/(sigma sqrt N)~L/N表示，有
\[
 P_\pi(S_N\le s-2L)
 \ge\overline\Phi_{\rm G}\left(z_N-\frac{L^2}{\sigma\sqrt N}\right)-2N^{-4}
 \ge\beta+cL^3/N.\tag{11}
\]
最后一步直接积分标准正态密度；位移趋0，密度至少phi(z_N)。
而(10)给≤beta+O(L^(5/2)/N)，与(11)矛盾。
因此不存在J≤e^-s的可行轮廓，log mathcalR_N≤s，得到(1)。第11轮已证的一般输入重排桥完成C_N的转移。

## 9. 本机软件诊断与证明角色

MATLAB R2026a使用解析Exp卷积
\[
 E_Y\Phi_{\rm G}((a-Y)/b)=\Phi_{\rm G}(a/b)-e^{-a+b^2/2}\Phi_{\rm G}(a/b-b)
\]
以及稳定erfcx等价式，对一步实际混合核直接检查残差。
诊断参数C_0=10^4、K=100，N=10^4、10^6、10^8，k=1、floor(N/2)、N，标准化状态0及±{1,3,6}sqrt(log N)。每点在c∈[0.05,5]做fminbnd并补查0、50、∞端点。
按W_k'归一化的最小采样残差分别为0.00793774135、0.00107160478、0.0000162565662，均为正。
这些参数没有被认证为本轮存在性证明中的足够大常数；采样不替代第3—5节的全参数证明。
Wolfram15符号核验Exp卷积ODE、均值—方差完成平方及Gaussian密度比恒等式，并以80位独立积分复算卷积。
复现：MATLAB `addpath('code'); sm_gaussian_barrier_checks;`；Wolfram `wolframscript -file code/gaussian_barrier_checks.wl`。
三项符号残差均为0，三次80位积分均吻合至45位；MATLAB Code Analyzer无报告。确定性算法，无随机种子；浮点值未作定向认证。本轮没有声称有限维的最优常数证书。

## 10. 剩余差距与下一轮

一般上下界的对数差距现为O((log n)²)，比此前显著缩小，但仍大于O(1)。真正log n系数、乘法等价及强矩形极值形态仍未闭合。
下一轮两项动作：① 将全带统一O(L/tau)移位换成状态相关的曲率补偿，研究能否把上侧总损失降为O(log n)；② 用平台轮廓的完整凸振幅结构检验可实现的对数级系数，并保持实际输入与反馈上界模型的区别。
本轮额外上界估计均在报告中直接证明；外部基础来源沿用sources.md。无新增优先权主张，无冗余已跟踪材料需清理。
