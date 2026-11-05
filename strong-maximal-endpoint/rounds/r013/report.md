# 第13轮：一般输入的完整三项展开

日期：2026-11-06（Asia/Shanghai）。状态：严格渐近定理；有限维Bellman结果为数值诊断。

## 1. 结论与范围

沿用Phi_n(t)=t(1+log^+t)^(n−1)、mu=log B=1+B^-1、sigma=sqrt(mu)，令L=log n。则
\[
\boxed{\log C_n=n\mu-\sigma\sqrt{nL}
 +\sigma[\log L+\log(\sigma\sqrt{2\pi})]\sqrt{n/L}
 +o(\sqrt{n/L}).}\tag{1}
\]
同一展开成立于第11轮的一般反馈上界mathcalR_n。实际输入下界来自第8轮Q_n，故并未把上界模型冒充真实强矩形算子。
本轮闭合一般第三项。尚未给出C_n的乘法渐近等价、有限维最优常数或极值形态，也未修复附件整个双曲示性族H_n的另一条证明。

## 2. 单步核及统一条件

第11—12轮已严格建立Q_n≤C_n≤mathcalR_n。mathcalR_n是所有有界非负非减h及可预测策略pi中J^-1的上确界，其中
\[
 J=\int_0^\infty e^{-t}\Phi_n(h(t))dt,\qquad E_\pi h(S_n)>1.
\]
S_n=ΣZ_i，每步在观察过去之后选择r_i∈[0,1]，其条件分布为Z_r=Y−(log r)J_r，Y~Exp(1)、J_r~Bern(r/(1+r))独立；r=0按Exp极限。
记m(r)、v(r)为均值和方差，r*=B^-1。第12轮的精确公式及紧性给固定C_m：
\[
1\le m(r)\le\mu<2,\quad v(r)\ge1,\quad v(r)\le V:=1+4/e^2,
\quad |v(r)-\mu|\le C_m\sqrt{\mu-m(r)}.\tag{2}
\]
核密度≤2e^-z；所以对整数k≥3，条件于任意过去，
\[
 E|Z_r-m(r)|^k\le2^{k-1}(E Z_r^k+m(r)^k)
 \le2^{k+1}k!\le\tfrac12 k!32^{k-2}v(r).\tag{3}
\]
其中EZ^k≤2k!、m^k≤2^k≤2k!，最后比率为16^(3−k)≤1。
同时正向条件中心mgf满足log E exp(eta(Z_r−m(r)))≤2V eta²，0≤eta≤1/2。

## 3. 停止损失，保留全部鞅增量

固定任意原策略pi，令D_i=Σ_{j≤i}(mu−m(r_j))。它递增，每步增量≤mu−1。
取A=2sigma+8sqrt(V)、D*=A sqrt(nL)。构造新策略piK：遵循pi直到首次D_i>D*的该步结束，此后所有步使用r*。
动作仍然可预测；停止后不是把增量置零。新策略路径上
\[
 D_n^K\le D^*+\mu-1.
\]
令M_K=Σ(Z_i−m(r_i^K))。由Cauchy–Schwarz，路径逐点有
\[
 |\langle M_K\rangle_n-\mu n|
 \le C_m\sum_i\sqrt{\mu-m(r_i^K)}
 \le C_m\sqrt{n(D^*+\mu-1)}=O(n^{3/4}L^{1/4}).\tag{4}
\]
这一界是可预测二次变差的几乎处处界，绝非无条件Var S_n的替代。

## 4. 相对鞅中偏差的合法应用

使用Fan--Grama--Liu，arXiv:1210.2198v3，第2节Theorems 2.1—2.2、式(14)。条件A1为
|E(xi_i^k|F_{i−1})|≤(k!/2)epsilon^(k−2)E(xi_i²|F_{i−1})，k≥3；A2为|<X>_n−1|≤delta² a.s.。
对xi_i=−(Z_i−m(r_i^K))/(sigma sqrt n)，(3)给epsilon_n=32/(sigma sqrt n)，(4)给delta_n²=C_v n^-1/4 L^1/4。常数独立于原策略。
原文结论是相对尾而非绝对Berry–Esseen误差，其对数误差受以下和控制：
\[
 x^3\epsilon_n+x^2\delta_n^2+(1+x)(\epsilon_n|\log\epsilon_n|+\delta_n).
\]
对0≤x≤2sqrt L，此和为O(kappa_n)，其中
\[
 \kappa_n=n^{-1/8}L^{5/8}=o(1).\tag{5}
\]
原文范围x≤alpha0/epsilon_n、delta_n≤alpha0对充分大n一致满足；注意delta_n是A2右端的平方根。
因此
\[
 P_K(M_K\le-d)=(1+O(\kappa_n))\overline\Phi(d/(\sigma\sqrt n)),
 \quad0\le d\le2\sigma\sqrt{nL}.\tag{6}
\]
终端分布无原子：最后一步条件核有密度，因此原文严格尾与此处闭尾无差别。

原始来源：[arXiv v3全文](https://arxiv.org/html/1210.2198v3)，其条件、版本与适用尺度逐项登记于sources.md。本轮并不声称文献首创性。

## 5. 从新策略回到任意原策略

两策略用相同随机数耦合，停止之前完全相同。在未停止事件上，若M_K≤−d，则S_n=nmu−D_n+M_K≤nmu−d，原策略也失败。
停止事件上，若原策略成功S_n>nmu−d，则原中心鞅M_pi>D*−d≥8sqrt(VnL)。
对u=8sqrt(VnL)，迭代正向中心mgf并取eta=u/(4Vn)≤1/2，得到P_pi(M_pi>u)≤exp(−u²/(8Vn))=n^-8。
故
\[
\begin{aligned}
 P_\pi(S_n\le n\mu-d)
 &\ge P_K(M_K\le-d)-P_\pi(\text{停止且成功})\\
 &\ge(1-O(\kappa_n))\overline\Phi(d/(\sigma\sqrt n))-n^{-8}.\tag{7}
\end{aligned}
\]
第一行可由“P_K(M_K≤−d)≤P_pi(失败且未停止)+P_pi(停止)”再减去原失败中的停止部分得到；在耦合下停止事件概率相同。
窗口最小正态尾为asymp n^-2/sqrt L，故n^-8可吸收到相对O(kappa_n)内。
这是所有策略失败概率的统一下界；没有声称每个原策略的S_n近似正态。固定r*的失败概率由同一文献以delta=0得到相应的相对正态等式。

## 6. 改善轮廓超额的首系数

设s=nmu−d、J≤e^-s，其中d~sigma sqrt(nL)。单调性严格给Phi_n(h(t))≤e^(t−s)。
令b=n^(3/4)，则t≤s−2L时h≤n^-2，t≤s时h≤1。
对s≤t≤s+b，设y=log h(t)>0；由
(n−1)y/(1+y)≤logPhi_n(h(t))≤t−s
反解y≤(t−s)/(n−1−(t−s))=O(n^-1/4)。因此一致有
\[
 (h(t)-1)_+\le(1+O(n^{-1/4}))(t-s)/n.\tag{8}
\]
写原策略S_n=nmu−D_n+M_pi，D_n≥0，E|M_pi|≤sqrt(nV)。所以E(S_n−s)_+≤d+sqrt(nV)。

远尾t>s+b分成h(t)<exp(4t/n)与相反区域。
前者用theta=c n^-1/4、eta=theta+4/n的正向mgf：
E[e^(4S_n/n)1_{S_n>s+b}]≤exp[−theta(b−d)+2Vn(theta+4/n)²+4mu]≤exp(−c1 sqrt n)，取固定足够小c>0。
后者由联合密度≤2^n exp(−Σz_i)，从而f_S(t)≤2^n e^-t t^(n−1)/(n−1)!，再用Phi_n(h)≥h(4t/n)^(n−1)，得到期望≤J 2^n(n/4)^(n−1)/(n−1)!≤e^-c2n。
这两个估计一致于策略和轮廓。
按基准1分解期望并使用严格可行性，得到
\[
 P_\pi(S_n\le s-2L)
 \le\frac{(1+O(n^{-1/4}))(d+\sqrt{nV})/n+e^{-c_1\sqrt n}+e^{-c_2n}}{1-n^{-2}}
 =(\sigma+o(1))\sqrt{L/n}.\tag{9}
\]

## 7. 比较Mills前因子，确定一般第三项

任意固定实数C，取
\[
 s=n\mu-\sigma\sqrt{nL}+\sigma\log L\sqrt{n/L}+C\sqrt{n/L},
 \qquad d_A=n\mu-s+2L.
\]
则x=d_A/(sigma sqrt n)满足x²=L−2log L−2C/sigma+o(1)。Mills公式和(7)给
\[
 P_\pi(S_n\le s-2L)
 \ge\left(\frac{e^{C/\sigma}}{\sqrt{2\pi}}+o(1)\right)\sqrt{L/n}.\tag{10}
\]
若C>sigma log(sigma sqrt(2pi))，(10)与(9)矛盾。于是所有充分大n都没有J≤e^-s的可行轮廓。
任意固定C高于该常数均被排除，证明一般上界的limsup。第8轮已构造的有界截断指数软轮廓给同一系数的liminf，并由固定维数深层极限实现为实际可积输入。因此Q_n≤C_n≤mathcalR_n夹逼证明(1)。

## 8. 本机计算、测试族与下一步

实际MATLAB R2026a Bellman递推扩展到1024步；步长0.05，≤256步的高精度诊断沿用第11轮0.0125网格。取q=nmu−sigma sqrt(n log n)，比较最小失败的网格近似与正态尾：

|n|反馈失败/正态尾|固定r*失败/正态尾|
|--:|--:|--:|
|64|0.664891316|0.751201033|
|128|0.705991920|0.778898159|
|256|0.747637357|0.808607580|
|512|0.788280355|0.837704294|
|1024|0.824543315|0.864698421|

反馈有限维数优于固定核；该差别不改变本轮三项展开。状态插值误差未认证，表格不是有限维常数或定理阈值的证书。
Wolfram15实际复核三个FGL误差项极限均为0、条件矩比率16^(3−k)、停止成功尾的log n指数8。复现：MATLAB `sm_bellman(.05,1024,fullfile(pwd,'rounds','r013','bellman_1024.json')); sm_controlled_gaussian_checks;`；Wolfram `wolframscript -file code/stopping_checks.wl`。
FGL存在性常数未数值化，不能声称1024已进入定理的充分大n范围。

下一轮两项动作：① 将正态尾方程barPhi(d_n/(sigma sqrt n))=d_n/n作为隐式中心，保留统一相对误差以改善余项；② 用本机符号展开及高精度求根独立核对后续逆对数项。一般乘法渐近等价与真实极值结构仍未解决。
