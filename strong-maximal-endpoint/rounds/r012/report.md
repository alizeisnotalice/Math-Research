# 第 12 轮：一般端点常数的首个次指数项

日期：2026-11-06（Asia/Shanghai）。状态：一般输入定理。证明使用第11轮全局重排归约、可预测条件倾斜和似然比鞅，不假设自适应总和具有独立和方差。

## 1. 主要结论

令B满足log B=1+B^-1，mu=log B，sigma=sqrt(mu)。存在绝对常数K和N，使n≥N时
\[
 \boxed{\log C_n\le n\mu-\sigma\sqrt{n\log n}+K\sqrt n.}\tag{1}
\]
结合第3轮的显式可积输入下界（或第8轮更强的软轮廓下界），得到
\[
 \boxed{\log C_n=n\mu-\sigma\sqrt{n\log n}+O(\sqrt n),\qquad
 \lim_{n\to\infty}\frac{n\mu-\log C_n}{\sqrt{n\log n}}=\sigma.}\tag{2}
\]
特别地C_n/B^n→0，因此C_n与B^n不具有双边常数可比性。
这确定了一般输入的首个次指数缺口系数；并不声称一般C_n已具有第8轮受限Q_n的完整三项展开，更不声称渐近等价。

## 2. 一般输入到可预测过程的桥

第11轮已证明，而非假设，
\[
 Q_n\le C_n\le\mathcal R_n.
\]
具体地，对有界紧支撑f，沿每个坐标的rising-sun恒等式给
sD≤∫_0^D f*+∫(f−s)_+。其一维偶递减实现g#(y)=f*(2|y|)达到此分布上界，因此
\[
 (M_i f)^*\le T f^*,\quad
 (Tg)(x)=\sup_{0\le r\le1}\frac{\int_0^xg+\int_0^{rx}g}{(1+r)x}.
\]
T保序并保持递减性，所以(M_Sf)*≤T^n f*。
中间坐标函数无需全局L1：有界盒支撑种子被H乘以各坐标区间极大函数的乘积支配，正水平集有限；下一未使用坐标仍有紧支撑，故纤维恒等式合法。
一般Phi_n可积f由min(f,m)1_[-m,m]^n递增逼近，严格halo和模数单调收敛。

置x=e^-u，T变为sup_r E h(u+Z_r)。有限时域反向归纳允许r_i为过去Z的Borel可预测函数；各步条件分布为Z_(r_i)。
对递减g的任意严格接收点x，缩放为g(xy)并截去y>1不改变从1出发的n步终端值，因为所有乘子e^-Z≤1。
这说明mathcalR_n可定义为：对有界非负非减h:[0,∞)→[0,∞)，
\[
 J=\int_0^\infty e^{-t}\Phi_n(h(t))dt>0,
 \quad\exists\pi:E_\pi h(S_n)>1,
 \qquad\mathcal R_n=\sup J^{-1}.
\]
无界轮廓可先截断，严格余量由单调收敛保留。本轮只需对这一定义给统一上界。

## 3. 单步参数与一致常数

写c=−log r，r∈(0,1]。Z_r=Y+cJ0，其中Y~Exp(1)、J0~Bernoulli(r/(1+r))独立；r=0按Exp(1)极限。
均值m(r)=1−r log r/(1+r)，方差v(r)=1+r(log r)^2/(1+r)^2。
均值在唯一r*=B^-1最大，m(r*)=v(r*)=mu；m''(r*)<0。
因此由局部Taylor和远处紧性，存在有限常数C_m使
\[
 |v(r)-\mu|\le C_m\sqrt{\mu-m(r)}\quad(0\le r\le1).\tag{3}
\]

对单步作向下倾斜e^-lambda Z，lambda≥0。令p_lambda=(1+e^(c(1+lambda)))^-1，则
\[
 m_r(\lambda)=\frac1{1+\lambda}+cp_\lambda,\qquad
 v_r(\lambda)=\frac1{(1+\lambda)^2}+c^2p_\lambda(1-p_\lambda).
\]
在所有r、lambda≥0上一致有
\[
 v_r(\lambda)\le V:=1+4/e^2,\qquad
 |v_r'(\lambda)|\le L_3:=2+27/e^3.\tag{4}
\]
因c^k e^-c的最大值为k^k/e^k；r=0的连续极限同样成立。
另外在0≤eta≤1/2，正向log mgf二阶导数≤4V，故
\[
 \log E e^{\eta Z_r}\le\mu\eta+2V\eta^2.\tag{5}
\]
这些都是单步一致界，未假设S_n的无条件方差为O(n)。

## 4. 逐步变测度的精确成本

固定小delta>0，目标均值t_delta=mu−delta。
在某步所选r满足m(r)≤t_delta时保持原律；否则用唯一lambda>0将倾斜均值降到t_delta。
令Delta=m(r)−t_delta∈(0,delta]。
倾斜均值严格下降至0，所以根存在唯一；且
\[
 \Delta=\int_0^\lambda v_r(a)da\ge\lambda/(1+\lambda),
 \quad\lambda\le\Delta/(1-\Delta)\le2\delta
\]
（delta<1/2）。由(3)—(4)，在[0,lambda]上
\[
 v_r(a)\ge\mu-\varepsilon_\delta,
 \qquad\varepsilon_\delta=C_m\sqrt\delta+2L_3\delta.
\]
取delta足够小使该下界正。
相对熵成本k(lambda)=−lambda m_r(lambda)−log E e^-lambda Z满足k'=lambda v_r(lambda)、k(0)=0，故
\[
 k(\lambda)=\int_0^\Delta\lambda(b)db
 \le\frac{\Delta^2}{2(\mu-\varepsilon_\delta)}
 \le\frac{\delta^2}{2(\mu-\varepsilon_\delta)}.\tag{6}
\]
未倾斜步成本为0，也满足同界。这里没有丢失自适应性：变测度时保留同一可预测决策函数r_i(z_1,...,z_(i−1))，仅更换当前条件增量核。

## 5. 路径似然比给统一低尾下界

令ell=log n，任取0≤d≤2sigma sqrt(n ell)，并设
\[
 q=n\mu-d-2\ell,\quad\omega=4\sqrt V,\quad
 \delta=(d+2\ell+\omega\sqrt n)/n.
\]
delta一致趋0，故充分大n时上述单步构造均合法。
原路径律为P，新条件倾斜路径律为Q。所有条件均值m_(Q,i)≤mu−delta，所有条件方差≤V。
所以Q下鞅M_Q=Σ(Z_i−m_(Q,i))满足E_Q M_Q²≤nV，且S_n≤n(mu−delta)+M_Q。
事件A={S_n≤q}因而有
\[
 Q(A)\ge1-V/\omega^2=15/16.
\]

令lambda_i为可预测倾斜参数、k_i为该步条件相对熵。路径密度精确满足
\[
 \log\frac{dQ}{dP}=\sum_i k_i-N_n,\qquad
 N_n=\sum_i\lambda_i(Z_i-m_{Q,i}).\tag{7}
\]
P、Q的条件核均在正半轴有正密度，因此相互绝对连续。
N_n是Q鞅，E_Q N_n²≤4V n delta²。Chebyshev给
Q(|N_n|≤8sqrt(Vn)delta)≥15/16。
此事件与A交集的Q概率至少7/8；在交集上，(6)—(7)给似然比上界。于是
\[
 \boxed{
 P(A)\ge\frac78\exp\left\{-\frac{n\delta^2}{2(\mu-\varepsilon_\delta)}
 -8\sqrt{Vn}\delta\right\}.}\tag{8}
\]
注意需要控制N_n的绝对值或下尾，而不能用N_n的上尾替代。

由于n delta² epsilon_delta=O(ell^(5/4)n^-1/4)=o(1)一致成立，展开delta并吸收固定常数，(8)推出某个绝对A0、N0下的统一界
\[
 \boxed{P_\pi(S_n\le n\mu-d-2\log n)
 \ge\exp\{-d^2/(2\mu n)-A_0(d/\sqrt n+1)\}.}\tag{9}
\]
它对所有可预测策略及给定整个d范围成立。此处没有使用独立和的Berry–Esseen或无条件方差近似。

## 6. 小模数强迫过小低尾

取s=nmu−d，假设某个可行(h,pi)有J≤e^-s。
单调性给Phi_n(h(t))≤e^(t−s)。因此h≤n^-2于t≤s−2ell；h≤1于t≤s；于s≤t≤s+n/4有
\[
 (h(t)-1)_+\le4(t-s)/n.
\]
最后一式由y=log h及(n−1)y/(1+y)≤logPhi_n(h)≤t−s直接反解。

P下写S_n=nmu−D_path+M_P，其中D_path=Σ(mu−m(r_i))≥0、M_P为条件中心鞅。
E|M_P|≤sqrt(nV)，故E(S_n−s)_+≤d+sqrt(nV)，不需要Var S_n=O(n)。
于是中间期望超额≤4(d+sqrt(nV))/n。

远尾S_n>s+n/4的一部分h<exp(4S_n/n)。取固定theta=1/(64V)，n足够大时eta=theta+4/n≤1/2。由(5)逐步条件化，
\[
 E[e^{4S_n/n}1_{S_n>s+n/4}]
 \le\exp\{-\theta(s+n/4)+n\mu(\theta+4/n)+2Vn(\theta+4/n)^2\}\le e^{-c_1n},
\]
因为s+n/4≥nmu+n/8，而theta/8−2Vtheta²>0。这里d=O(sqrt(nell))。
另一部分h≥exp(4S_n/n)：每步条件密度f_r(z)≤2e^-z，即使动作依赖历史，联合密度仍≤2^n exp(−Σz_i)。对单纯形积分得
\[
 f_{S_n}(t)\le2^ne^{-t}t^{n-1}/(n-1)!.
\]
用Phi_n(h)≥h(4t/n)^(n−1)，高振幅期望≤
\[
 J\frac{2^n(n/4)^{n-1}}{(n-1)!}\le e^{-c_2n},
\]
因为括号n次根趋e/2<B。
逐点以基准1分解，
\[
 E h(S_n)\le1-(1-n^{-2})P(A)+4E(S_n-s)_+/n+E[h(S_n)1_{S_n>s+n/4}].
\]
综合以上，而E_pi h(S_n)>1，强迫
\[
 P_\pi(S_n\le s-2\ell)\le
 \frac{4(d+\sqrt{nV})/n+e^{-cn}}{1-n^{-2}}
 \le C_0\sqrt{\ell/n}.\tag{10}
\]
所有常数一致于h和可预测策略。

## 7. 比较两个低尾界完成一般上界

现在取d=sigma sqrt(nell)−K sqrt(n)，其中K是固定充分大的常数；n足够大时d位于(9)的范围。
(9)右边的对数至少
\[
 -\tfrac12\ell+(K/\sigma-A_0\sigma)\sqrt\ell-O_K(1).
\]
选K使括号系数正，则此概率除以sqrt(ell/n)趋于无穷，与(10)矛盾。
所以充分大n时不存在J≤exp(−nmu+sigma sqrt(nell)−K sqrt(n))的可行轮廓，
\(\log\mathcal R_n\le n\mu-\sigma\sqrt{n\ell}+K\sqrt n\)。
由第2节一般输入桥C_n≤mathcalR_n得到(1)。第3轮已证明的显式实际输入下界给匹配首项及O(sqrt n)两侧范围，证明(2)。
第8轮还提供更强下侧项sigma[loglogn+log(sigma sqrt(2pi))]sqrt(n/logn)，但本轮未在一般上界匹配它。

## 8. 本机计算与证明角色

MATLAB R2026a对delta=0.05到0.001、每个delta的501个动作网格及局部标量优化，实际求解倾斜根并比较相对熵闭式与正积分∫lambda v(lambda)dlambda：

|delta|数值最大成本/delta²|lambda/delta（该点）|
|--:|--:|--:|
|0.05|0.4012239652|0.8128745971|
|0.02|0.3950761536|0.7941800628|
|0.01|0.3930737833|0.7881384998|
|0.005|0.3920811495|0.7851521094|
|0.002|0.3914882756|0.7833711273|
|0.001|0.3912910996|0.7827792630|

其极限校准为1/(2mu)=0.39109414714；采样不是(6)的一致证明，后者由精确积分和(3)—(4)证明。
Wolfram15符号复核m'=-v、k'=lambda v和v'显式公式，三个残差均为0。
Code Analyzer无报告。复现：MATLAB `result=sm_tilt_checks;`，Wolfram `code/tilt_checks.wl`。
第11轮的三档Bellman网格作为独立机制诊断，不承担本轮定理证明。

## 9. 当前范围与下一轮

一般首个次指数缺口系数sigma现在闭合；C_n与B^n双边可比性的猜想已被(2)排除。
仍未闭合：一般C_n在sqrt(n/logn)尺度的完整第三项、O(sqrt n)上界余项的锐化、真正渐近等价及极值结构。
附件所称整个双曲示性族H_n的全空间三项公式仍需独立修补；本轮证明不以它为前提。
没有声称公开首创性。

下一轮两项动作：① 锐化反馈低尾的似然比估计或构造带Mills前因子的Bellman屏障，争取一般第三项；② 继续用非锚定/不等振幅可积输入检验上界余项及极值形态，保持真实模数与接收测度账本。
项目无冗余已跟踪文件需要清理。
