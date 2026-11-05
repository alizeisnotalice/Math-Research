# 第15轮：负指数目标的精确控制与固定比例偏移

日期：2026-11-06（Asia/Shanghai）。状态：指数目标的严格机制定理；真实分位数收益系数仍待证明。

## 1. 本轮结论与适用范围

对任意lambda>0，定义
\[
 M_c(\lambda)=E e^{-\lambda Z_{e^{-c}}}
 =\frac{1+e^{-(1+\lambda)c}}{(1+\lambda)(1+e^{-c})},\quad c\in[0,\infty].
\]
存在唯一c_lambda∈(0,∞)使
\[
 e^{(1+\lambda)c_\lambda}-(1+\lambda)e^{c_\lambda}-\lambda=0.
\]
它是M_c的唯一全局最小点，而且
\[
\boxed{\inf_\pi E_\pi e^{-\lambda S_n}=e^{-n\lambda c_\lambda}.}\tag{1}
\]
可预测动作允许完整依赖过去，但这里的下确界由每一步固定使用c_lambda达到，且达到等号的策略必须每一步如此（几乎处处）。

这个严格机制障碍只针对负指数终值。它不证明真实尾概率、非减轮廓或C_n的最优策略没有反馈收益。
第11轮对真实两步尾目标已证明反馈严格优于所有固定比例；负指数是递减函数，也不是端点定义允许的非减轮廓。
本轮未改善第14轮一般C_n的误差界，保留其已证结果。

## 2. 唯一全局极小点与最优值

令t=e^c≥1、a=1+lambda。直接求导得
\[
 \partial_c\log M_c(\lambda)
 =\frac1{1+e^c}-\frac a{1+e^{ac}}
 =\frac{t^a-at-(a-1)}{(1+t)(1+t^a)}.
\]
分子F(t)=t^a−at−(a−1)满足F(1)=−2(a−1)<0、F'(t)=a(t^(a−1)−1)>0于t>1，且F(t)→∞。
故仅有一个根t_lambda>1；导数在其前负、其后正。端点c=0与c=∞的M均为1/a，排除了端点全局极小值。
在根处t^a+1=a(t+1)，所以
\[
 M_{c_\lambda}(\lambda)
 =\frac{t^a+1}{a\,t^{a-1}(t+1)}=t^{-(a-1)}=e^{-\lambda c_\lambda}=:k_\lambda.\tag{2}
\]

另有与已核验锐Lp根的精确代数联系：令p=1+1/lambda、A=exp(lambda c_lambda)，则
\[
 (p-1)A^p-pA^{p-1}-1=0.
\]
这正是sources.md中Grafakos--Montgomery-Smith一维非中心Lp范数的根方程，所以k_lambda=A^-1。
这里仅作代数识别；(1)的证明不依赖引用任何Lp范数定理。

## 3. 可预测策略的精确损失账本

令S_i=Σ_{j≤i}Z_j，S_0=0。可预测动作c_i条件于过去时，
\[
 E(e^{-\lambda S_i}\mid\mathcal F_{i-1})
 =e^{-\lambda S_{i-1}}M_{c_i}(\lambda)
 \ge k_\lambda e^{-\lambda S_{i-1}}.
\]
迭代即得(1)，固定c_lambda实现等号。更精确地，递推展开给
\[
\boxed{E e^{-\lambda S_n}-k_\lambda^n
 =\sum_{i=1}^n k_\lambda^{n-i}
 E\left[e^{-\lambda S_{i-1}}(M_{c_i}(\lambda)-k_\lambda)\right].}\tag{3}
\]
所有项非负，且指数权重严格正；结合唯一极小点可知等号当且仅当每一步c_i=c_lambda几乎处处。
这些是指数加权的动作损失，不能去掉权重后作为普通概率下的累计方差或累计均值亏损。

## 4. 小倾斜的二阶动作偏移与三阶收益

令mu=log B，满足e^mu(mu−1)=1。除以lambda后的驻点方程可以写为
\[
 F(c,\lambda)=e^c\left(\frac{e^{\lambda c}-1}{\lambda}-1\right)-1=0.
\]
它在lambda=0按可去奇点解析延拓，F(c,0)=e^c(c−1)−1；F(mu,0)=0、F_c(mu,0)=mu e^mu>0。
解析隐函数定理与泰勒系数匹配给
\[
 \boxed{c_\lambda=\mu-\frac\mu2\lambda
 +\frac{\mu(9-\mu)}{24}\lambda^2+O(\lambda^3).}\tag{4}
\]
可直接代入F展开核验每个系数，余项由F在固定邻域的解析性保证，软件不是这项保证的前提。

写K(c,lambda)=log M_c(lambda)。在c=mu处，均值和方差的精确导数为
\[
 m''(\mu)=-\frac{\mu-1}{\mu},\qquad v'(\mu)=\mu-1.
\]
从K=−lambda m(c)+lambda²v(c)/2+O(lambda³)作联合Taylor，或者将(4)代入(2)，得到
\[
 \boxed{K(c_\lambda,\lambda)-K(\mu,\lambda)
 =-\frac{\mu(\mu-1)}8\lambda^3+O(\lambda^4).}\tag{5}
\]
计算时不能只保留K的首阶−lambda m(c)；方差的一阶动作变化恰与均值的二阶损失共同确定偏移。
在mu附近固定小邻域内K_cc=−lambda m''(c)+O(lambda²)，且m''严格为负。因此存在固定a0>0使小lambda下
K(c,lambda)−K(c_lambda,lambda)≥a0 lambda(c−c_lambda)²，只要c、c_lambda都在该邻域。这是(3)可用的局部强制性，而非一般尾目标的稳定性定理。

## 5. 对数尺度的线索与尚缺步骤

取L=log n、lambda_n~sqrt(L)/(sigma sqrt n)，则(4)预测c_lambda−mu~−(sigma/2)sqrt(L/n)，(5)给n步对数Laplace收益量级L^(3/2)/sqrt n。
如果能够把这种精度的**两输入比较**转到实际下尾，再除以尾概率对阈值的对数导数lambda_n，就会出现O(L)分位数收益；形式系数为(mu−1)/8。
这一步本轮没有证明。FGL在iid情形给的相对误差正是O(L^(3/2)/sqrt n)，与所比较的Laplace收益同阶，不能用来确定这个系数；只给相对O(1/L)余项的普通鞍点公式则更不足够。
需要参数可微的余项界、两输入共同误差的抵消，或另一个精确比较恒等式。也不能把分位数改善直接作为端点下界改善，仍须结算实际Orlicz模数及严格接收平均。

## 6. 本机软件验证

MATLAB R2026a用expm1形式的驻点方程求根，并用K_c的积分计算(5)的收益，避免两个近似相等对数相减。lambda从0.5到0.001，积分AbsTol=1e-18、RelTol=1e-10，根TolX=1e-13：

|lambda|(c_lambda−mu)/lambda|收益/lambda³|
|--:|--:|--:|
|0.1|−0.600860421|−0.0400374777|
|0.01|−0.635148567|−0.0440263442|
|0.001|−0.638821247|−0.0444531270|

理论极限分别为−mu/2、−mu(mu−1)/8。Wolfram15符号核对(4)、(5)及m''、v'；在lambda=0.01、0.001独立用40/80位求根，稳定至30位，并核验最优值恒等式。

另外，MATLAB对beta_n=sigma sqrt(log n/n)的真实Gamma–二项CDF做55步分位数二分，在c∈[0.3,2.3]用fminbnd搜索，TolX=1e-8。以下是固定c候选相对c=mu的分位数收益：

|n|收益|收益/log n|
|--:|--:|--:|
|256|0.03303831771|0.00595802714|
|1024|0.06123926961|0.00883495906|
|4096|0.09344619239|0.0112345299|

所有求解状态exitflag=1，最大分位数残差1.14e-15；驻点方程残差至多5.72e-13，积分与直接对数差的偏差至多6.33e-17。具体残差记录于checks.json。搜索不是全局最优性证书，也不认证上述形式极限；这里是同维固定比例候选比较，不是与反馈最优常数比较。
确定性算法，无随机种子；双精度结果没有定向区间保证。Code Analyzer无报告。
复现：MATLAB `addpath('code'); sm_laplace_checks;`；Wolfram `wolframscript -file code/laplace_checks.wl`。

## 7. 下一轮两项动作

① 对Gamma–二项混合的动作参数导数建立可核验的尾比较，争取把第5节的O(log n)收益线索升级为有证明的机制结论。
② 用(3)的指数加权损失构造或检验尾目标的Bellman屏障，明确单一指数控制为什么不足以捕获有限维反馈收益。
一般C_n的乘法渐近等价、O(log n)以下的项及真实极值结构仍未闭合。无已跟踪冗余材料需清理；不主张公开优先权。
