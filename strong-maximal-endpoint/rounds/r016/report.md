# 第16轮：真实下分位数的对数级改善

日期：2026-11-06（Asia/Shanghai）。状态：有界偏移等比例族的严格定理；不宣称一般端点常数的对数级项。

## 1. 结论

令Z_c=Y+cJ_c，Y~Exp(1)、J_c~Bern(p_c)、p_c=(1+e^c)^-1，二者独立。n项独立和的CDF、密度为F_{n,c}、f_{n,c}。
记L=log n、rho=sqrt(L/n)、beta_n=sigma rho，q_c为F_{n,c}(q_c)=beta_n的唯一下分位数（充分大n时0<beta_n<1）。
mu=log B、sigma=sqrt(mu)、p_*=(mu−1)/mu。对每个固定T，以下收敛一致于|t|≤T：
\[
\boxed{\frac{q_{\mu+t\rho}-q_\mu}{L}
 \longrightarrow-\frac{p_*t^2}{2}-\frac{p_*\sigma t}{2}.}\tag{1}
\]
特别地，明确指定c_n=mu−sigma rho/2便有
\[
\boxed{q_{c_n}-q_\mu=\frac{\mu-1}{8}\log n+o(\log n).}\tag{2}
\]
在任何包含−sigma/2于内部的固定t区间中，最优候选的归一化偏移趋于−sigma/2。
这证明第15轮的一项真实尾目标预测，而不是用Laplace收益推断尾收益。
不声称全局所有c的最优者必在此窗口；也未把分位数改善直接登记为C_n或Q_n的同量级改善，后者还需严格平均与实际模数结算。

## 2. 精确CDF与分位数导数

第15轮的负Laplace核M_c(lambda)满足精确恒等式
\[
 \partial_c M_c(\lambda)=p_c[M_c(\lambda)-e^{-\lambda c}].
\]
F_{n,c}的Laplace变换为M_c(lambda)^n/lambda，故
\[
 \partial_cF_{n,c}(q)=np_c[F_{n,c}(q)-F_{n-1,c}(q-c)],\quad
 q_c'=np_c\frac{F_{n-1,c}(q_c-c)-F_{n,c}(q_c)}{f_{n,c}(q_c)}.\tag{3}
\]
反演合法性可直接从有限Gamma–二项混合验证：对c的紧集，逐项微分的CDF导数有一致界，其Laplace积分绝对收敛；Laplace唯一性再由连续性给逐点等式。
取n足够大，Gamma密度在0处的相应阶导数连续，所有使用的微分均合法。q_c>0处密度严格正，隐函数微分合法。

## 3. 自包含的统一Fourier引理

将单步分布按e^-lambda Z倾斜，参数c在mu的固定小闭邻域I内、0≤lambda≤lambda0，lambda0为固定小正数。倾斜后Y速率为1+lambda，J参数为(1+exp((1+lambda)c))^-1，仍独立。
设中心化一步特征函数为phi，方差vtilde。所有固定阶绝对矩一致有界，且vtilde有正下界。连续指数因子给全局界
\[
 |\phi(u)|\le(1+u^2/(1+\lambda)^2)^{-1/2}.\tag{4}
\]
N项倾斜居中和的密度g_N因而可Fourier反演。对j=0,1,2，
\[
 \|g_N^{(j)}\|_\infty
 \le(2\pi)^{-1}\int|u|^j|\phi(u)|^Ndu
 =O(N^{-(j+1)/2}).\tag{5}
\]
最后一步可用u=(1+lambda)x代换后的Beta积分；常数一致于I及lambda范围。

还需要中心处更好的奇导数界。存在统一epsilon,b,C>0，使|u|≤epsilon时
log phi(u)=−vtilde u²/2+O(|u|³)、Re log phi(u)≤−b u²。因此
\[
 |\phi(u)^N-e^{-Nvtilde u^2/2}|
 \le C N|u|^3e^{-b'Nu^2}.\tag{6}
\]
这是用指数差的积分恒等式及两端实部的共同负二次界得到，取更小epsilon即可。
于|u|>epsilon，(4)给Fourier积分及其j≤2加权积分均为O(e^-b''N)：抽出N/2次幂在epsilon处的指数界，剩余幂的可积积分一致有界。Gaussian尾同样如此。
在g_N(0)反演中积分(6)，误差O(N∫|u|³e^-b'Nu²du)=O(N^-1)。
在g_N'(0)中先减去Gaussian反演的奇项（其积分为0），误差为O(N∫|u|⁴e^-b'Nu²du)=O(N^-3/2)。故
\[
 g_N(0)=(2\pi Nvtilde)^{-1/2}+O(N^{-1}),
 \qquad g_N'(0)=O(N^{-3/2}).\tag{7}
\]
这里不需要更精细的Edgeworth公式；尤其没有把Berry–Esseen分布误差当作密度导数误差。

## 4. 适用于临界分位数的倾斜参数

对c=mu+t rho、|t|≤T，一步均值m(c)=mu+O(rho²)、方差v(c)=mu+O(rho)，Bernstein常数一致。
第13轮核验的FGL定理在此独立同分布、确定方差情形给相对正态尾，一致于标准化距离0≤x≤2sqrt L。
比较x=(1±epsilon)sqrt L的尾与beta_n，Mills界给
\[
 x_c:=\frac{nm(c)-q_c}{\sqrt{nv(c)}}=\sqrt L(1+o(1))\tag{8}
\]
一致于窗口。因为两端的正态尾分别比beta_n大、小时都有幂次量级的严格分离。

置N=n−1、a=q_c−c。充分大n时0<a<Nm(c)，负倾斜后均值从m(c)严格降到0，所以有唯一lambda>0使N步倾斜和均值为a。
由m_c(lambda)=m(c)−v(c)lambda+O(lambda²)，结合(8)，得到
\[
 \lambda=\rho/\sigma\,(1+o(1)).\tag{9}
\]
所有误差一致于|t|≤T，且lambda最终位于第3节固定范围内。

## 5. 从密度导数到CDF参数导数

恢复原N步密度，写f=f_{N,c}。倾斜恒等式为
\[
 f(s)=M_c(\lambda)^N e^{\lambda s}g_N(s-a).
\]
(5)、(7)严格给
\[
 \frac{f'(a)}{f(a)}=\lambda+O(N^{-1}),\qquad
 \frac{|f''(s)|}{f(a)}\le C e^{\lambda(s-a)}(\lambda^2+N^{-1})\quad(s\in\mathbb R).\tag{10}
\]
第二式来自lambda²g+2lambda g'+g''及lambda/sqrt N≤(lambda²+1/N)/2；g_N(0)有cN^-1/2下界。原密度在负半轴按0延拓，N足够大时依然具有所需连续导数。

由卷积F_{n,c}(q_c)=E F_{N,c}(a−(Z_c−c))，Taylor公式和(10)给
\[
 \frac{F_{n,c}(q_c)-F_{N,c}(a)}{f(a)}
 =c-m(c)+\frac{E(Z_c-c)^2}{2}\{\lambda+O(N^{-1})\}
 +O(\lambda^2+N^{-1}).\tag{11}
\]
余项的绝对值由C E[|Z_c−c|³e^{lambda|Z_c−c|}](lambda²+1/N)控制；该期望一致有界，因为c在紧集，Z_c有一致小正指数矩。
对密度本身作同一Taylor，利用E(Z_c−c)=m(c)−c=O(rho)，得到
\[
 f_{n,c}(q_c)/f(a)=1+O(rho\lambda+\lambda^2+N^{-1})=1+O(L/n).\tag{12}
\]
这一步保留了中心化一阶矩；若只用E|Z_c−c|则会丢掉所需精度。

在窗口内c−m(c)=(c−mu)+O(rho²)、E(Z_c−c)²=mu+O(rho)。将(9)—(12)代入(3)，得到一致导数渐近
\[
 \frac{q_c'}{n\rho}=-p_*\left(t+\frac\sigma2\right)+o(1),
 \qquad c=\mu+t\rho.\tag{13}
\]
其中(11)的O(lambda²+1/n)乘n仅为O(L)=o(nrho)。p_c−p_*=O(rho)也只产生更小误差。
从0到t积分：n rho²=L，故(13)证明(1)。二次多项式严格最大点为t=−sigma/2，值p_*sigma²/8=(mu−1)/8，完成(2)。

## 6. 本机计算和适用范围

MATLAB R2026a以众数两侧12个标准差的二项权重递推、精确Gamma CDF及55步二分，比较指定偏移c_n与mu。此处“精确CDF”指分布公式，不表示浮点认证。

|n|指定偏移的分位数收益/log n|
|--:|--:|
|256|−0.00666507846|
|1024|−0.0000154272081|
|4096|0.00461694217|
|1e4|0.00693404908|
|1e5|0.01141505985|
|1e6|0.01455258045|
|1e7|0.01688810327|

理论极限(mu−1)/8=0.03480806785，趋近缓慢；小n收益可为负，不能删除这些反例性诊断。
省略二项质量以binocdf估计，从约1.69e-28到3.54e-33；并非定向截断证书。n≤4096还用c步长1e-4的中心差分检验(3)，残差约2e-9。算法确定性，无随机种子；双精度的大数相减误差未区间认证。
Wolfram15以80位对n=8,32,64的Gamma–二项混合直接微分，与(3)独立对照，均吻合至50位；单步Laplace导数恒等式符号残差为0。MATLAB Code Analyzer无报告。
复现：MATLAB `addpath('code'); sm_quantile_derivative_checks;`；Wolfram `wolframscript -file code/quantile_derivative_checks.wl`。最小输出见本轮JSON。

## 7. 差距与下一轮两项动作

本轮在真正CDF目标上证明了明确的正对数级改善，闭合第15轮相应预测，但尚未证明全局固定比例优化或可预测反馈的对数级项。
一般C_n仍只有第14轮的隐式根误差界；本轮没有改变它。
下一轮：① 将同一密度导数方法用于截断指数平均A_c(q)=E[e^((S_n−q)/n)1_{S_n≥q}]，同时结算输入模数，争取真正软轮廓下界的对数级改善；② 复核指数加权动作损失能否控制反馈相对固定比例的剩余收益。
来源：FGL使用范围见sources.md；第3节所需Fourier估计已在本文直接证明。无新增优先权主张，无冗余跟踪文件需清理。
