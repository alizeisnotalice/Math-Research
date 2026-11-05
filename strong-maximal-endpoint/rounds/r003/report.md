# 第 3 轮：指数软轮廓与 log log n 级严格收益

日期：2026-11-06（Asia/Shanghai）。状态：严格下界定理，两名独立子代理审查通过；没有确定一般输入的第二项。

## 1. 严格改进

沿用 mu=log B、sigma=sqrt(mu)、log B=1+1/B。对每个整数 n≥1024，
\[
\boxed{\log C_n\ge n\mu-\sigma\sqrt{n(\log n-2\log\log n)}-\log(2n).}\tag{1}
\]
因此存在 epsilon_n→0，使
\[
\log C_n\ge n\mu-\sigma\sqrt{n\log n}
+\sigma\log\log n\sqrt{n/\log n}+\epsilon_n\sqrt{n/\log n}.\tag{2}
\]
新增的正 log log n 项超过附件单层展开的固定第三项；本轮不依赖附件全局 H_n 上界。
首个负项系数仍为 sigma。若要声称优于整个双曲示性族，还须补齐该族全局上界；不作公开首创宣称。

## 2. 自包含的倾斜尾界

固定 r=B^-1=e^-mu，Z=Y+mu J，其中 Y~Exp(1)、J~Bernoulli(r/(1+r)) 独立。
E Z=Var Z=mu。令 T=ΣZ_i、W=mu−Z、S_n=ΣW_i。
对 lambda≥0，
\[
K(\lambda)=\log E e^{\lambda W}
=\lambda\mu-\log(1+\lambda)+\log(1+e^{-\mu(1+\lambda)})-\log(1+e^{-\mu}).
\]
记 p_lambda=(1+e^(mu(1+lambda)))^-1<1/2，微分得到
\[
K^{(3)}(\lambda)=-2/(1+\lambda)^3-\mu^3p_\lambda(1-p_\lambda)(1-2p_\lambda)<0.
\]
K(0)=K'(0)=0，K''(0)=mu，故 K(lambda)≤mu lambda²/2 对所有 lambda≥0 成立。

以 dP_lambda/dP=exp(lambda S_n−nK(lambda)) 倾斜后，Y_i~Exp(rate=1+lambda)、J_i~Bernoulli(p_lambda) 仍独立。
S_n 的密度 f_lambda 是平移、反射 Gamma(n,rate=1+lambda) 密度的混合。
Gamma 众数和 Stirling 阶乘下界给出（n≥2）
\[
\|f_\lambda\|_\infty\le(1+\lambda)(n-1)^{n-1}e^{-(n-1)}/(n-1)!
\le(1+\lambda)/\sqrt{2\pi(n-1)}.
\]
对 d>0 反倾斜并积分：
\[
P(S_n\ge d)=e^{nK-\lambda d}\int_0^\infty e^{-\lambda y}f_\lambda(d+y)dy
\le e^{nK-\lambda d}\|f_\lambda\|_\infty/\lambda.
\]
取 d=sigma sqrt(n) z、lambda=z/(sigma sqrt(n))，得到所有 z>0 的界
\[
\delta=P(T<n\mu-\sigma\sqrt n z)
\le\frac{(1+\lambda)\sigma\sqrt n}{z\sqrt{2\pi(n-1)}}e^{-z^2/2}.\tag{3}
\]
T 有密度，事件端点无影响。1/z 前因子来自密度界，不需要引用中偏差定理。

## 3. 软轮廓、余量与模数

取 H_n(u)=exp(u/n)1_(u≥0)、q=nmu−sigma sqrt(n)z。Jensen 给完整指数的期望≥exp((nmu−q)/n)，删去 T<q 部分的损失≤delta，因此
\[
E H_n(T-q)\ge e^{\sigma z/\sqrt n}-\delta.\tag{4}
\]
只要 delta<sigma z/sqrt(n)，平均就严格超过1。
设 ell=log n、z²=ell−2log ell。由 e^(-z²/2)=ell/sqrt(n)，(3) 除以 sigma z/sqrt(n) 的上界为
\[
R_n=(1+\lambda)\sqrt{n/(n-1)}\frac{\ell}{\sqrt{2\pi}(\ell-2\log\ell)}.
\]
当 n≥1024，sigma>1，lambda≤sqrt(ell/n)<0.083，sqrt(n/(n−1))<1.001。
ell/(ell−2log ell) 在 ell>e 上递减；ell>6.9、log(6.9)<1.932 使此比值<2.28。
由 sqrt(2pi)>2.5，R_n<1.083×1.001×2.28/2.5<1。
log/exp 数值分界由 soft_constants.py 的精确分数 Taylor 余项验证；pi>25/8 是标准初等界。
同时 q>0：sqrt(ell/n)<0.083，且 0.083 sigma<sigma²=mu。

固定维数的 Orlicz 模数极限为
\[
I_n=\int_0^\infty e^{-(1-1/n)u}(1+u/n)^{n-1}du
=\sum_{k=0}^{n-1}\frac{(n-1)_k n^{-k}}{(1-1/n)^{k+1}}\le\frac{n^2}{n-1}\le2n.
\]
这里 (n−1)_k 是下降阶乘，逐项≤n/(n−1)。

为保持实际函数有界，可先取 H_(n,M)=min(H_n,M)。固定 n,q 后，单调收敛保证某个有限 M 仍有严格平均>1，模数积分≤I_n。
取实际输入 f_L=H_(n,M)(S(x)−L)1_[-1,1]^n，所选矩形在 E_(L-q) 内均给严格见证。
随后固定 n,M,q 再令 L→∞。第2轮支配收敛的证明适用于一切有界、非减、左侧最终为零的轮廓，得到 C_n≥exp(q)/I_n≥exp(q)/(2n)，证明 (1)。
亦可用有限阶梯从下逼近 H_(n,M) 保持余量。没有交换跨维极限，也没有用原子代替可积输入。

## 4. 展开与实际验证

令 b=log ell，有
\[
\sqrt{n(\ell-2b)}=\sqrt{n\ell}-b\sqrt{n/\ell}
+O(b^2\sqrt n/\ell^{3/2}).
\]
余项除以 sqrt(n/ell) 为 O(b²/ell)→0，log(2n)/sqrt(n/ell)→0，故 (1) 推出 (2)。

MATLAB 对 g(t)=exp(a(t−q)/n)1_(t≥q)，a∈{0.25,0.5,0.75,0.9,1}，固定 r=B^-1，计算精确正倾斜 Gamma–二项卷积。
二分选择 q，使对数期望大于1e-10。一般 a 的模数由精确有限和
\(I_n(a)=\sum_{k=0}^{n-1}(n-1)_k(a/n)^k/(1-a/n)^{k+1}\) 求得，评分为 q−log I_n(a)。
以下均为数值候选；定理 (1) 则独立严格成立。

|n|a|log候选|标准化缺口|
|--:|--:|--:|--:|
|4|1|3.343980712|0.7515978|
|8|1|7.343248058|0.7072084|
|16|0.9|15.894152361|0.6848320|
|32|0.9|33.800179369|0.6752099|
|64|0.9|70.826812332|0.6739283|
|128|0.9|146.730613231|0.6786565|
|256|0.9|301.375682351|0.6877182|
|512|0.9|615.027101566|0.6997482|
|1024|1|1249.034148846|0.7135257|
|2048|1|2527.511212204|0.7265008|
|4096|1|5100.021838826|0.7398924|

MATLAB R2026a 实际运行到4096维，无随机采样；Code Analyzer 无报告。
Wolfram 15 对 n=16,256,1024,4096 以60位复算，均保留正余量，对数评分差最大2.77e-12。
Wolfram 符号核验三阶导数、均值与方差恒等式；50位数值积分核对 n=4,16,64、a=1/2,1 的有限和。
这些数值是实现检查，不是定向证书；严格论证在第2—3节，常数检查采用精确有理数。

## 5. 当前差距与下一轮

当前严格范围是 (2) 的下界与 log C_n≤nmu+mu/2+o(1) 的一般上界。
(2) 是带该展开的下界序列，不是声称 log C_n 自身有此展开。
尚未改变首个负项系数、确定一般第二项或证明 C_n 与 B^n 可比。

下一轮两项动作：① 为有限阶梯建立定向算术证书；② 分析完整凸轮廓与指数轮廓的机制边界，并继续研究一般上界近阈值核心。
复现：MATLAB `rows=sm_exponential;`；Wolfram `code/softproof_checks.wl`、`code/exponential_check.wl`；Python `code/soft_constants.py`。
本轮无冗余已跟踪文件需要清理。
