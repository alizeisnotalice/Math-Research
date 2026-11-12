# 第183轮　一般源的见证内可加阈值及乘积扰动

2026-11-12。固定中央紧区间 \(I\) 与 \(D>0\)，在第169轮相同的例外价格 \(|E|e^{B_n(p)}n^{-D}\) 下，对任意源、剩余点的每件严格成功矩形，令 \(f=1_E|_R\)、\(r=\mathbb E_R f\)、\(g_i=\mathbb E_R[f\mid X_i]-r\)、\(W=\sum_i\mathbb E_R g_i^2\)。记 \(A=\sum_i g_i/\sqrt W\)、\(u=\Phi^{-1}(1-p)\)。在上述范围内一致有
\[
 P_R(\{f=1\}\triangle\{A>u\})=o(1),\qquad A\Rightarrow N(0,1).
\]
阈值函数由真实的一维条件密度给出，可随见证变化。对任意固定 \(\lambda\)，独立地以因子 \(e^{\lambda g_i}/\mathbb E e^{\lambda g_i}\) 倾斜各轴，则入源概率为 \(\bar\Phi(u-\lambda\sqrt W)+o(1)\)。

**稀疏随机分块测出大跳跃。** 由[第169轮](../r169/report.md)，\(r-p=o(1)\)、\(W-\varphi(u)^2=o(1)\)。独立以 \(\theta=1/\log n\) 抽各轴，记指标 \(B_i\)、集合 \(S\)，且与 \(X\) 独立。按[第179轮](../r179/report.md)统一于分块比例的分位筛选，对固定分块先取坏接收外包，再平均并删除坏概率超过 \(n^{-2}\) 的点。取 \(\eta=C_{I,D}\log n/\sqrt n\)，并给各例外集多留两次对数幂，总代价如上。空组、全组可直接算坏，其概率指数小于 \(n^{-2}\)。对好分块，\(b_S=\mathbb E[f\mid X_S]\) 的分布在分位 \(L^1\) 距离内距
\(\Phi((tZ-u)/\sqrt{1-t^2})\) 至多 \(\eta\)，其中 \(Z\) 标准正态、\(t^2=(|S|+t_0)/(n+2t_0)\)。目标均值为 \(p\)，所以有一件好分块即得 \(|r-p|\le\eta\)。二项集中给 \(t^2/\theta\to1\) 依概率成立。由于 \(\eta/\sqrt\theta\to0\)，分位耦合及正态分布函数在 \(-u\) 处的一阶展开表明，沿任意 \(u\to u_*\) 的序列，\((b_S-r)/\sqrt\theta\Rightarrow N(0,\varphi(u_*)^2)\)。这里概率同时包括抽轴与盒内抽样，泰勒展开可先限制 \(|Z|\le M\) 再令 \(M\to\infty\)。

以 \(W_j\) 记各阶正交能量。剔除高阶后的均方误差为
\[
 \mathbb E\big|(b_S-r)-\sum_i B_i g_i\big|^2
 =\sum_{j\ge2}\theta^j W_j\le\theta^2/4 .
\]
故 \(Y=\theta^{-1/2}\sum_iB_i g_i\) 具有同一正态极限。其特征函数在频率 1 的模对数满足
\[
 -\log|\mathbb E e^{iY}|=\theta\sum_i\mathbb E[1-\cos(g_i/\sqrt\theta)]+o(1).
\]
确实，每个因子是 \(1+z_i\)，其中 \(z_i=\theta(\mathbb E e^{ig_i/\sqrt\theta}-1)\)，\(\max|z_i|\le2\theta\)，由中心化及二阶 Taylor 余项得 \(\sum |z_i|\le W/2\)。左端趋于 \(\varphi(u_*)^2/2\)，而 \(W/2\) 同趋此值。对任意固定 \(a>0\)，当 \(|g_i|>a\) 且 \(n\) 充分大时，\(g_i^2/2-\theta(1-\cos(g_i/\sqrt\theta))\ge g_i^2/4\)；其余取值上差非负。因此 \(\sum_i\mathbb E[g_i^2 1_{|g_i|>a}]\to0\)。这一步从随机细分导出 Lindeberg 条件，单独的最大单轴方差趋零不足以代替它。因各 \(g_i\) 独立且 \(W\) 有一致正下界，Lindeberg 中心极限定理给出所述 \(A\) 的极限。

**从一阶矩恢复事件。** 正交性给 \(\mathbb E[fA]=\sqrt W\)。写 \(J=1_{A>u}\)，则
\[
 \mathbb E[|A-u|\,|J-f|]=\mathbb E[(A-u)_+]-\sqrt W+ur=o(1).
\]
这里 \(\mathbb E A^2=1\) 保证一次矩一致可积，正态正部期望为 \(\varphi(u)-u\bar\Phi(u)\)。在 \(|A-u|>h\) 上除以 \(h\)，其余由正态连续性控制，再令 \(h\downarrow0\)，得到对称差趋零。以上子列论证对中央参数一致。

**入源概率。** 因 \(|g_i|\le1\)，刚证的 Lindeberg 条件给 \(\sum_i\mathbb E|g_i|^3=o(1)\)（按固定小门槛截断再趋零），最大单项方差也趋零。对固定实数 \(z\) 展开乘积矩母函数，得 \(\log\mathbb E e^{z\sum g_i}=z^2 W/2+o(1)\)。所以倾斜后的 \(A\) 沿每个收敛参数子列趋于均值 \(\lambda\sqrt{W_*}\)、方差 1 的正态；取 \(z=2\lambda\) 还给出密度二阶矩一致有界。由 Cauchy–Schwarz，前段对称差在新律下仍趋零，入源公式成立。固定足够大的正 \(\lambda\) 可把入源概率推到任意预设的 \(1-\delta\) 以上（充分大维数）。

一般系数差仍为 \((\mu-1)/8\)。所得是任意近极值见证内的定性结构，并无趋零速率、维数增长的倾斜参数或跨矩形共用的坐标函数。下一步要估计倾斜前后各坐标改变的概率对低维事件的影响，并研究相交见证的一维条件密度如何共同约束原源的接收体积。
