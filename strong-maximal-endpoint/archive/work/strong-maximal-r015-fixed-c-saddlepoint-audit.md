# r015：固定比例下尾分位数的微小优化与所需 CDF 精度

日期：2026-11-06。本笔记只审计等比例 \(c=-\log r\) 的下尾分位数比较；不把其结论外推到一般输入、非等比例或自适应策略。

设
\[
 Z_c=Y+cJ_c,\qquad Y\sim\operatorname{Exp}(1),\quad
 \Pr(J_c=1)=\frac1{1+e^c},
\]
且令 \(m(c)=\mathbb EZ_c\)、\(v(c)=\operatorname{Var}Z_c\)。记
\(\beta_n=\sigma\sqrt{\log n/n}\)，而 \(q_c(\beta_n)\) 是 \(\sum_{i=1}^nZ_c^{(i)}\) 的下 \(\beta_n\)-分位数。

## 1. 局部优化确实预测 \(O(\log n)\) 的增益

负 Laplace 变换为
\[
 K_c(\lambda)=\log\mathbb Ee^{-\lambda Z_c}
 =-\log(1+\lambda)+\log(1+e^{-(1+\lambda)c})-\log(1+e^{-c}). \tag{1}
\]
对 \(\lambda>0\)，其关于 \(c\) 的唯一局部最小元满足
\[
 e^{(1+\lambda)c}-(1+\lambda)e^c-\lambda=0. \tag{2}
\]
在 \(\mu=1+e^{-\mu}\) 处展开得
\[
 c_\lambda=\mu-\frac\mu2\lambda+O(\lambda^2),\qquad
 K_{c_\lambda}(\lambda)-K_\mu(\lambda)
 =-\frac{\mu(\mu-1)}8\lambda^3+O(\lambda^4). \tag{3}
\]

同一个结论可从均值和方差的普通局部展开看出。Wolfram 15 实际计算得到
\[
 m''(\mu)=-\frac{\mu-1}{\mu},\qquad v'(\mu)=\mu-1.
\]
若 \(z\sim\sqrt{\log n}\) 是相应正态下尾参数，则先忽略 Edgeworth 修正的分位数近似为
\[
 F_n(c)=nm(c)-\sqrt{nv(c)}z.
\]
其局部最大元及改善为
\[
 c_{\rm opt}-\mu=-\frac{\sigma z}{2\sqrt n}+o\left(\sqrt{\frac{\log n}{n}}\right),
 \qquad F_n(c_{\rm opt})-F_n(\mu)=\frac{\mu-1}{8}z^2+o(\log n). \tag{4}
\]
故候选增益的主项正是 \((\mu-1)\log n/8\)。

Wolfram 数值直接核验 (3)：对 \(\lambda=.01,.005,.002\)，
\((K_{c_\lambda}-K_\mu)/\lambda^3\) 分别为
\(-.0440263,-.0442628,-.0444055\)，而
\(-\mu(\mu-1)/8=-.04449\ldots\)。

## 2. 为什么 Fan--Grama--Liu 与普通鞍点余项都不够

临界倾斜参数 \(\lambda\asymp\sqrt{\log n/n}\)。由 (3)，两种比例的对数尾概率只相差
\[
 n\lambda^3\asymp\frac{(\log n)^{3/2}}{\sqrt n}. \tag{5}
\]
这正是所预测的 \(O(\log n)\) 分位数位移：在该分位数处密度为 \(\asymp(\log n)/n\)，而 CDF 的绝对差为
\[
 \beta_n\,n\lambda^3\asymp\frac{(\log n)^2}{n}. \tag{6}
\]

因此 CDF 比较必须具有
\[
 \text{相对误差}=o\left(\frac{(\log n)^{3/2}}{\sqrt n}\right), \tag{7}
\]
等价地，绝对误差 \(o((\log n)^2/n)\)。FGL 的相对误差正处于
\(O((\log n)^{3/2}/\sqrt n)\) 级，不能确定 (4) 的系数。只给出
\(O(1/\log n)\) 相对余项的普通首阶 saddlepoint 近似更远远不够，不能作为本问题的证据。

## 3. 足够强、且看起来可证明的统一展开

充分条件是：存在 \(c\) 的固定紧邻域 \(I\ni\mu\)，对 \(c\in I\)、
\(0\le z\le2\sqrt{\log n}\) 一致有一项 Edgeworth 相对尾展开
\[
 \begin{aligned}
 \Pr\left\{\frac{\sum Z_c-nm(c)}{\sqrt{nv(c)}}\le-z\right\}
 ={}&\Phi(-z)+\frac{\kappa_3(c)}{6v(c)^{3/2}\sqrt n}
 (1-z^2)\phi(z)\\
 &+\phi(z)\,O\left(\frac{(1+z)^5}{n}\right), \tag{8}
 \end{aligned}
\]
or any equivalent saddlepoint formula whose *relative* remainder is
\(O((1+z)^6/n)\). At \(z\asymp\sqrt{\log n}\), the remainder is
\(O((\log n)^3/n)=o(n\lambda^3)\), so (8) resolves the comparison that FGL cannot.

The displayed Edgeworth term itself cannot be discarded: its relative size is
\(\asymp z^3/\sqrt n\), exactly the scale in (5). However, after it is retained, its variation between
\(c=\mu\) and \(c=\mu+O(\lambda)\) is only
\(O((\log n)^{3/2}/\sqrt n)\) in the *quantile* correction and does not alter the \((\mu-1)\log n/8\) gain in (4).

The model is favorable for proving (8), but this note does not claim that proof completed. Its moment generating function is analytic for \(\Re s>-1\), uniformly for \(c\in I\), and its characteristic function has the uniformly nonlattice continuous factor \((1-it)^{-1}\). These are the standard ingredients for a triangular, compact-parameter Edgeworth/saddlepoint theorem. What remains is to write or cite a theorem with the uniform moderate-deviation remainder in (8), including its dependence on the compact parameter \(c\). A leading Bahadur--Rao formula without this remainder is insufficient.

## 4. Consequence conditional on (8)

Under (8), solve the CDF equation defining \(q_c(\beta_n)\) by an implicit-function expansion. The first Edgeworth term must be kept, but is smooth in \(c\). Optimizing the resulting expression yields
\[
 \max_c q_c(\beta_n)-q_\mu(\beta_n)
 =\frac{\mu-1}{8}\log n+o(\log n), \tag{9}
\]
and the maximizing parameter obeys
\[
 c_{\rm opt}=\mu-\frac\sigma2\sqrt{\frac{\log n}{n}}
 +o\left(\sqrt{\frac{\log n}{n}}\right). \tag{10}
\]
Until a proof of (8), (9)--(10) remain a well-supported local saddlepoint prediction, not a certified improvement.
