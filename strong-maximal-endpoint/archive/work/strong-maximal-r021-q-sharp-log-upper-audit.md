# r021：独立比例向量与单调轮廓常数 \(Q_n\) 的 sharp \(\log n\) 上界审计

日期：2026-11-06。这里的 \(Q_n\) 仅指独立坐标比例向量与非减标量轮廓构成的模块常数。本文不把任何结论用于真实强极大算子常数 \(C_n\)；从 \(Q_n\) 到 \(C_n\) 的推广仍需要独立的一般输入归约。

下界候选（r019 加 r017）为
\[
 \log Q_n\ge n\mu-d_n+
 \left\{\frac{3-\mu}{6}+\frac{\mu-1}{8}\right\}\log n+o(\log n)
 =n\mu-d_n+\frac{9-\mu}{24}\log n+o(\log n). \tag{1}
\]
本笔记审计匹配上界的可行路线。

## 1. 从点态 Orlicz 预算到尾部不等式

令 \(h\) 非减且
\(J=\int e^{-t}\Phi_n(h(t))dt\le e^{-s}\)。单调性给出点态预算
\[
 \Phi_n(h(t))\le e^{t-s}. \tag{2}
\]
在 \(t<s\) 时，\(h(t)\le e^{t-s}\)，所以对任意独立比例向量的和 \(S\)，
\[
 \mathbb E(1-h(S))_+
 \ge F_S(s)-\int_{-\infty}^s e^{t-s}f_S(t)dt. \tag{3}
\]
式 (3) 中的积分不能只由 \(e^tf_S(t)\) 单调性粗估；那会多出一个长度因子。正确做法是在使
\(s\) 成为均值的负指数倾斜下写
\[
 \frac{f_S(s-u)}{f_S(s)}=e^{-\lambda u}
 \frac{g_\lambda(s-u)}{g_\lambda(s)}.
\]
三角阵列的倾斜 Fourier 局部极限定理给出 \(g_\lambda(s)\asymp n^{-1/2}\)、
\(\sup g_\lambda=O(n^{-1/2})\)，故该积分为 \(O(f_S(s)/(1+\lambda))\)。在临界窗口它是
\(O(\log n/n)\)，只影响阈值的 \(O(1)\) 项。

在 \(t>s\) 的中段，反演 (2) 给出
\[
 h(t)-1\le\frac{t-s}{n}+O\left(\frac{(t-s)^2}{n^2}\right)
 \quad(0\le t-s\le n^{3/4}). \tag{4}
\]
结合远尾的正矩母函数/密度分割，若
\(M=\mathbb ES=n\mu-D\)，则
\[
 \mathbb E(h(S)-1)_+
 \le\frac{((n\mu-s)-D)_+}{n}+O\left(\frac{\log n}{n}\right). \tag{5}
\]
因此由 \(\mathbb Eh(S)\ge1\)，令 \(d=n\mu-s\)，得到基础约束
\[
 F_S(s)\le\frac{(d-D)_+}{n}+O\left(\frac{\log n}{n}\right). \tag{6}
\]
正式证明 (5) 时必须保留 \(D\)：把 \(S-s=d-D+M_0\) 中的 \(D\) 直接丢掉会破坏 sharp 系数。

## 2. 近最优比例向量的三角阵列尾展开

设 \(d_i=\mu-m(r_i)\)、\(D=\sum d_i\)、\(V=\sum v(r_i)\)。在近最优区，
\[
 V=\mu n+\Delta V,\qquad |\Delta V|\le C\sqrt{nD},\qquad
 K_3=\sum\kappa_3(r_i)=n\mu(3-\mu)+O(\sqrt{nD}). \tag{7}
\]
固定 \(s=n\mu-d\)，令 \(y=d-D\) 和 \(x=y/\sqrt V\)。需要的充分 CDF 引理是
\[
 F_S(s)=\Phi(-x)+\frac{K_3}{6V^{3/2}}(1-x^2)\phi(x)+O(n^{-1}), \tag{8}
\]
在 \(x\le2\sqrt{\log n}\) 且所有近最优向量上一致。它可由每一坐标的指数特征函数因子给出：负倾斜后连续部分提供 \((1+t^2)^{-n/2}\) 型 Fourier 外区衰减，内区用三阶累积量 Taylor 展开。此处必须是绝对 \(O(n^{-1})\) 余项；FGL 的相对尾误差不足以确定 \(\log n\) 系数。

在近最优情形 \(D<d\) 下，把 (8) 代入 (6) 并解隐式方程，得到
\[
 d\ge d_n-a\log n+D+\frac{\lambda_n}{2}\Delta V+o(\log n),
 \quad a=\frac{3-\mu}{6},\quad
 \lambda_n=\frac{d_n}{\mu n}\sim\frac1\sigma\sqrt{\frac{\log n}{n}}. \tag{9}
\]
式 (9) 的符号很关键：\(D\) 是均值损失，而负的方差变化可通过
\((\lambda_n/2)\Delta V\) 回补。

## 3. 坐标局部优化与 \((\mu-1)/8\) 系数

对每一个坐标，考察
\[
 H_\lambda(r)=\mu-m(r)+\frac\lambda2(v(r)-\mu). \tag{10}
\]
当 \(\lambda\downarrow0\)，其全局最小元位于唯一最大均值点 \(r=B^{-1}\) 的邻域；离开该邻域，
\(\mu-m(r)\) 有固定正缺口，不能被 \(O(\lambda)\) 项补偿。用
\[
 m''(\mu)=-\frac{\mu-1}{\mu},\qquad v'(\mu)=\mu-1
\]
在 \(c=-\log r\) 变量展开，得到
\[
 \min_rH_\lambda(r)=-\frac{\mu(\mu-1)}8\lambda^2+O(\lambda^3),
 \qquad c_{\rm min}=\mu-\frac\mu2\lambda+O(\lambda^2). \tag{11}
\]
把 (11) 对 \(n\) 个坐标求和并代入 \(n\lambda_n^2=(\log n)/\mu+o(\log n)\)，有
\[
 D+\frac{\lambda_n}{2}\Delta V
 \ge-\frac{\mu-1}{8}\log n+o(\log n). \tag{12}
\]
由 (9)--(12) 推出所需的条件上界
\[
 \log Q_n\le n\mu-d_n+
 \left\{a+\frac{\mu-1}{8}\right\}\log n+o(\log n)
 =n\mu-d_n+\frac{9-\mu}{24}\log n+o(\log n). \tag{13}
\]

## 4. 尚需补齐的严格环节

本路线的局部优化 (10)--(12) 是完整且全局排除远离 \(B^{-1}\) 参数的。匹配 (13) 尚依赖三个可单独写成引理的分析步骤：

1. 在 (5) 中保留随机/可预测均值缺口 \(D\) 的高侧预算，而不是以 \(\mathbb E|S-M|\) 的粗界替代。
2. 对近最优三角阵列证明 (8) 的统一绝对 \(O(n^{-1})\) Fourier--Edgeworth 余项，以及倾斜密度比界以完成 (3)。
3. 由 (6)、(8) 到 (9) 的隐式方程展开须明确处理 \(D\) 较大时的排除；(10) 的全局最小化表明其不可能改善系数，但需要把该排除先接入尾部不等式。

因此 (13) 是一条结构完整、与 r019 下界匹配的 sharp-\(\log n\) 上界蓝图，尚不是已完成定理。它仅涉及 \(Q_n\)，不应替代一般输入 \(C_n\) 的上界工作。
