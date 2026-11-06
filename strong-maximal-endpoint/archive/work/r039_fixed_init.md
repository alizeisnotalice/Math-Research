# r039：固定短时全空间初始化审计

日期：2026-11-06。结论：所提固定 \(k\ge4\) 初始化命题可严格闭合；但它只是
time-\(k\) 的上界初始化，不能替代后续全空间 Bellman 单步比较。

令 \(a=a_N\to\mu\) 为匹配动作，\(b=m(a)-a\)，
\(\delta=\mu-m(a)\). 在固定紧邻域内 \(\delta\le Cb^2\)。记
\(f=f_{k,a}\), \(v=K_a^k\phi_N\), \(J=J_k^a\)，并令
\[
 I(r)=v(r)\exp\{[\alpha J(r)+s_0f(r)]/v(r)\}
       +\Lambda k b^2/N .\tag{1}
\]
存在有限 \(s_0,\Lambda\)（只依赖 \(k,\alpha\) 和预先固定的动作紧邻域）使充分大
\(N\) 时 \(\mathcal U_{k,N}\le I\) 全空间成立。

## 1. 左支与零附近

若 \(r\le0\)，增量非负，故精确地
\[
 \mathcal U_{k,N}(r)=1-r/N+k\mu/N,
 \quad v(r)=1-r/N+km(a)/N,
 \quad f(r)=J(r)=0.
\]
选 \(\Lambda\ge\sup\delta/b^2\) 即吸收差 \(k\delta/N\)。

取固定 \(0<r_0<\inf a/2\)，且 \(r_0<\rho_N\)（充分大 \(N\)）。此时
\(\phi_N(x)\le1-x/N\) 对所有可能的 \(x=r-S\le r_0\)，从而
\[
 \mathcal U_{k,N}(r)\le L(r):=1-r/N+k\mu/N,quad 0<r<r_0.
\]
而
\[
 L(r)-v(r)=k\delta/N+D(r),
 \quad D(r)=E\{[1-e^{-(r-S)}-(r-S)/N]1_{S<r}\}.
\]
在 \(r<a/2\)，只有所有 Bernoulli 位移为零的项贡献 \(D\)。故由一次积分
\(D(r)\le Cr^{k+1}\), \(f(r)\asymp r^{k-1}\)，一致地
\[
 0\le D(r)/f(r)\le Cr^2.\tag{2}
\]
又 \(I-v\ge\alpha J+s_0f\ge s_0f\)。选择 \(s_0\ge Cr_0^2\) 和上述
\(\Lambda\) 即完成零附近；没有发生“均值差除以小 \(f\)”的非法步骤。

本机 WolframScript 15.0 实际积分核验了 (2) 的小 \(r\) 级数。以 \(q=1-p_a\)，
在 \(k=4,5,8\) 分别得到
\[
 D/f=(1-1/N)r^2/20+O(r^3),\quad
 (1-1/N)r^2/30+O(r^3),\quad
 (1-1/N)r^2/72+O(r^3).
\]

## 2. 中间紧区间

在 \([r_0,R]\)，\(f_{k,a}\) 对 \(r,a\) 连续且严格正（零 Bernoulli 项给一致正下界）。
同时 \(r\ge0\) 时 \(\phi_N(r-S)\le1+S/N\)，所以
\(\mathcal U_{k,N}\le1+k\mu/N\le1+\mu\)。先固定足够大的 \(R\)，再增大有限
\(s_0\)，使 \(s_0f\) 覆盖这个紧区间上的有界差；(1) 的指数项至少增加
\(s_0f\)。常数对充分大 \(N\) 一致（把 \(1/N=0\) 加入紧参数集即可）。

## 3. 右尾（必须用指数矩屏障）

不能以常数全局屏障收尾：\(J/v\) 的正确增长是线性而非二次。由 r036 的精确式和
\(C_k\ge0\)，对足够大 \(r\)
\[
 J/f\ge c r^2,qquad f/v\ge c/r,qquad J/v\ge\kappa r,\tag{3}
\]
其中常数对 \(a\) 的紧邻域一致。第二个不等式来自有限 Gamma--Bernoulli 混合的
统一首项界
\[
 f(r)\ge c e^{-r}r^{k-1},\qquad
 v(r)\le C e^{-r}r^k\quad(r\ge R).\tag{4}
\]
取 \(0<\kappa_0\le\alpha\kappa\) 且 \(\theta=1-\kappa_0/2\in(0,1)\)。充分大
\(N\) 时 \(\theta\ge1/N\) 且全局 \(\phi_N(x)\le e^{-\theta x}\)。故任意控制
\[
 \mathcal U_{k,N}(r)\le M_*(\theta)^k e^{-\theta r}.\tag{5}
\]
另一方面 (3)--(4) 给指数项下的 \(I(r)\ge c e^{-(1-\kappa_0)r}r^k\)，
它以严格较慢的衰减压住 (5)，只须把 \(R\) 再增大一次。这样所有常数仍仅依赖
\(k,\alpha\)。

这也修正了“用 \(O^-\) 的常数尾上界收尾”的错误直觉：常数不能被衰减的 \(I\)
支配；必须使用 (5)。

## 4. 范围

由三区拼接，(1) 是固定短时的真值上界。又全局 \(O_k\) 支配真值，
\(\min(O_k,I)\) 仍支配真值，适合作为之后的初始化数据。它没有证明
\(T I_{t-1}\le I_t\) 于所有后续状态，故不构成 global one-step closure，也不产生
新的 sharp 一般上界结论。
