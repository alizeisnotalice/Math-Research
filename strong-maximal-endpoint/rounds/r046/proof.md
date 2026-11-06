# 第46轮：实际模数预算下的逐块外带吸收

日期：2026-11-07。以下均为固定参数的独立块；不包括逐步反馈选参。
令 \(\mu=1+e^{-\mu}\)，\(L=\log N\)，
\(\bar\Phi(z_N)=\sqrt{\mu/N}z_N\)，\(d_N=\sqrt{\mu N}z_N\)。
取固定 \(m\)，整数 \(b_j\asymp N\)，\(\sum b_j=N\)。
块 \(T_j\) 是 \(b_j\) 个独立
\(Z_{c_j}=\operatorname{Exp}(1)+c_j\operatorname{Bern}((1+e^{c_j})^{-1})\)
之和，不同块独立，\(c_j\ge0\) 在抽样前确定。记
\[
 m(c)=1+\frac c{1+e^c},\qquad
 \nu(c)=1+\frac{c^2e^c}{(1+e^c)^2}.
\]
有 \(m(c)\le m(\mu)=\mu\)、\(\nu(\mu)=\mu\)。

## 定理及范围

设任意非负可测函数 \(H\) 满足实际预算
\[
 J(H)=\int_{\mathbb R_+^m}e^{-\sum t_j}\Phi_N(H(t))\,dt\le e^{-s},
 \quad s=N\mu-d_N+O(L),\quad
 \Phi_N(h)=h(1+\log^+h)^{N-1}.
\]
对每个固定 \(K>0\)，存在固定 \(A\)，使得一致地
\[
 E_{\bf c}\left[H(T)1_{\{\exists j:\ |T_j-b_jm(c_j)|>A\sqrt{b_jL}\}}\right]
 \le N^{-K} \tag{1}
\]
在充分大的 \(N\) 成立。无需 \(H\) 单调、有界或分离变量。
若额外 \(|c_j-\mu|\le K_0\sqrt{L/N}\)，可将中心换为 \(b_j\mu\)，
且**任意 \(A>\sqrt{2\mu K}\) 都可用**。
低于该严格阈值的普遍结论不成立；等号边界本轮不判定。
所有一致常数可依赖固定 \(m,K,K_0\)、块比例下界及预算中的 O 常数。

## 一次性 Fenchel 结算

写 \(U=\sum t_j\)，\(d\xi=e^{-U}dt\)，
\(R=e^U\prod f_{b_j,c_j}(t_j)\)。单增量密度满足
\(q_c(t)\le2e^{-t}1_{t>0}\)，故
\[
 R\le2^N\prod_j\frac{t_j^{b_j-1}}{\Gamma(b_j)}
 \le2^N\left(\frac{eU}{N-m}\right)^{N-m}.\tag{2}
\]
第二步用 \(k!\ge(k/e)^k\) 和权重 \(k_j=b_j-1\) 的 AM–GM；
这里保留 Gamma 分母，避免遗漏块尺度。

令 \(E\) 为(1)中的整个并事件，\(\lambda=e^sN^{-(K+2)}\)。
由(2)和 \(s=N\mu+o(N)\)，
\((R/\lambda)^{1/(N-1)}\le C_K(1+U/N)\)。
对于 \(h=e^y>1\)，
\[
 \Phi'_N(h)=(1+y)^{N-2}(N+y)\ge y^{N-1}.
\]
因此共轭函数对 \(a\ge0\) 满足
\(\Phi_N^*(a)\le a\exp(a^{1/(N-1)})\)。
当最优点为折点 \(h=1\) 或 \(h=0\) 时同一界仍成立。
Young–Fenchel 只在整个并事件用一次，得到
\[
 E[H1_E]\le\lambda J(H)+\lambda\int_E\Phi_N^*(R/\lambda)d\xi
 \le N^{-(K+2)}+E[e^{C_K(1+U/N)}1_E].\tag{3}
\]
没有按块或角区域重复花费预算。

## 均匀概率控制与锐常数

对固定 \(\theta_0<1\)，\(Z_c\) 在 \(|\theta|\le\theta_0\) 上的指数矩、
前若干阶倾斜累积量对所有 \(c\ge0\) 一致有界：Bernoulli 项的导数由
\(c^re^{-(1-\theta_0)c}\) 控制。因而
\[
 \log E e^{\theta(Z_c-m(c))}\le C\theta^2\quad(|\theta|\le\theta_0),
 \qquad \log E e^{\theta Z_c}\le\mu\theta+C\theta^2\quad(0\le\theta\le\theta_0).
 \tag{4}
\]
第二条只用于正 \(\theta\)；不将 \(m(c)\le\mu\) 用于负号。
Bernstein 给 \(P(E)\le C_mN^{-cA^2}\)。对任意固定 \(p>1\)，(4)又给
\(E\exp(pC_K U/N)=O_{p,K}(1)\)。Hölder 将(3)右端第二项压到
\(O(P(E)^{1-1/p})\)，先取 \(p\)，再充分增大 \(A\)，即得(1)。

在近 \(\mu\) 窗口，\(m'(\mu)=0\)，故各中心漂移为 \(O(L)\)，
\(\nu(c_j)=\mu+o(1)\)。在偏差 \(A\sqrt{b_jL}\) 处对(4)保留二次系数，
三次余项为 \(O(L^{3/2}/\sqrt N)=o(1)\)，得
\[
 P(\exists j:|T_j-b_j\mu|>A\sqrt{b_jL})
 \le N^{-A^2/(2\mu)+o(1)}.
\]
当 \(A^2/(2\mu)>K\) 时选有限 \(p\) 使
\((1-1/p)A^2/(2\mu)>K\)，代入(3)；严格指数余量吸收固定常数。
这证明所述锐阈值的充分方向。

## 必要方向及常数宽度障碍

取偶数 \(N\)，\(b=N/2\)，两块均取 \(c=\mu\)，\(\tau=N\mu-d_N\)。
给定 \(A<B<\sqrt{2\mu K}\)，令
\[
 r_1=b\mu+B\sqrt{bL},\quad r_2=\tau-r_1>0,\quad
 H=1_{\{t_1\ge r_1,t_2\ge r_2\}}.
\]
直接积分给 \(J(H)=e^{-(r_1+r_2)}=e^{-\tau}\)，且 H 完全位于半径 A 的外带。
而
\[
 EH=P(T_1\ge r_1)P(T_2\ge r_2)
   =N^{-B^2/(2\mu)+o(1)}(1-o(1)).\tag{5}
\]
这里第二因子用 Chebyshev 即趋于1。第一因子无需引用未经核验的中偏差定理：
令单块 log-MGF 为 \(b\Lambda(\theta)\)，按均值
\(r_1+a\sqrt b\) 倾斜，其中固定 a 足够大。倾斜后方差为
\(b(\mu+o(1))\)，Chebyshev 在 \([r_1,r_1+2a\sqrt b]\) 给固定正质量。
反向似然的 log 为 \(-B^2L/(2\mu)+O(\sqrt L)+o(L)\)；
普通 Chernoff 给匹配上界，证明(5)。其衰减严格慢于 \(N^{-K}\)。
若将 \(r_1\) 改成 \(b\mu+a\sqrt b\)，则相同精确预算下
\(EH\to\bar\Phi(a/\sqrt\mu)>0\)，排除常数宽度的普遍删除。
这些 H 都只有0、1值，EH<1，不是端点可行见证或 C_N 新下界。

## 平移与缺口

把 H 在正象限外延为0，对任意固定 \(u\in\mathbb R^m\)，换元得
\(J(H(u+\cdot))\le e^{\sum u_j}J(H)\)。只要
\(s-\sum u_j=N\mu-d_N+O(L)\)，同一证明适用。
本定理仅删除外带，没有控制中心管内连续矩形选参的熵、共享侧边关联或锐严格余量。
