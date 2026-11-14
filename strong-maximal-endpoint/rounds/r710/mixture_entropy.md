# 实际接收归属产生的矩形混合熵

设正有限体积的轴平行矩形 $R_i\subset\mathbb R^n$ 构成有限目录，$\Omega=\bigcup_iR_i$，$Q=\mathrm{Unif}(\Omega)$。把 $\Omega$ 可测分割为 $U_i\subseteq R_i$，令 $P=\sum_i Q(U_i)\mathrm{Unif}(R_i)$。存在绝对常数 $C$，对每一种这样的归属都有
\[
 D(P\Vert Q)=\log|\Omega|-h(P)\le \mu n-\tfrac12\log n+C,\qquad \mu=1+e^{-\mu}. \tag{1}
\]
$h$ 是相对 Lebesgue 测度的微分熵。下文还给出有限目录和接收归属使左侧达到右侧至有界误差。本命题只对由接收归属产生的混合给出上界；实际中位数 LP 中优化的行权没有这个约束。

## 逐秩熵比较

设 $p=dP/dx$，$p_\delta=(1-\delta)p+\delta 1_\Omega/|\Omega|$，$0<\delta<1$。选 $M$ 使 $g=(\log p_\delta+M)1_\Omega$ 在 $\Omega$ 内严格为正，外延零。它是有界紧支撑非负函数。定义递减重排 $k^*(v)=\inf\{y:|\{k>y\}|\le v\}$。每个负责点上的 $M_Sg$ 至少为其负责盒上的平均；从而
\[
 \int\log p_\delta\,dP+M\le |\Omega|^{-1}\int_0^{|\Omega|}(M_Sg)^*(v)\,dv.
\]
[第11轮](../r011/report.md)的非负输入全局重排比较在秩 $v$ 以上界 $\sup_\pi\mathbb E g^*(v e^{-S_n})$ 控制被积函数，$S_n=\sum Z_{c_j}$ 允许可预测选动作。$Z_c=W+cJ_c$，$W\sim\mathrm{Exp}(1)$ 与 $J_c\sim\mathrm{Bern}(1/(1+e^c))$ 条件独立；$c=\infty$ 的离散项为零。因 $S_n\ge0$，对 $0<v<|\Omega|$ 所有被测试的秩仍在 $\Omega$ 的体积内，故精确地 $g^*(ve^{-S_n})=M+\log p_\delta^*(ve^{-S_n})$。常数 $M$ 在两侧相消。

对固定 $v,\pi$，在 $y\ge0$ 上的函数 $k_v(y)=ve^{-y}p_\delta^*(ve^{-y})$ 的积分为 $\int_0^v p_\delta^*(x)dx\le1$。与归一化后 $k_v$ 作 KL 比较给
\[
 \mathbb E\log p_\delta^*(ve^{-S_n})\le-\log v+\mathbb ES_n-h(S_n). \tag{2}
\]
$S_n$ 有有界密度和有限一阶矩（最后一步的独立指数卷积见下文），故熵有限；对固定 $v,\delta$，相应各项可积。允许每秩不同策略，不要求选取同时达到上确界的策略。下证 $\mathbb ES_n-h(S_n)\le\mu n-\tfrac12\log n+C$ 一致成立；由于 $|\Omega|^{-1}\int_0^{|\Omega|}-\log v\,dv=1-\log|\Omega|$，积分得 (1)。最后 $\delta\downarrow0$ 时 $\int \log p_\delta\,dP\to\int\log p\,dP$，有限目录的联合端点原子只有有限种密度，零密度单元不贡献。

## 受控和的均值与熵

写 $m(c)=1+c/(1+e^c)$、$\Delta(c)=\mu-m(c)\ge0$、$D=\sum_j\mathbb E\Delta(c_j)$。以每步同一个独立均匀数 $q_j$，用分位函数耦合真实 $Z_{c_j}$ 与固定动作增量 $Z^*_j=Z_\mu(q_j)$。对只依赖历史增量的策略这是逐步模拟；额外随机化可按相同条件核及额外独立随机量逐步生成。所有 $Z^*_j$ 独立同分布。设 $q_0=(1-e^{-c})/(1+e^{-c})$。积分指数混合密度可得分位函数的两支
\[
 Z_c(q)=\begin{cases}-\log[1-(1+e^{-c})q],&0<q\le q_0,\\
 -\log[(1+e^{-c})(1-q)/2],&q_0<q<1.\end{cases}
\]
在公共接缝取同值。$c$ 处于 $\mu$ 的固定紧邻域时，第一支的对数真数至少 $e^{-c}$，两支关于 $c$ 的导数一致有界。固定 $q$ 的接缝最多一次，故 $|Z_c(q)-Z_\mu(q)|\le C|c-\mu|$。$m$ 在 $\mu$ 有唯一非退化极大，局部 $|c-\mu|^2\le C\Delta(c)$。邻域外 $\Delta$ 有固定正下界（含无穷动作），每个增量均有一致二阶矩。于是该耦合在所有动作上给 $\mathbb E_q|Z_c(q)-Z_\mu(q)|^2\le C\Delta(c)$。

记 $S^*=\sum Z^*_j$；两增量差的条件均值为 $-\Delta(c_j)$。扣除条件均值后，鞅差正交及上述二阶矩给
\[
 \mathbb E|S_n-S^*|\le D+C\sqrt D .
\]
令 $K=\lfloor S_n\rfloor$、$K^*=\lfloor S^*\rfloor$。后者由 $\Gamma(n,1)$ 与独立离散量之和整分箱；Stirling 给 Gamma 密度峰高 $\le C/\sqrt n$，所以 $H(K^*)\ge\tfrac12\log n-C$。又 $H(K^*)\le H(K)+H(K^*-K)$。任意整数变量 $J$ 满足 $H(J)\le C+\log(1+\mathbb E|J|)$，可先编码符号，再与同均值的非负几何分布比较。故
\[
 H(K)\ge\tfrac12\log n-\log(1+D+C\sqrt D)-C.
\]
最后一步可按其真实条件律写 $S_n=V+W$，其中 $V=S_{n-1}+c_nJ_{c_n}$ 与最终独立指数变量 $W$ 独立。这是关于边际律的恒等式，无需在之前的分位耦合上沿用同一个 $W$。密度 $f(x)=\int_{v\le x}e^{-(x-v)}dP_V(v)\le1$ 满足 $f(x+w)\ge e^{-w}f(x)$（$w\ge0$）。令 $p_j=P(K=j)$，在 $j\le x<j+1$ 上，积分 $w\in[0,1]$ 给 $f(x)\le (p_j+p_{j+1})/(1-e^{-1})$。因此
\[
 h(S_n)\ge H(K)+\log(1-e^{-1})-\sum_{p_j>0}p_j\log(1+p_{j+1}/p_j)\ge H(K)-C
\]
（对数的凹性使最后的和不超过 $\log2$）。代入 $\mathbb ES_n=\mu n-D$，由 $-D+\log(1+D+C\sqrt D)$ 有界即得 (2) 所需估计。

## 有限矩形目录上的锐性

沿用[第212轮](../r212/report.md)的几何。固定 $n\ge2$、$T=n^4$，在 $[-1,1]^n$ 取 $A=\{L(x)>T\}$，$L(x)=\sum_i-\log|x_i|$。均匀接收点 $X\sim Q_A$ 处，按 $x_i$ 的符号反射区间 $[-r|x_i|,|x_i|]$，$r=e^{-\mu}$，乘成实际矩形，整盒都在 $A$ 中。令 $Y$ 为其均匀盒样本（不以任何入源事件条件化）。给定 $X$ 后，$L(Y)-L(X)$ 的律恒为上述 $S^*$，均值 $\mu n$，方差 $O(n)$。$W_0=L(X)-T\ge0$ 的密度为 $q(w)=c_T e^{-w}(1+w/T)^{n-1}$，其中 $1-(n-1)/T\le c_T\le1$，故 $W_0$ 有一致有界二阶矩。对 $H=L(Y)-T$ 有 $H\overset d=W_0+S^*$，右侧两项独立，因此 $\mathbb EH=\mu n+O(1)$、$\operatorname{Var}H=O(n)$。映射 $L-T$ 的数据处理及固定方差正态最大熵给
\[
 D(P_Y\Vert Q_A)\ge D(P_H\Vert q)=\mathbb EH-h(H)-\log c_T-(n-1)\mathbb E\log(1+H/T)
 \ge\mu n-\tfrac12\log n-C.
\]
最后一项的绝对费用最多 $(n-1)\mathbb EH/T$。

每个象限取二进网格，保留向上舍入绝对值所得顶点在 $A$ 的单元；收集这些顶点的动作盒，令 $\Omega_j$ 为实际盒的完整并。单元含于其动作盒且 $\Omega_j\subset A$。$A$ 的几乎所有点都是各坐标非零的严格内点，网格足细后其顶点也在 $A$，故好单元覆盖的体积趋于 $|A|$，$Q_j=\mathrm{Unif}(\Omega_j)\to Q_A$ 全变差。在好单元用其顶点盒负责，剩余点按第一件包含它的盒分配；余量体积趋零。在几乎每个严格内点，顶点盒的各边收敛至该点的动作盒，边长为正，故盒均匀律依全变差收敛。以有界连续测试积分和支配收敛得相应有限混合 $P_j\Rightarrow P_Y$。KL 对一对弱收敛概率律联合下半连续（由其有界连续测试变分表示），所以 $\liminf_jD(P_j\Vert Q_j)\ge D(P_Y\Vert Q_A)$。逐维取足够细网格，即得有限目录的匹配下例。
