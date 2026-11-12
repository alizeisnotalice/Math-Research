# 第302轮　用弱矩传递任意截面的中央界

2026-11-12。设 $0\le f\le1$ 可积，主坐标 $d$ 维、附加坐标 $k$ 维，$n=d+k$。采用[第301轮](../r301/report.md)假设：任一中央紧区间 $J$ 上，几乎处处的截面均有 $|\{M_x f_z>q\}|\le C_J e^{B_d(q)}\int f_z$，常数统一。则对 $p\in I\Subset(0,1)$、充分大 $d$、任意 $k\ge1$ 一致有
\[
 \log\frac{|\{M_S f>p\}|}{\int f}\le B_n(p)+O_I(1+k/\sqrt d).
\]
常数依赖稍大区间的截面常数。由此对独立变化的中心、权重和主盒得到第300轮的自由维数范围 $k=o(\sqrt d\log d)$；只需逐截面上界。下文的弱矩取上确界，从而省去第301轮积分时的 $\sqrt s$。

**精确传播幂包络。** 对全空间函数令 $W_s(g)=\sup_{t>0}t^s|\{g>t\}|$。其重排满足 $g^*(v)\le W_s(g)^{1/s}v^{-1/s}$。由[第11轮](../r011/report.md)全空间逐轴重排比较，幂函数 $v^{-1/s}$ 的单轴放大系数恰为
\[
 L_s=(1-1/s)^{-1}\sup_{0\le r\le1}\frac{1+r^{1-1/s}}{1+r}.
\]
所以 $W_s(M_z g)\le L_s^{sk}W_s(g)$，对每个有限 $s\ge2$ 成立（截断正函数再单调穷尽即可）。这是在主、辅全部坐标的乘积测度上使用重排。写 $u=1/s$，分子中的 $r^{1-u}=r-ur\log r+O(u^2)$，因为 $r^{1-u}\log^2r$ 在 $0\le u\le1/2$ 一致有界。于是 $s\log L_s=1+\max_{r}[-r\log r/(1+r)]+O(1/s)=\mu+O(1/s)$。最后等号由极值方程 $-\log r=1+r$ 及 $\mu=1+e^{-\mu}$。

**与中央体积曲线相切。** 写 $a=M_x f_z$、$L=-B'_d(p)\asymp_I\sqrt d$。取 $s=c_I\sqrt d$、$\delta=s/L$、$b=p-\delta$、$g=(a-b)_+/\delta$。则 $M_S f>p$ 迫使 $M_z g>1$。我们证明 $W_s(g)\le C'_I e^{B_d(p)}\int f$。对 $q=b+\delta t<1$，固定 $t>0$ 后按截面积分。令 $x=(q-p)/\delta=t-1$。当 $b<q\le p+\varepsilon$，相对基线的对数指数至多
\[
 s\log(1+x)+B_d(q)-B_d(p)
 \le s[\log(1+x)-x]+C_I\sqrt d\,\delta^2x^2 .
\]
在 $-1<x\le1$，亏损 $x-\log(1+x)\ge x^2/8$；选 $c_I$ 小可吸收余项。在 $x\ge1$，亏损至少 $(1-\log2)x$；选固定 $\varepsilon$ 小使 $\delta x\le\varepsilon$ 时余项可吸收。各常数取在预先选定的稍大紧区间内。对 $p+\varepsilon<q<1$ 改用 $p+\varepsilon$ 的截面上界，比基线少 $c'_I\varepsilon\sqrt d$，而 $s\log t\le c_I\sqrt d\log(1+C_I/c_I)$。先定 $\varepsilon$ 再缩小 $c_I$，使后者小于前者；同时 $b$ 留在该紧区间。$q\ge1$ 对应严格水平集为空。这证明弱矩界。收取 $sk\log L_s$ 后即得 $B_d(p)+\mu k+O_I(1+k/\sqrt d)$，与 $B_n(p)$ 的差计入同一误差。

**对一般输入的精确缺口。** 证明实际只用积分后的截面分布 $\int|\{M_x f_z>q\}|dz$。如果这一分布的统一上界带有额外因子 $e^{r_d}$，结论增加 $r_d$。这允许部分截面超过基线，只要按实际质量求和仍满足估计。一般端点上下系数差仍为 $(\mu-1)/8$。下一步检验不同分块方向的截面积分能否联合压低 $r_d$，以及能否由真实矩形的附加坐标约束改善 $k/\sqrt d$。
