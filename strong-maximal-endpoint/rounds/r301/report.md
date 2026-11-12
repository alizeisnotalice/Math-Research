# 第301轮　只按各截面的中央水平集计价

2026-11-12。设 $0\le f(x,z)\le1$ 可积，$x\in\mathbb R^d,z\in\mathbb R^k$，$n=d+k$。若几乎每条截面在每个中央紧区间 $J$ 一致满足
\[
 |\{M_x f_z>q\}|\le C_J e^{B_d(q)}\int f_z,\qquad q\in J,
\]
则对每个 $I\Subset(0,1)$、$p\in I$、$2\le s\le d^{1/4}$，有
\[
 \log\frac{|\{M_S f>p\}|}{\int f}\le B_d(p)+\mu k+O_I(1+k/s)+\tfrac12\log s.
\]
常数亦依赖稍大紧区间的 $C_J$。假设只涉及逐截面总体积，无需不同截面共享排序、峰位、权重、支撑盒或中心。特别地 $k=d^{o(1)}$ 时取 $s=\max(2,k)$（充分大 $d$），即得 $B_n(p)+o(\log n)$。这包含每条截面有独自的任意非负深度权重、中心及轴平行支撑盒的单调深度输入：单条先由第285轮同质量替换，再由第111轮全空间加权定理与仿射不变性得到假设。同一主盒上异常轴的可能位置可以遍及全部 $d$ 轴。与第300轮相比，本轮对自由维数的无损范围较小。

**选门槛使两项的一阶导数相消。** 设 $a(x,z)=M_x f_z(x)$，用有理矩形保证可测。$M_S f\le M_z a$。写 $L=-B_d'(p)\asymp_I\sqrt d$，$\delta=(s-1)/L$，$b=p-\delta$。若 $M_z a>p$ 则 $M_z (a-b)_+>\delta$。第275轮的逐轴 $L^s$ 模数给总接收体积至多
\[
 e^{\mu k+O(k/s)}\int ((a-b)_+/\delta)^s\,dx\,dz .
\]
以下固定 $z$，写 $V(q)=|\{a(\cdot,z)>q\}|$。层积分为 $\int_b^1(s/\delta)(1+(q-p)/\delta)^{s-1} V(q)\,dq$。

**均匀积分界。** 在 $[b,p+\varepsilon]$ 代入截面上界，$x=(q-p)/\delta$，提取 $e^{B_d(p)}\int f_z$。Taylor 余项至多 $C_I\sqrt d\,\delta^2 x^2$，线性项为 $-(s-1)x$。在 $-1<x\le1$，$x-\log(1+x)\ge c x^2$（例如 $c=1/8$），故余项由负二次项的一半吸收，因 $\sqrt d\,\delta^2/(s-1)=O_I(s/\sqrt d)$。在 $x\ge1$，$x-\log(1+x)\ge(1-\log2)x$；当 $\delta x\le\varepsilon$，余项至多 $C'_I\varepsilon(s-1)x$。先取固定 $\varepsilon>0$ 足够小，使这些门槛同属固定中央紧区间且余项可吸收。因此该积分至多
\[
 C_I s\left(\int_{-1}^1 e^{-c(s-1)x^2}\,dx+\int_1^\infty e^{-c(s-1)x}\,dx\right)\le C'_I\sqrt s .
\]
对 $q>p+\varepsilon$，用 $V(p+\varepsilon)$，原层积分在此段至多 $V(p+\varepsilon)((1-b)/\delta)^s$。前者相对基线节省 $c_I\sqrt d$ 的对数体积，后者对数为 $O_I(s\log d)=o(\sqrt d)$，故远端被吸收。积分 $z$ 即得定理。阈值的选择明确给出半对数项的来源；本轮没有证明此项对真实算子必要。

全体一般源尚缺逐截面中央假设；若右端统一多出 $e^{r_d}$，本论证将 $r_d$ 原样带入结论。一般端点系数差仍为 $(\mu-1)/8$。下一步检验能否直接使用真实穿点辅助矩形以去掉 $\sqrt s$ 积分费用，并研究近最优源是否在适当分块上满足逐截面假设。
