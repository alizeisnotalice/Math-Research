# 第174轮　允许更换全部真实见证后的覆盖障碍

2026-11-12。固定中央紧区间 \(I\Subset(0,1)\)。在[第173轮](../r173/report.md)取参数 \(a=1\) 的源 \(E=F\times S\) 与目标接收区 \(U=U_0\times V\) 上，有 \(\log(|U|/|E|)\ge B_n(p)-216\mu\log n-O_I(1)\)。对任何满足 \(t_n\sqrt n\to\infty\)、\(0<t_n\le1\) 的序列，任取有限个阈值 \(p\) 的真实成功矩形，任意排序，令 \(A_j\) 为第 \(j\) 个矩形在全空间的新增部分。则一致有
\[
 |U|^{-1}\sum_{|A_j|\ge t_n|R_j|}|A_j\cap U|=o_I(1).
\]
矩形不必属于原先构造的族，也不必包含在 \(U\) 中。特别覆盖了每个 \(t_n=n^{-a}\), \(0<a<1/2\)。

**基底门槛带。** 用 \(b=M_d1_F\)，\(m=\lceil216\log n\rceil\)，\(d=n-m\)，\(\beta=e^{-m/72}\)。第173轮取等权双曲十字，给 \(|U_0|\ge |F|\exp(B_d(p/(1-\beta))-C_I)\)。[第111轮全空间上界](../r111/central-cross111.md)对同一实际截断给 \(|\{b>h\}|\le |F|\exp(B_d(h)+C_I)\)。对小的 \(\delta>0\)，中央正态分位导数的绝对值有一致正下界，而 \(\sqrt n\beta=o(1)\)，故
\[
 |\{b>p\}|\le C_I|U_0|,\qquad
 |\{b>p+\delta\}|\le C_I e^{-c_I\sqrt n\delta}|U_0|.
\]
在 \(U\) 内保留 \(W=(U_0\cap\{b\le p+\delta\})\times (V\setminus S)\)，其补集相对体积至多 \(C_I e^{-c_I\sqrt n\delta}+\beta\)。

**对替代见证计费。** 写每个真实矩形为 \(R_j=P_j\times Q_j\)。只考虑新增占比至少 \(t=t_n\) 且触及 \(W\) 的盒。盒内平均是两个实际平均之积。因 \(P_j\) 含一点满足 \(b\le p+\delta\)，其平均至多 \(p+\delta\)，故 \(|Q_j\cap S|/|Q_j|>p/(p+\delta)\)。同时 \(P_j\subset\{b>p\}\)（边界可忽略）。记 \(\epsilon=\delta/(p+\delta)\)、\(D=\{b>p\}\times V\)，则 \(|R_j\setminus D|\le\epsilon|R_j|\)，且 \(|R_j\cap W|\le\epsilon|R_j|\)。从 \(|A_j\cap D|\ge(t-\epsilon)|R_j|\) 及新增区域互不相交得
\[
 \sum_j |A_j\cap W|\le \frac{\epsilon}{t-\epsilon}|D|
 \le C_I\frac{\epsilon}{t-\epsilon}|U|.
\]
未触及 \(W\) 的盒和其余在 \(U\setminus W\) 中的新增部分合计按互不相交只收一次补集体积。令 \(\delta=\sqrt{t/\sqrt n}\)，则 \(\delta/t\to0\)、\(\sqrt n\delta\to\infty\)，得到结论。所有界独立于盒数和排序。

一般系数差仍为 \((\mu-1)/8\)。第171轮即使对固定少数坐标也使用约 \(\log n/\sqrt n\) 的新增占比尺度，该尺度已被上述全见证障碍覆盖。此族仍付出了明确的 \(O(\log n)\) 维数代价。下一步需估计达到 \(B_n(p)+o(\log n)\) 时这种门槛带与分布分离能否共存，并比较按新增区域的实际源密度计费与单按其占比筛选的差异。
