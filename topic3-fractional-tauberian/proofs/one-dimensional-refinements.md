# 一维精细渐近：证明与计算记录

本笔记仅记录第\ `one:all-orders`、`one:transition` 和
`one:critical-window` 的推导与可复现计算；正文中的结论不依赖浮点数值。

## 固定参数的全阶展开

令 \(t=1-\alpha\)、\(a=1/\beta\)、\(q=1/(1-\beta)\)。精确闭式给出

\[
A_\beta(t)=1-(1-t)^a,
\qquad D_\beta(t)=(1-t)^{-q}-1.
\]

故第 \(j\) 阶系数分别是

\[
(-1)^{j+1}\binom aj,\qquad \frac{(q)_j}{j!}.
\]

对任何紧参数区间 \(K\Subset(0,1)\)，这两个函数在
\([0,t_0]\times K\) 上有有界的 \(N+1\) 阶导数。正文列出的 Lagrange
界给出 \(O_{K,N,t_0}(t^{N+1})\) 的统一余项；因此这里的全阶展开并非逐个固定
\(\beta\) 的形式计算。

## 转换曲线

中心两分支相等当且仅当

\[
e^{-\ell/\beta}+e^{\ell/(1-\beta)}=2,
\qquad \ell=-\log(1-t).
\]

左边关于 \(\beta\) 严格递增。在 \(0<t<1/2\) 时，其在 \(0+\) 的极限小于
2，而在 \(1-\) 发散，所以根唯一。令
\(\beta=1/2+p_1\ell+p_2\ell^2+p_3\ell^3+p_4\ell^4+\cdots\)，
逐阶比较给出

\[
p_1=-\frac12,\quad p_2=0,\quad p_3=-\frac16,\quad p_4=0.
\]

从而

\[
\beta_*(t)=\frac12-\frac t2-\frac{t^2}{4}-\frac{t^3}{3}
-\frac{3t^4}{8}+O(t^5).
\]

这也说明临界曲线向 \(\beta<1/2\) 一侧移动；在精确临界点 \(\beta=1/2\)，
外部捕获分支已经在二阶占优。

## 实际符号和数值复核

- `verification/one_dimensional.wl` 以 Wolfram Language 15 对系数和临界窗口展开作符号复算；
- `verification/one_dimensional.py` 以 80 位 `mpmath` 和独立的 SciPy `brentq` 求根比较；
- Python 输出中的 ``abs(mp-scipy)`` 是两种数值方法的一致性检查，不能代替唯一性证明。
