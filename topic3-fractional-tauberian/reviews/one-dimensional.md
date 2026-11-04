# 一维精细渐近独立复核

审查对象：`paper/sections/06-one-dimensional.tex` 与
`proofs/one-dimensional-refinements.md`。本审查核对闭式适用域、级数余项的参数一致性、中心两分支的转换曲线及临界窗口；不重写原章节。

## 结论

一维非中心与中心闭式的证明路线正确；分支方程、转换曲线前四个 $t$ 系数和临界窗口的四阶系数均正确。Wolfram Language 15 与 SymPy 1.14 分别实际展开，给出

\[
\begin{aligned}
A_{1/2+\lambda t}(t)
 &=2t+(-4\lambda-1)t^2+(8\lambda^2+6\lambda)t^3\\
 &\quad+\left(-16\lambda^3-20\lambda^2-\frac{4\lambda}{3}\right)t^4+O(t^5),\\
D_{1/2+\lambda t}(t)
 &=2t+(4\lambda+3)t^2+(8\lambda^2+10\lambda+4)t^3\\
 &\quad+\left(16\lambda^3+28\lambda^2+\frac{52\lambda}{3}+5\right)t^4+O(t^5).
\end{aligned}
\]

并且令 $\ell=-\log(1-t)$ 后，转换根满足
\[
\beta_*(\ell)=\frac12-\frac\ell2-\frac{\ell^3}{6}+O(\ell^5),
\]
故
\[
\beta_*(t)=\frac12-\frac t2-\frac{t^2}{4}-\frac{t^3}{3}-\frac{3t^4}{8}+O(t^5).
\]
高精度数值求根与该四阶展开一致：$t=10^{-2}$ 的绝对误差约为 $5.9\times10^{-11}$，$t=10^{-1}$ 时约为 $6.8\times10^{-6}$。这是符号推导的交叉检查，未作为证明依据。

## 必须修正

1. **全阶命题的常数闭式缺少有效域。** 命题 `one:all-orders` 中 $A_\beta,D_\beta$ 的 Taylor 展开及其余项，对 $K\Subset(0,1)$ 与任意 $0<t\le t_0<1$ 都正确；但紧随其后的
   \[
   C_{\beta,\mathcal Q}(1-t)=1+2\sum_{j=1}^N\frac{(q)_j}{j!}t^j+O(t^{N+1})
   \]
   调用了定理 `one:uncentered`，后者还要求
   \(1-t>2^{-\beta}\)。为使结论关于 $\beta\in K=[\beta_-,\beta_+]$ 一致，应把该常数展开和“中心常数为两级数最大值”限定为
   \[
   0<t\le t_0<1-2^{-\beta_-}.
   \]
   或在命题中明确区分“代数函数的全域展开”和“由闭式导出的 Tauberian 常数展开”。

2. **临界窗口的常数等式需要补一个最终有效域。** 命题 `one:critical-window` 的渐近系数及 $O_\Lambda(t^5)$ 一致性正确，但
   \[
   C_{1/2+\lambda t,\mathcal Q_c}(1-t)=\max\{A_{1/2+\lambda t}(t),D_{1/2+\lambda t}(t)\}
   \]
   应表述为：对每个 $\Lambda>0$，存在 $t_\Lambda>0$，使 $0<t<t_\Lambda$ 和 $|\lambda|\le\Lambda$ 时成立。这里 $t_\Lambda$ 同时保证 $1/2+\lambda t\in(0,1)$ 及 $1-t>2^{-(1/2+\lambda t)}$。目前的“当 $t\downarrow0$”表达可被理解为这一含义，但把存在量词写出才能使闭式的调用无歧义。

3. **生产排版缺陷。** `one:pairwise`、`one:centered` 等处有字面串 `qquad`，缺少反斜杠；这会在 PDF 中输出错误文本。它不影响推导，但应在编译前修复。

## 已核对的逻辑点

- 在 $\alpha>2^{-\beta}$ 下，严格见证给出 $s_I>\alpha^{1/\beta}>1/2$，因此所有见证区间两两相交；再通过任意两区间并的长度界控制整个光晕，量词有效。
- 转换函数
  \[
  F_\ell(\beta)=e^{-\ell/\beta}+e^{\ell/(1-\beta)}-2
  \]
  的导数严格为正；在 $0<t<1/2$ 时，$F_\ell(0+)=e^\ell-2<0$ 且 $F_\ell(1-)=+\infty$，所以代数转换根存在且唯一。它成为 Tauberian 常数的分支转换还须添加上述闭式有效域。
- 余项常数 $M^A_{K,N,t_0}$ 与 $M^D_{K,N,t_0}$ 确实控制 $N+1$ 阶导数；后者取 $\beta_+$ 给出对所有 $\beta\in K$ 的统一上界。
- 在 $\beta=1/2$，外部捕获分支为 $(1-t)^{-2}-1=2t+3t^2+4t^3+\cdots$，故在一阶并列后，二阶由外部机制严格占优。
