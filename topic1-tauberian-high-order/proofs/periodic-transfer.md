# 有界半径周期见证向有限测度输入的定量传递

本引理为周期模型认证提供必要的有限测度桥梁；它不证明候选模型的捕获比例，也不处理没有半径上界的数值搜索。

令 $D=[-1/2,1/2)^n$，$A\subset D$ 可测，$a=|A|<1$，并令
\[
 E_\infty=\mathbb R^n\setminus\bigcup_{k\in\mathbb Z^n}(k+A).
\]
设 $0<t<1$，$B\subset A$ 为可测集，$b=|B|$，每个 $x\in B$ 都有 $0<r_x\le R<\infty$ 满足
\[
 |E_\infty\cap Q(x,r_x)|>(1-t)|Q(x,r_x)|.
\]
忽略零测集不改变以下结论。取整数 $q,N$，满足 $q\ge R$ 和 $N>2q$，定义
\[
 D_N=[-1/2,N-1/2]^n,\qquad E_N=E_\infty\cap D_N.
\]
则 $|E_N|=N^n(1-a)$，且
\[
 \frac{|\{M_c\mathbf1_{E_N}>1-t\}|}{|E_N|}
 \ge 1+\left(1-\frac{2q}{N}\right)^n\frac b{1-a}. \tag{T1}
\]

证明：对 $k\in\{q,\ldots,N-q-1\}^n$ 与 $x\in B$，$k+x$ 到 $D_N$ 每个面的距离至少为 $q$，所以全部见证方体包含在 $D_N$ 内。周期性使 $k+B$ 上的严格平均不变。这些集合互不交且位于 $E_N$ 的补集，总体积为 $(N-2q)^n b$。Lebesgue 密度定理又给出 $E_N$ 几乎处处包含在该严格水平集。两项相加即得 (T1)。

由 Bernoulli 不等式，截断损失至多
\[
 \left[1-\left(1-\frac{2q}{N}\right)^n\right]\frac b{1-a}
 \le\frac{2nq}{N}\frac b{1-a}. \tag{T2}
\]
因此 $C_{c,n}(1-t)\ge1+b/(1-a)$。若 $b/(1-a)=O_n(t)$、$q=O_n(1)$，要使截断误差为 $o(t^{1+1/n})$，可以取 $N(t)t^{1/n}\to\infty$。给定允许误差 $\eta>0$，显式选择
\[
 N>\max\left\{2q,\frac{2nqb}{(1-a)\eta}\right\}
\]
即可使截断误差小于 $\eta$。

边界说明：半开基本胞与闭截断盒仅相差零测面；严格平均的裕量不因截断降低，因为所有已用见证完整留在盒内。本引理不把采样点比例升级为 $b$ 的严格体积下界。数值采样不能提供其所需的可测捕获区域证明。
