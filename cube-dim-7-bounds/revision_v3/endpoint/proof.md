# 第三次修订：参数端点、内部极小值存在性与维数单调性

## 范围与证明状态

本附件只分析所给的一参数上界，不使用旧下界，也不进行网格搜索。
下列结论均有解析证明：

* `U_1=2, U_2=4, U_3=8`，三者均不在有限参数处取得，只能由 `t↓1` 逼近。
* 对每个整数 `n≥4`，`U_n<2^n`，而且全局极小值在某个有限 `t>1` 处取得。
* 对每个 `n≥4`，该函数至少有一个局部极大点和一个全局极小点。因此，不能声称只有一个驻点。
* 中心立方体 Hardy–Littlewood 极大算子的最优弱 `(1,1)` 常数满足 `C_n≤C_{n+1}`。

**未证明：**`n≥4` 的全局极小点唯一性；本附件不把内部极小值存在性表述成唯一性。

## 1. 定义与一般性质

设整数 `n≥1`，`t>1`，并定义

\[
A_n(t)=(1+t^{-1})^n,
\]
\[
B_n(t)=1-\left(\frac{1+t^{-1}}2\right)^n+n\log t
-\frac n2\int_1^t
\left(1+\frac1{2t}-\frac1{2\rho}\right)^{n-1}\frac{d\rho}{\rho},
\]
\[
F_n(t)=(\sqrt{A_n(t)}+\sqrt{B_n(t)})^2,
\qquad U_n=\inf_{t>1}F_n(t).
\]

对 `1≤ρ≤t`，积分中括号的值属于 `[(1+t^{-1})/2,1]`。
因此积分非负且不超过 `log t`，从而

\[
B_n(t)\ge 1-\left(\frac{1+t^{-1}}2\right)^n+
\frac n2\log t>0,\qquad B_n(t)<1+n\log t.
\tag{1}
\]

于是 `F_n` 在 `(1,∞)` 光滑，连续延拓到 `t=1` 后有
`F_n(1)=2^n`；同时 `(1)` 给出 `F_n(t)→∞` 当 `t→∞`。

直接求导，在 `t=1` 处积分的变上限项等于 `1`，参数导数项的积分区间退化，故

\[
B_n(1)=0,\qquad B_n'(1)=\frac n2+n-\frac n2=n.
\]

于是，当 `ε↓0` 时

\[
F_n(1+\varepsilon)=2^n+2\sqrt{n2^n}\sqrt\varepsilon+O(\varepsilon).
\tag{2}
\]

特别地，右端点附近有 `F_n(t)>2^n`。这里 `B_n'(1)` 指公式在闭端点的右导数。

## 2. 一维、二维、三维的端点极小值

令 `x=1/t∈(0,1)`，`d=1−x`，`L=−log x>0`。以下 `B_n` 均在 `t=1/x` 取值。

### 一维

\[
B_1=\frac{d+L}{2}>\frac d2.
\]

而

\[
(\sqrt2-\sqrt{1+x})^2
=\frac{d^2}{(\sqrt2+\sqrt{1+x})^2}<\frac d2.
\]

所以 `F_1(t)>2` 对每个 `t>1` 成立，结合端点极限即 `U_1=2`。

### 二维

精确积分给出

\[
B_2=\frac54-x-\frac{x^2}4+\left(1-\frac x2\right)L.
\]

因此

\[
B_2-d^2=\frac{d(5x+1)}4+\left(1-\frac x2\right)L>0.
\]

由于 `√A_2=1+x`，所以 `√A_2+√B_2>(1+x)+d=2`。
这证明 `F_2(t)>4` 对全部 `t>1` 成立，故 `U_2=4`。

### 三维

精确积分给出

\[
B_3=\frac{35-18x-15x^2-2x^3+(24-24x-6x^2)L}{16}.
\]

我们证明更强的下界

\[
B_3>\frac92(1-x)^2.\tag{3}
\]

记

\[
P=-37+126x-87x^2-2x^3,\qquad Q=24-24x-6x^2.
\]

则 `16(B_3−(9/2)d²)=P+QL`。使用对全部 `0<x<1` 成立的基本对数界

\[
\frac{2(1-x)}{1+x}\le-\log x\le\frac{1-x}{x}.\tag{4}
\]

左边可由 `h(x)=−log x−2(1−x)/(1+x)` 满足 `h(1)=0`、
`h'(x)=−(1−x)²/[x(1+x)²]≤0` 得出；右边直接由
`−log x=∫_x^1 du/u≤(1−x)/x` 得出。

当 `Q≥0` 时，`(4)` 左侧给出

\[
P+QL\ge P+Q\frac{2d}{1+x}
=\frac{d(2x^3+79x^2+4x+11)}{1+x}>0.
\]

当 `Q<0` 时，`(4)` 右侧给出

\[
P+QL\ge P+Q\frac d{x}
=\frac{d(2x^3+83x^2-61x+24)}x>0.
\]

最后的严格正性由 `x>0` 及二次式 `83x²−61x+24` 的正首项系数和负判别式
`61²−4·83·24=−4247<0` 保证。由此 `(3)` 得证。

另一方面，由积分形式的均值估计

\[
2\sqrt2-(1+x)^{3/2}
=\int_{1+x}^2\frac32\sqrt u\,du
\le\frac3{\sqrt2}(1-x).
\]

结合 `(3)`，得到 `√B_3>2√2−(1+x)^{3/2}`，故 `F_3(t)>8`，即 `U_3=8`。

上述三维数下，严格不等式、`t→∞` 时发散以及紧区间连续性还说明：
若 `F_n(t_j)→2^n`，则必有 `t_j→1`；不存在另外的有限或无穷端参数极小序列。

## 3. 每个 n≥4 均存在内部全局极小点

只需证明同一个参数 `t=2` 对所有 `n≥4` 满足 `F_n(2)<2^n`。

当 `n=4` 时，精确积分得

\[
A_4(2)=\frac{81}{16},\qquad
B_4(2)=\frac{1841+72\log2}{768}
<\frac{1913}{768}<\frac{49}{16}.
\]

这里仅用了 `log2<1`。所以
`√A_4(2)+√B_4(2)<9/4+7/4=4`，即 `F_4(2)<16`。

对 `n≥5`，记

\[
\delta_n=2^{n/2}-(3/2)^{n/2}>0,\qquad \ell=\log2.
\]

由 `(1)`，`B_n(2)<1+nℓ`。而

\[
\delta_5^2=\frac{1267}{32}-18\sqrt3
>\frac{1267}{32}-18\cdot\frac74
=\frac{259}{32}>6>1+5\ell.
\]

其中 `√3<7/4` 可直接平方核对。又有

\[
\delta_{n+1}=\sqrt2\,\delta_n+
(\sqrt2-\sqrt{3/2})(3/2)^{n/2}>\sqrt2\,\delta_n,
\]

以及 `2(1+nℓ)>1+(n+1)ℓ`。归纳得到
`δ_n²>1+nℓ` 对全部 `n≥5` 成立。因此

\[
\sqrt{A_n(2)}+\sqrt{B_n(2)}
<(3/2)^{n/2}+\sqrt{1+n\ell}<2^{n/2},
\]

即所需的严格见证。由于 `F_n(1)=2^n>F_n(2)`，`F_n(t)→∞` 当 `t→∞`，
可先排除 `1` 的充分小邻域与充分大的尾部，再在余下紧区间取最小值。
这证明对每个 `n≥4`，`U_n<2^n` 且存在 `t_n∈(1,∞)` 满足 `F_n(t_n)=U_n`。

由 `(2)`，可取 `s∈(1,2)` 使 `F_n(s)>2^n`。延拓后的 `F_n` 在 `[1,2]`
有最大值，而最大值严格大于两个端点值，故某个最大值点位于 `(1,2)`，是驻点。
它的值大于 `2^n`，与值小于 `2^n` 的全局极小点不同。因此，每个 `n≥4`
至少有两个不同驻点。这一论证不涉及全局极小点的唯一性。

## 4. 中心立方体弱 (1,1) 常数的维数单调性

采用半边长参数，定义

\[
M_ng(x)=\sup_{r>0}\frac1{(2r)^n}\int_{[-r,r]^n}|g(x+y)|\,dy,
\]

并令 `C_n` 为使 `λ|{M_ng>λ}|≤C_n‖g‖_1` 对所有 `g∈L¹(Rⁿ)`、`λ>0`
成立的最优常数。若 `C_{n+1}=∞`，待证单调性显然成立。否则，先取维数 `n+1`
的任意可用弱型常数 `K`，最后对 `K` 取下确界。

固定任意 `g∈L¹(Rⁿ)`、`r_0>0` 和 `R>r_0`，定义

\[
f_R(x,z)=g(x)\mathbf1_{[-R,R]}(z),\qquad
M_n^{\le r_0}g(x)=\sup_{0<r\le r_0}
\frac1{(2r)^n}\int_{[-r,r]^n}|g(x+y)|\,dy.
\]

`f_R∈L¹(R^{n+1})`，且 `‖f_R‖_1=2R‖g‖_1`。若 `|z|≤R−r_0` 且 `r≤r_0`，
则 `[z−r,z+r]⊂[-R,R]`，从而严格的乘积积分恒等式给出

\[
\frac1{(2r)^{n+1}}\int_{[-r,r]^{n+1}}
|f_R(x+y,z+w)|\,dy\,dw
=\frac1{(2r)^n}\int_{[-r,r]^n}|g(x+y)|\,dy.
\]

因此

\[
\{M_n^{\le r_0}g>\lambda\}\times[-R+r_0,R-r_0]
\subseteq\{M_{n+1}f_R>\lambda\}.
\]

由维数 `n+1` 的弱型估计及 Tonelli 定理，

\[
\lambda\,2(R-r_0)|\{M_n^{\le r_0}g>\lambda\}|
\le K\,2R\|g\|_1.
\]

固定 `r_0` 后令 `R→∞`，得到
`λ|{M_n^{≤r_0}g>λ}|≤K‖g‖_1`。随后令 `r_0↑∞`；截断极大函数单调上升到
包含**全部尺度**的 `M_ng`，其严格超水平集单调增至 `{M_ng>λ}`，由测度从下连续性
得到 `λ|{M_ng>λ}|≤K‖g‖_1`。对可用 `K` 取下确界即 `C_n≤C_{n+1}`。

这个证明对一般 `L¹` 输入成立，无需紧支集、光滑性或固定支撑尺度。
极限顺序是先 `R→∞`，再 `r_0→∞`，从而避免把未截断的全部半径放进有限内区。

## 5. 实际运行的精确核验与复现

配套 `symbolic_check.py` 以本机 SymPy 对有限二项式展开作精确积分，核对：
`B_2,B_3,B_4(2)` 的精确式、三维分情况的两个因式分解、判别式、`n≥4` 见证中的
有理数/根式恒等式，以及 `n=1,…,8` 的端点值和端点导数。后者只是对一般解析导数论证的额外核验。

实际执行：

```sh
PYTHONPATH=/Users/zhengzhihao/Documents/ChatGPT/bnupt/tmp/y_pendulum_full/deps \
python3 cube-dim-7-bounds/revision_v3/endpoint/symbolic_check.py
```

输出：

```text
PASS: exact B2/B3/B4 formulas, logarithm-bound factorizations,
      discriminant, n>=4 witness certificates, B(1), and Bprime(1).
No grid search or floating-point certification was used.
```

普通已安装 SymPy 的环境可直接运行 `python3 symbolic_check.py`。
本机 MATLAB 调用也曾尝试，但连接返回“无法解析名称 `jupyter.Idler`”，未取得计算结果；
因此不把 MATLAB 当作本附件的验证证据。所用证明本身是上述解析不等式，精确符号运行只核验代数恒等式。

## 6. 可直接纳入论文的英文 TeX 片段

以下片段假定正文已经定义 `A_n,B_n,F_n,U_n`，并已加载 `amsmath,amssymb`。

```tex
\begin{proposition}\label{prop:parameter-endpoints}
For the parameter bound $U_n=\inf_{t>1}F_n(t)$, one has
$U_1=2$, $U_2=4$, and $U_3=8$. In each of these dimensions,
$F_n(t)>2^n$ for every $t>1$, and the infimum is approached
only as $t\downarrow1$.
For every integer $n\ge4$, one has $F_n(2)<2^n$, and $U_n$
is attained at some finite parameter $t_n>1$.
Moreover, $F_n$ has at least two distinct stationary points when
$n\ge4$. These assertions do not establish uniqueness of a global
minimizer or certify any particular numerical candidate as a global minimizer.
\end{proposition}

\begin{proof}
The factor in the integrand defining $B_n(t)$ belongs to
$[(1+t^{-1})/2,1]$ for $1\le\rho\le t$. Consequently,
\begin{equation}\label{eq:Bn-elementary-bounds}
0<1-\left(\frac{1+t^{-1}}2\right)^n+
\frac n2\log t\le B_n(t)<1+n\log t.
\end{equation}
Thus $F_n$ is smooth on $(1,\infty)$, extends continuously to
$F_n(1)=2^n$, and tends to infinity as $t\to\infty$.
Set $x=t^{-1}\in(0,1)$, $d=1-x$, and $L=-\log x$.
Direct integration gives
\begin{align*}
B_1&=\frac{d+L}{2},\\
B_2&=\frac54-x-\frac{x^2}{4}
       +\left(1-\frac x2\right)L,\\
16B_3&=35-18x-15x^2-2x^3+(24-24x-6x^2)L.
\end{align*}
Here and below $B_n$ in formulas involving $x$ means $B_n(1/x)$.
In dimension one,
\[
B_1>\frac d2>
\frac{d^2}{(\sqrt2+\sqrt{1+x})^2}
=(\sqrt2-\sqrt{1+x})^2.
\]
In dimension two,
\[
B_2-d^2=\frac{d(5x+1)}4+
\left(1-\frac x2\right)L>0.
\]
These inequalities imply $F_1>2$ and $F_2>4$.

For dimension three, write
\[
16\left(B_3-\frac92d^2\right)=P+QL,
\quad
P=-37+126x-87x^2-2x^3,
\quad Q=24-24x-6x^2.
\]
The elementary logarithmic bounds
\[
\frac{2d}{1+x}\le L\le\frac d{x}
\]
follow, respectively, by differentiating
$-\log x-2(1-x)/(1+x)$ and by integrating $1/u\le1/x$
over $[x,1]$. If $Q\ge0$, then
\[
P+QL\ge\frac{d(2x^3+79x^2+4x+11)}{1+x}>0.
\]
If $Q<0$, then
\[
P+QL\ge\frac{d(2x^3+83x^2-61x+24)}{x}>0,
\]
where $83x^2-61x+24>0$ because its discriminant is
$61^2-4\cdot83\cdot24=-4247<0$.
Thus $B_3>(9/2)d^2$. Since
\[
2\sqrt2-(1+x)^{3/2}
=\int_{1+x}^2\frac32\sqrt u\,du
\le\frac3{\sqrt2}d,
\]
we obtain $F_3>8$. Together with $F_n(1)=2^n$, this proves
$U_n=2^n$ for $n=1,2,3$. Strict inequality on compact
subintervals of $(1,\infty)$ and divergence at infinity show that
any minimizing sequence must tend to $1$.

At $t=2$, direct integration in dimension four yields
\[
A_4(2)=\frac{81}{16},\qquad
B_4(2)=\frac{1841+72\log2}{768}
<\frac{1913}{768}<\frac{49}{16},
\]
so $F_4(2)<(9/4+7/4)^2=16$.
For $n\ge5$, let
$\delta_n=2^{n/2}-(3/2)^{n/2}$ and $\ell=\log2$.
The initial estimate is
\[
\delta_5^2=\frac{1267}{32}-18\sqrt3
>\frac{1267}{32}-18\cdot\frac74
=\frac{259}{32}>6>1+5\ell.
\]
Also,
\[
\delta_{n+1}=\sqrt2\,\delta_n+
(\sqrt2-\sqrt{3/2})(3/2)^{n/2}>\sqrt2\,\delta_n,
\qquad 2(1+n\ell)>1+(n+1)\ell.
\]
Induction gives $\delta_n>\sqrt{1+n\ell}$ for all $n\ge5$.
Using \eqref{eq:Bn-elementary-bounds}, we conclude that
\[
\sqrt{A_n(2)}+\sqrt{B_n(2)}
<(3/2)^{n/2}+\sqrt{1+n\ell}<2^{n/2}.
\]
Hence $F_n(2)<2^n$ for every $n\ge4$.
Continuity at $1$ and divergence at infinity now place a global
minimizer in a compact subinterval of $(1,\infty)$.

Finally, differentiating the defining formula at $t=1$ gives
$B_n'(1)=n$, and therefore
\[
F_n(1+\varepsilon)
=2^n+2\sqrt{n2^n}\sqrt\varepsilon+O(\varepsilon)
\qquad(\varepsilon\downarrow0).
\]
For $n\ge4$, the continuous extension of $F_n$ to $[1,2]$
therefore has an interior maximum with value greater than $2^n$.
This stationary point is distinct from a global minimizer, whose
value is less than $2^n$.
\end{proof}

\begin{proposition}\label{prop:dimension-monotonicity}
Let $C_n$ denote the optimal weak $(1,1)$ constant of the
centered maximal operator over axis-parallel cubes in $\mathbb R^n$.
Then $C_n\le C_{n+1}$.
\end{proposition}

\begin{proof}
There is nothing to prove if $C_{n+1}=\infty$.
Let $K$ be any admissible weak $(1,1)$ constant in dimension $n+1$.
Fix $g\in L^1(\mathbb R^n)$, $\lambda>0$, a scale cutoff $r_0>0$,
and $R>r_0$. Define
\[
f_R(x,z)=g(x)\mathbf1_{[-R,R]}(z),\qquad
M_n^{\le r_0}g(x)=\sup_{0<r\le r_0}
\frac1{(2r)^n}\int_{[-r,r]^n}|g(x+y)|\,dy.
\]
Then $\|f_R\|_1=2R\|g\|_1$.
If $|z|\le R-r_0$ and $0<r\le r_0$, the interval
$[z-r,z+r]$ lies in $[-R,R]$, so
\[
\frac1{(2r)^{n+1}}\int_{[-r,r]^{n+1}}
|f_R(x+y,z+w)|\,dy\,dw
=\frac1{(2r)^n}\int_{[-r,r]^n}|g(x+y)|\,dy.
\]
Consequently,
\[
\{M_n^{\le r_0}g>\lambda\}\times[-R+r_0,R-r_0]
\subseteq\{M_{n+1}f_R>\lambda\}.
\]
The dimension-$(n+1)$ weak estimate implies
\[
\lambda\,2(R-r_0)|\{M_n^{\le r_0}g>\lambda\}|
\le 2RK\|g\|_1.
\]
First let $R\to\infty$ with $r_0$ fixed. Then let
$r_0\uparrow\infty$; the truncated maximal functions increase
to the maximal function over all positive radii, and their strict
superlevel sets increase to $\{M_ng>\lambda\}$.
Continuity of measure from below yields
$\lambda|\{M_ng>\lambda\}|\le K\|g\|_1$.
Taking the infimum over admissible $K$ proves the claim.
The argument applies to arbitrary $L^1$ data, without a support
or regularity assumption.
\end{proof}
```
