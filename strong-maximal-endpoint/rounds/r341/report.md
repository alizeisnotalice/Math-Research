# 第341轮　临界抽稀下全部中央门槛的统一同源费用

- **A 当前覆盖。** 任意可测等权深度轮廓，中央误差有界（第325轮）。
- **本轮缩小的缺口。** 第340轮的临界抽稀模型在全部 $p\in I\Subset(0,1)$ 上有同时有效的双侧同源质量估计，包含随维数移动的门槛。给出可用确定性标签对交积检验的下界。
- **下一条未证不等式。** 任意有限成功盒族需要在缺少二项层结构时也能控制高重数的有效命中质量，达到 $B_n(p)+o_I(\log n)$。
- **一般系数差。** 仍为 $(\mu-1)/8$。

沿用[第340轮](../r340/report.md)的 $a=\alpha n,B=\alpha+(1-\alpha)l,s^2=(1-\alpha)l(1-l)$、$j_0=\lfloor nB+ts\sqrt n\rfloor$、$\theta=1/\binom{j_0}a$、$\lambda_j=\theta\binom ja$；$W$ 为实际覆盖重数，$V=\int W$。令 $j_p=\lfloor nB+s\sqrt n\,\Phi^{-1}(1-p)\rfloor$。以概率趋于一，同一个保留族对全部 $p\in I$ 有
\[
 \frac{m(p)}V\asymp_{I,\alpha,l,t}
 \begin{cases}
 1/(\sqrt n\lambda_{j_p}),&j_p\ge j_0,\\
 (1+j_0-j_p)/\sqrt n,&j_p<j_0.
 \end{cases}
\]
特别对固定 $p<1-\Phi(t)$，实际接收价格满足
\[
 \log(|U|/m(p))=h'(B)s[\Phi^{-1}(1-p)-t]\sqrt n+\tfrac12\log n+O(1),
 \qquad h'(B)=\log(B/(B-\alpha)).
\]
因为 $|U|/V\to\Phi(t)>0$（依概率），且相邻 $\lambda$ 的比为 $j/(j-a)$，在中央层取对数并相加只产生 $O(1)$ 的线性化误差。

**确定性上界和下界接口。** 记 $H=a+\operatorname{Bin}(n-a,l)$、$V_0=\binom na\theta l^a$。每件实际盒中 $J$ 的律都为 $H$。由 Berry–Esseen、中央正态密度下界，可固定整数 $K$ 和 $c>0$，同时有 $\Pr(H\ge j_p-K)>p$ 和 $\Pr(H>j_p+K)\le p-c/\sqrt n$。用 $U\cap\{J\ge j_p-K\}$ 作一个源即可严格成功。$J=j$ 层的均匀体积恰为 $V_0\Pr(H=j)/\lambda_j$（$j\ge a$），该层实际覆盖体积还至多为 $V\Pr(H=j)$。以上分别从环境体积和盒体积之和计价；单层概率最多 $C/\sqrt n$。当 $V\asymp V_0$ 时，$j\ge j_0$ 按 $1/\lambda_j$ 几何求和，其下按层数计费，就得到所列两种上界。$K$ 个相邻中央层的 $\lambda$ 比值均有统一常数界。

对下界，令 $j=j_p+K$。任何成功源都须在 $J\le j$ 支付至少 $cV/\sqrt n$ 的重数积分。若某整数 $D\ge1$ 满足
\[
 \int_{J\le j}(W-D)_+\le cV/(2\sqrt n),
\]
则必有 $m(p)\ge cV/(2D\sqrt n)$。这条接口适用于任意确定性保留族。左端亦可由 $D^{-1}\int_{J\le j}W(W-1)$ 支配，后一积分精确等于对不同保留短边标签 $A,A'$ 的有序对求和，项为 $l^b\Pr(b+\mathrm{Bin}(n-b,l)\le j)$，$b=|A\cup A'|$。故无需假设各盒对同一源的命中彼此独立。

**随机族同时满足溢出界。** 对每个中央整数 $j$ 取 $D_j=\lceil C_*\max(1,\lambda_j)\rceil$。如第340轮，溢出期望除以 $V_0$ 至多
\[
 D_j^{-1}\sum_{i\le j}\Pr(H=i)\lambda_i\le C/(C_*\sqrt n).
\]
选 $C_*$ 使此项至多 $c/(8\sqrt n)$。改变一个标签的溢出量变化至多盒体积 $l^a$，与 $D_j$ 无关；相对 $V_0$ 的方差至多 $1/(\theta\binom na)$，指数小。Chebyshev 后对至多 $n+1$ 个整数层取并，同时 $V/V_0=1+o(1)$，得所有所需下界 $c'/[\sqrt n\max(1,\lambda_{j_p})]$。

**补足 $j_p<j_0$ 时的线性层数。** 任一源有 $m(p)\ge pV-V+|U|$。在此模型 $|U|/V=\Pr(H\le j_0)+O(1/\sqrt n)$ 以概率趋一：第340轮期望换元的系数在 $i\le j_0$ 与一之差至多 $\lambda_i$（二项并界的二阶 Bonferroni），在 $i>j_0$ 至多 $1/\lambda_i$。两端几何求和并用 $\sup_i\Pr(H=i)\le C/\sqrt n$；随机波动仍指数小。故 $m(p)/V\ge p-\Pr(H>j_0)-O(1/\sqrt n)$。在 $j_p$ 至 $j_0$ 的固定标准化窗口，每个 $H$ 格点质量至少 $c_I/\sqrt n$（Stirling），且 $p=\Pr(H>j_p)+O_I(1/\sqrt n)$。结合上一段已有的 $c'/\sqrt n$ 即得完整下界。[精确检验](check.py)直接比较一族非随机小盒的截断重数对积分与标签对公式。

下一步应从确定性标签交积中提取足够的局部重数界，或以实际源的高重数份额取代二阶上界；此处的全门槛双侧式依赖上述随机族。
