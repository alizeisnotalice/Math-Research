# r736：固定共同源给出的负责 LF 面与容量顶点表示

本稿只处理有限真实矩形目录的精确外层表示，不给中央容量估计，顶点枚举与容量证书见同目录报告。记 $X=\bigcup_{i=1}^kR_i$，其非空成员原子为 $A_j$，体积 $m_j>0$，$\Omega=\sum_jm_j$，$v_i=|R_i|>0$。固定 $0<q<1$。负责权多面体为
\[
\Pi=\left\{\pi\ge0:\sum_i\pi_i=1,\quad
\sum_{i\in I}\pi_i\le\frac{|\bigcup_{i\in I}R_i|}{\Omega}
\quad(I\subseteq\{1,\ldots,k\})\right\}.
\]
令 $d_j(\pi)=\sum_{i:A_j\subset R_i}\pi_i/v_i$，并定义真实 NP 价格
\[
b_q(\pi)=\min\left\{\sum_jm_jg_j:
0\le g_j\le1,\quad\sum_jm_jd_j(\pi)g_j\ge q\right\}.
\]
写 $s=\max_{\pi\in\Pi}b_q(\pi)$，$\Pi_{\rm LF}=\{\pi\in\Pi:b_q(\pi)=s\}$。

取 r716/r732 下包络源端对偶给出的一个**共同最优源** $f=(f_j)$：
\[
0\le f_j\le1,\qquad \sum_jm_jf_j=s,\qquad
\min_{\pi\in\Pi}\sum_i\pi_i a_i(f)\ge q,\qquad
a_i(f)=v_i^{-1}\sum_{j:A_j\subset R_i}m_jf_j. \tag{1}
\]
若使用证书变量 $\ell_j$，条件 $\ell_j\le a_i(f)$ 对每个含 $A_j$ 的行、以及 $\sum_jm_j\ell_j\ge q\Omega$，通过负责运输分割保证 (1) 的最后一个不等式。以下不要求共同源唯一，也不要求其原子坐标均为零或一。

## 1. LF 权的线性充要条件

置 $D_{\max}=\max_i v_i^{-1}$。定义 $(\pi,\eta)$ 的有界多面体
\[
\mathcal P_f=\left\{\begin{array}{l}
\pi\in\Pi,\quad 0\le\eta\le D_{\max},\quad
\sum_i\pi_i a_i(f)=q,\\
d_j(\pi)\le\eta\quad(f_j=0),\\
d_j(\pi)\ge\eta\quad(f_j=1),\\
d_j(\pi)=\eta\quad(0<f_j<1).
\end{array}\right\}. \tag{2}
\]
则有精确投影表示
\[
\boxed{\qquad \Pi_{\rm LF}=\operatorname{proj}_{\pi}\mathcal P_f.\qquad} \tag{3}
\]

**必要性。** 对任意 LF 权，(1) 给出 $P_\pi(f)\ge q$，而源预算是 $s=b_q(\pi)$。若命中率严格大于 $q$，把 $f$ 乘以 $q/P_\pi(f)<1$，得到预算严格小于 $s$ 的可行 NP 源，矛盾。因此 $P_\pi(f)=q$，且 $f$ 为该权的 NP 最优检验。这就是任意 LF 权与任意共同最优源的交叉鞍点性质。

有限 NP 检验按密度 $d_j$ 从高到低填充：若某较低密度原子已取正量，而较高密度原子尚未取满，把等量空间体积从前者移向后者会增加命中率；再把所得源整体略缩小即降低预算。故 $f$ 的非零部分不能位于任何尚未填满的严格更高密度原子之下。等价地存在阈值 $\eta$ 满足 (2)。所有密度在 $[0,D_{\max}]$ 内，故阈值可选在同一区间。该论证允许密度相等、源分数坐标、零密度原子及退化阈值。

**充分性。** 设 (2) 成立。对任意 $0\le g\le1$，逐原子有
\[
(d_j-\eta)(g_j-f_j)\le0,
\]
故
\[
P_\pi(g)-q\le\eta\left(\sum_jm_jg_j-s\right). \tag{4}
\]
事实上满足 (2) 的阈值必有 $\eta>0$：若 $\eta=0$，任何 $d_j>0$ 的原子只能有 $f_j=1$，于是 $P_\pi(f)=\sum_jm_jd_j=1$，与 $q<1$ 矛盾。因而任何命中率至少为 $q$ 的 $g$，由 (4) 都有预算至少为 $s$。同时 $f$ 达到该预算，故 $b_q(\pi)=s$，即 $\pi$ 为 LF 权。

所以 (3) 是线性充要表示，而非只列必要 KKT 条件。当原子体积、行体积与所选共同源 $f$ 都为有理数时，(2) 及其投影均为有理多面体；一般实参数时为实多面体。无需假设 $\eta$ 唯一，也无需把所有分数原子误作独立阈值。

## 2. 几何容量的凸性、连续性与顶点穷尽

令 $Q=\operatorname{Unif}(X)$。对固定标签权定义
\[
C(\pi)=e^{\Lambda(\pi)}
=\sup_{F\ge0,\ QF=1}\sum_i\pi_i
\exp\left(\frac1{v_i}\int_{R_i}\log F\right). \tag{5}
\]
零值按几何平均的连续延拓解释。把 $F$ 在每个成员原子上替换成算术平均，不改变 $QF$，并由对数 Jensen 不减任何行几何平均。因此 (5) 等价于有限维最大值
\[
C(\pi)=\max_{x_j\ge0,\ \sum_jm_jx_j=\Omega}
\sum_i\pi_i\prod_{j:A_j\subset R_i}x_j^{m_j/v_i}. \tag{6}
\]
可行域紧，每个乘积连续（各出现的指数严格正），故最大值确实取得。若原定义仅允许严格正函数，令 $x\mapsto(1-\varepsilon)x+\varepsilon\boldsymbol1$ 即得相同上确界。

(6) 是线性于 $\pi$ 的函数族的上确界，故 $C$ 凸。每个行几何平均由算术平均控制为至多 $\Omega/v_i$，从而
\[
|C(\pi)-C(\pi')|\le\Omega D_{\max}\|\pi-\pi'\|_1. \tag{7}
\]
特别地，即使部分标签权为零，$C$ 仍连续；并有 $1\le C\le\Omega D_{\max}$。不需要断言 $\Lambda$ 本身凸。

有限多面体中每一点是顶点的凸组合。由凸性，任一点的 $C$ 不超过其组合中最大的顶点值；因此
\[
\max_{\pi\in\Pi_{\rm LF}}\Lambda(\pi)
=\log\max_{\pi\in\operatorname{vert}(\Pi_{\rm LF})}C(\pi)
=\log\max_{(\pi,\eta)\in\operatorname{vert}(\mathcal P_f)}C(\pi). \tag{8}
\]
最后一式也可直接对扩展多面体上的凸连续函数 $(\pi,\eta)\mapsto C(\pi)$ 应用同一论证。扩展顶点的投影不必本身是投影多面体的顶点，且可能重复；全部扩展顶点的投影仍足以穷尽外层最大值。式 (8) 只说至少一个顶点取最大，不说全部最大点必须是顶点。

由此，固定 $\pi$ 的 r731/r732 内层容量证书可以与精确 LF 多面体的有限顶点表示衔接。本轮另以精确枚举与独立重建实现有限目录外层认证；近中央统一估计仍未证明。

**本机小例核验。** 对真实区间 $R_1=[0,1]$、$R_2=[0,2]$，取 $q=1/2$，两个单位成员原子上共同源均为 $f_j=1/2$、预算 $s=1$。写 $\pi=(a,1-a)$，Hall 为 $0\le a\le1/2$，两密度为 $a+(1-a)/2$ 与 $(1-a)/2$。本机 Wolfram 精确执行 `Reduce[0<=a<=1/2 && 0<=eta<=1 && a+(1-a)/2==eta && (1-a)/2==eta,{a,eta},Reals]`，实际返回 `a == 0 && eta == 1/2`，与唯一 LF 权全部放在大区间一致。此检查没有枚举一般目录的顶点。
