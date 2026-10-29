# 无约束权弱型信息不能改进测度放大上界

固定 n≥1，m=2^n。E 为有限正测度可测集，F={M_c 1_E>1-t}\E。
以下仅否证一条丢掉高密度几何信息的证明路线，不是否证真实高阶上界。

## 2. 一个看似可行的转移公式

对 \(c>0\) 令
\[
 d\mu_c=(1+c\mathbf1_F)\,dx,\qquad f=\mathbf1_E.
\]
这是局部有限测度；不要求总测度有限。因为 \(E\cap F=\varnothing\)，有 \(\|f\|_{L^1(\mu_c)}=|E|\)。由 (1.1)，在 \(F\) 上
\[
 M_{c,\mu_c}f>
 \frac{1-t}{1+ct}=:\tau_{c,t}.\tag{2.1}
\]
确实，见证方体满足
\[
 \mu_c(Q)=|Q|+c|F\cap Q|<(1+ct)|Q|,
 \qquad \int_Q f\,d\mu_c=|E\cap Q|>(1-t)|Q|.
\]

定义无约束的密度比弱型常数
\[
 W_n(c):=\sup_{1\le w\le1+c}\|M_{c,w\,dx}\|_{L^1(wdx)\to L^{1,\infty}(wdx)}.
\tag{2.2}
\]
这里的上确界包括所有可测权、所有非零 \(L^1\) 输入与所有正阈值。用它估计 (2.1) 得
\[
 \frac{|F|}{|E|}
 \le W_n(c)\frac{1+ct}{(1+c)(1-t)}.\tag{2.3}
\]
若有 \(W_n(c)\le m-a_nc^{-1/(n+1)}\)，取 \(c\asymp t^{-(n+1)/n}\) 形式上会产生目标负修正。**这个候选稳定性不等式是假的。**下面给出定量反例，并证明 (2.3) 在此路径中完全不能改进 \(U_n\)。

## 3. 完整背景下的反例：\(W_n(c)\ge m-O_n(c^{-1})\)

### 命题 3.1

对每个 \(n\ge1\)、\(c>0\)、\(0<\delta<1/3\)，有
\[
 W_n(c)\ge
 B_{n,\delta}(c):=
 \frac{m(1+c)}{c+b_{n,\delta}},
 \qquad b_{n,\delta}=\left(\frac{1+\delta}{\delta}\right)^n.
\tag{3.1}
\]
因而
\[
 0\le m-W_n(c)
 \le\frac{m(b_{n,\delta}-1)}{c+b_{n,\delta}}
 \le\frac{m(b_{n,\delta}-1)}c.\tag{3.2}
\]
再令 \(\delta\uparrow1/3\)，得到更方便的公式
\[
 \boxed{\quad W_n(c)\ge \frac{m(1+c)}{c+4^n}.\quad}\tag{3.3}
\]

### 完整证明

令
\[
 A_v=v+[-\delta,\delta]^n\quad(v\in\{-1,1\}^n),\qquad
 A=\bigcup_vA_v,
\]
\[
 E_\varepsilon=[-\varepsilon,\varepsilon]^n,
 \quad 0<\varepsilon<1-3\delta,
 \quad R=1+\delta+\varepsilon,
 \quad w_c=1+c\mathbf1_A.
\]
缩小 \(\varepsilon\) 无损；上述条件已经保证 \(E_\varepsilon\cap A=\varnothing\)。记
\[
 V=(2\delta)^n,
 \qquad D_\varepsilon=(2R)^n+cV.
\]
对任意 \(x\in A_v\)，同一个半径 \(R\) 满足以下三项。

1. \(E_\varepsilon\subset Q(x,R)\)：逐坐标有 \(|x_i|+\varepsilon\le1+\delta+\varepsilon=R\)。
2. \(A_v\subset Q(x,R)\)：同块内点距至多 \(2\delta<R\)。
3. 若 \(u\ne v\)，存在一个坐标符号相反，于是在该坐标 \(A_v\) 与 \(A_u\) 的最小距离为 \(2-2\delta\)。而
   \[
   R=1+\delta+\varepsilon<2-2\delta
   \]
   正是 \(\varepsilon<1-3\delta\)。所以 \(Q(x,R)\cap A_u=\varnothing\)。

因此 \(\mu_c(Q(x,R))=D_\varepsilon\)，并且对归一化输入
\[
 g_\varepsilon=|E_\varepsilon|^{-1}\mathbf1_{E_\varepsilon},
 \qquad \|g_\varepsilon\|_{L^1(\mu_c)}=1,
\]
有 \(M_{c,\mu_c}g_\varepsilon(x)\ge D_\varepsilon^{-1}\) 对所有 \(x\in A\) 成立。

为保留严格水平集，固定任意 \(0<\tau<D_\varepsilon^{-1}\)，则
\[
 A\subset\{M_{c,\mu_c}g_\varepsilon>\tau\},
 \qquad \tau\mu_c\{M g_\varepsilon>\tau\}
 \ge\tau m(1+c)V.
\]
先令 \(\tau\uparrow D_\varepsilon^{-1}\)，再令 \(\varepsilon\downarrow0\)，即得
\[
 W_n(c)\ge\frac{m(1+c)(2\delta)^n}
 {(2(1+\delta))^n+c(2\delta)^n}
 =\frac{m(1+c)}{c+b_{n,\delta}}.
\]
所有极限都用于下界的上确界，不声称某个正 \(\varepsilon\) 在等号阈值达成。将分式从 \(m\) 中相减得到 (3.2)，随后 \(\delta\uparrow1/3\) 给出 (3.3)。证毕。

**输入已经是指标函数的常数倍。**将 \(g_\varepsilon\) 换成 \(\mathbf1_{E_\varepsilon}\)，同时把阈值乘以 \(|E_\varepsilon|\)，弱型商不变；所以即使把 (2.2) 的输入限为指标函数，也不能获得任何 \(m-a_nc^{-p}\)、\(p<1\) 的统一上界。

### 推论 3.2：排除所有 \(p<1\) 的密度比亏损

若有固定 \(a_n>0\)、\(p<1\) 使大 \(c\) 时
\(W_n(c)\le m-a_nc^{-p}\)，结合 (3.2) 得
\(a_nc^{1-p}\le m(b_{n,\delta}-1)\)，当 \(c\to\infty\) 矛盾。

### 反例为什么没有否证真实二阶上界

反例的原始 Lebesgue 高密度条件完全失效。若 \(x\in A\)，任何与 \(E_\varepsilon\) 相交的中心方体半径至少为 \(1-\delta-\varepsilon\)。故
\[
 M_c\mathbf1_{E_\varepsilon}(x)
 \le\left(\frac{\varepsilon}{1-\delta-\varepsilon}\right)^n
 \longrightarrow0.
\tag{3.4}
\]
它只能否证**丢掉 (1.1) 后**的无约束权范数加强。它不支持“真实 \(R_n(t)\to0\)”的结论。

## 4. 精确障碍：最优无约束权范数也只能恢复 \(U_n\)

### 命题 4.1

对每个 \(n\ge1\) 和 \(0<t\le4^{-n}\)，
\[
 \boxed{\quad
 \inf_{c>0}W_n(c)\frac{1+ct}{(1+c)(1-t)}
 =\frac{mt}{1-t}.
 \quad}\tag{4.1}
\]

证明：由 (3.3)，每个 \(c>0\) 的左侧被积表达式至少为
\[
 \frac{m(1+ct)}{(c+4^n)(1-t)}.
\]
它与 \(mt/(1-t)\) 的差恰为
\[
 \frac{m(1-4^nt)}{(c+4^n)(1-t)}\ge0.\tag{4.2}
\]
另一方面，Aldaz 所有测度弱型上界给 \(W_n(c)\le m\)，而 (3.3) 给 \(\lim_{c\to\infty}W_n(c)=m\)。令 \(c\to\infty\)，(2.3) 的右侧趋向 \(mt/(1-t)\)，完成证明。

这项障碍的范围很明确：辅助测度是 \((1+c\mathbf1_F)dx\)，使用的弱型信息只有“权在 \([1,1+c]\) 内”，并使用 (2.1) 的阈值。它不排除使用高密度见证、输入与权的耦合、或几何局部化的更精细论证。

## 5. 零背景原测度路线的更强障碍

保留上面的 \(E_\varepsilon,A\)，但改用有限测度
\[
 d\nu_{c,\varepsilon}=(\mathbf1_{E_\varepsilon}+c\mathbf1_A)dx.
\]
同一方体有 \(\nu_{c,\varepsilon}(Q)=|E_\varepsilon|+cV\)。以同样的严格阈值极限计算，固定任意 \(c>0\) 时，弱型常数至少为
\[
 \frac{mcV}{|E_\varepsilon|+cV}\longrightarrow m.
\]
结合所有测度的上界，说明：**即使密度只取 \(0,1,c\) 三个值，而且输入限为指标函数，对每个固定 \(c>0\) 的此类测度取上确界，弱型常数仍恰为 \(m\)。**所以原 \(\mathbf1_E+c\mathbf1_F\) 辅助测度的“分片常密度”本身完全不足以降低弱型常数。
