# 第21轮独立审查：低阶背景与分离奇点

## 结论

候选的近端点精确公式成立。宜把背景条件写成：\(b>0\) a.e.、\(b\in L^1_{\rm loc}\)，且在每个主奇点 \(z_j\) 满足
\[
\varepsilon^{n-\tau_j}b(z_j+\varepsilon\,\cdot)\longrightarrow0
\quad\hbox{in }L^1_{\rm loc}.\tag{R21.0}
\]
局部本质有界是该条件的充分条件。

设 \(A_t=\max_jC^{\rm out}_{K,h_j}(1-t)\) 和
\(B_t=C^{\rm out}_{K,b}(1-t)\)。若 \(E\) 是 \(b+v\) 的合格输入，正性和有限总分母保证 \(b(E),h_j(E)\) 都正且有限。对严格光晕 \(H\) 逐项加和，得到
\[
 \int_H(b+v)\le\max\{B_t,A_t\}\int_E(b+v).\tag{R21.1}
\]
另一方面，(R21.0)、异点分离以及有界输入和紧严格光晕的缩放给
\(C^{\rm out}_{K,b+v}(1-t)\ge C^{\rm out}_{K,h_j}(1-t)\)。
由第20轮，\(A_t-1\sim S_*t^{\theta_*}\)，其中
\(\theta_*=\min_j\tau_j/n\) 且 \(S_*>0\)。若
\(B_t-1=o(t^{\theta_*})\)，则充分小 \(t\) 后 \(B_t<A_t\)，故
\[
C^{\rm out}_{K,b+v}(1-t)=A_t.\tag{R21.2}
\]
这也给出一阶、已适用二阶、以及共同中心多指数近 \(1\) 的不变性。

若 \(0<a\le b\le B<\infty\)，无权严格高水平溢出界给
\[
C^{\rm out}_{K,b}(1-t)-1\le (B/a)B_Kt/(1-t)=o(t^{\theta_*}),
\]
所以这是直接可用的充分条件。

## 弱型范数

同样的逐项加和和局部缩放给
\[
\max_jW^{\rm out}_{K,h_j}\le W^{\rm out}_{K,b+v}
\le\max\{W^{\rm out}_{K,b},\max_jW^{\rm out}_{K,h_j}\}.\tag{R21.3}
\]
一般不等号可严格。背景“有界时有限”必须指全局上下有界
\(0<a\le b\le B<\infty\)；若仅上有界而允许趋零，不能据此断言弱型范数有限。

## 全阈值反例与算术核验

取 \(n=1\)、\(K=[-1,1]\)、\(v=|x|^{-1/2}\)、
\[
b=1+64\mathbf1_{(1,2)},\qquad E=(0,1),\qquad\alpha=1/4.
\]
直接优化中心区间的交长度给严格光晕 \((-1,2)\)（模端点）：在左边
\(x=-a\) 的最大平均为 \(1/[2(1+a)]\)，在右边 \(x>1\) 的最大平均为 \(1/(2x)\)。故
\[
\frac{\mu_{b+v}((-1,2))}{\mu_{b+v}((0,1))}
=\frac{69+2\sqrt2}{3}.
\]
纯幂权有 \(W^{\rm out}_{K,v}\le2(1+\sqrt2)\)，因而
\(C^{\rm out}_{K,v}(1/4)\le8(1+\sqrt2)\)。前者严格更大，等价于
\(45>22\sqrt2\)，其平方为 \(2025>968\)。所以近端点公式不能延伸为所有阈值的公式。

同一输入还给
\[
W^{\rm out}_{K,b+v}\ge(69+2\sqrt2)/12>2(1+\sqrt2)\ge W^{\rm out}_{K,v},
\]
差的分子为 \(45-22\sqrt2\)，平方差为 \(1057\)。等指数 \(m=2\) 的临界因子为 \(1/2\)，故临界留数也严格不同。

## 显式端点区间

纯幂权的 (H.9)（\(\tau=1/2,D_0=2\)）给
\[
C^{\rm out}_{K,v}(1-t)-1\ge\frac{\sqrt t}{\sqrt{1+2t}-\sqrt t}\ge\frac{\sqrt t}{1+t}.
\]
该背景有 \(a=1,B=65,B_K=2\)，故背景增量至多 \(130t/(1-t)\)。对
\(0<t\le2^{-16}\)，有 \(130\sqrt t(1+t)<1-t\)，并且下界构造的
\(3D_0h<R\) 自动成立。因此该明确区间已有
\[
C^{\rm out}_{K,b+v}(1-t)=C^{\rm out}_{K,v}(1-t),
\]
但弱型范数和留数仍不同。

## 可重合奇点的安全合并推论

设有限和的 \(z_j\) 可重合，另有全局上下有界背景 \(0<a\le b_0\le B\)。令
\(\tau_*=\min_j\tau_j\)，并在每个位置合并最小次数项
\[
H_z=\sum_{z_j=z,\ \tau_j=\tau_*}c_jh_j.
\]
每个 \(H_z\) 仍满足 (H.1)，次数为 \(\tau_*-n\)，且不同位置已分离。将所有高次数项和 \(b_0\) 归入背景。逐项上界给
\[
C^{\rm out}_{K,b}(1-t)-1\le\max\{O(t),\max_{\tau_j>\tau_*}O(t^{\tau_j/n})\}=o(t^{\tau_*/n}).
\]
在主点，同位置高次数项缩放为 \(\varepsilon^{\tau_j-\tau_*}h_j\)，异点项局部有界，\(b_0\) 有界；故 (R21.0) 成立。于是充分小 \(t\) 时
\[
C^{\rm out}_{K,w}(1-t)=\max_zC^{\rm out}_{K,H_z}(1-t).
\]
首项必须取 \(\max_zS(H_z)\)，而非合并前单分量的最大值；这解释第20轮碰撞例。二阶项要求所有首项并列的 \(H_z\) 满足相应的正则性和 \(q\) 窗口。全部阈值的精确最大公式仍仅应在原始奇点两两分离、无背景时宣称。

## 仓库 BG.1--BG.12 终审

已实际阅读 paper/background_extension.tex 及 paper.tex 的输入位置。BG.1--BG.10
与前述核验一致，标签均在被引用前定义；background_extension 在 mixture_extension
之后输入，故没有循环或前置引用。BG.11--BG.12 的紧支撑扰动也成立：
\[
\int|w_\varepsilon-u|=64\sqrt\varepsilon,
\]
而 \(E=(0,\varepsilon)\) 的严格光晕为
\((-\varepsilon,2\varepsilon)\) 模端点，并且令 \(q=\sqrt\varepsilon\le1\) 后，
\[
3(3q+66+2\sqrt2)-(q+2)(69+2\sqrt2)
=(60+2\sqrt2)(1-q)\ge0.
\]
所以 (BG.12) 的下界正确。MX.2 对常数背景还给
\(W^{\rm out}_{c,u}=W^{\rm out}_{c,v}\)，从而该段确实证明弱型范数与固定
阈值常数不连续；建议在 BG.12 后补出这一个引用句，使这一点无需读者回溯。

对父代理所问的 \(b_0\equiv0\)：可安全放宽为
\[
b_0\equiv0\quad\hbox{或}\quad0<a\le b_0\le B.
\]
若存在高次数项，则它们的和正 a.e.，仍是满足 BG.1--BG.2 的背景；若没有
高次数项，则所有项都已被合并为互异中心的 \(G_z\)，直接应用 MX.2（事实上
得到全部阈值的等式）。建议在 grouped-singularities 推论后以一句分支说明，
不要在证明中把零背景代入需要 \(b>0\) a.e. 的 BG 定理。


# 第 21 轮独立审查：同点角向合并的 \(S\) 单调性

## 结论

候选引理成立，且等号刻画正确。应在正式证明中明确 cone measure
\(\sigma_K\) 已归一化为概率测度，以及所有被相加的角向密度非零且为非负；
本文的 (H.1) 正性假设满足这些要求。

令 \(h(r\omega)=r^{\tau-n}a(\omega)\)，其中
\(\omega\in\partial K\)，并取
\[
p=\frac{n}{n-\tau}>1.
\]
锥坐标的 Jacobian 可写成
\[
dx=n|K|\,r^{n-1}dr\,d\sigma_K(\omega).
\]
因此
\[
\kappa(h)=\frac{n|K|}{\tau}\,\mathbb E_\sigma a,\qquad
D(h)=|K|\,\mathbb E_\sigma a^p.
\]
代回 (H.2) 立即得到精确而无额外常数的式子
\[
S(h)=\frac{\|a\|_{L^p(\sigma_K)}}{\mathbb E_\sigma a},
\qquad
\kappa(h)S(h)=\frac{n|K|}{\tau}\|a\|_{L^p(\sigma_K)}.\tag{R21A.1}
\]

对于同次数的有限正和 \(H=\sum_i c_ih_i\)，令
\(A_i=c_ia_i\)。Minkowski 不等式给
\[
\kappa(H)S(H)
=\frac{n|K|}{\tau}\left\|\sum_iA_i\right\|_p
\le\sum_i\frac{n|K|}{\tau}\|A_i\|_p
=\sum_i\kappa(c_ih_i)S(h_i).\tag{R21A.2}
\]
注意 \(S(c_ih_i)=S(h_i)\)。再除以正的
\(\kappa(H)=\sum_i\kappa(c_ih_i)\)，得到
\[
S(H)\le\frac{\sum_i\kappa(c_ih_i)S(h_i)}
{\sum_i\kappa(c_ih_i)}
\le\max_iS(h_i).\tag{R21A.3}
\]

等号分析有两层。第二个不等式取等要求所有正权重分量都有同一最大
\(S(h_i)\)；第一层 Minkowski 取等（\(1<p<\infty\)）要求全部非零
\(A_i\) 在 \(L^p\) 中为同一非负函数的正标量倍。因此
\[
S\!\left(\sum_ic_ih_i\right)=\max_iS(h_i)
\quad\Longleftrightarrow\quad
h_i\text{ 两两 a.e. 成正比例}.
\tag{R21A.4}
\]
单个分量时此条件按惯例是空真。正性是必要的；若允许符号改变，Minkowski
的严格等号叙述和权模型都不适用。

所以 BG.7 中同位置、同最小次数项必须先合并：合并不会提高端点首项系数；
若角向形状不成比例，则严格降低。这是第 20 轮碰撞例的结构性解释。

## 方体算例

在方体 cone probability 下，\(\mathbb E\theta_1=0\)，且
\[
\mathbb E\theta_1^2=\frac{n+2}{3n}.
\]
可由 \(X=R\Theta\) 的锥分解复核：均匀 \(X\in[-1,1]^n\) 满足
\(\mathbb EX_1^2=1/3\)、\(\mathbb ER^2=n/(n+2)\)。
取 \(p=2\)（即 \(\tau=n/2\)）及
\[
a_\pm=1\pm\delta\theta_1,\qquad0<\delta<1,
\]
则
\[
S(a_\pm)=\sqrt{1+\delta^2\frac{n+2}{3n}}>1,
\qquad
S(a_++a_-)=S(2)=1.
\]
本机 Python 对 \(n=1,2,3,8\) 的有理矩
\(1,2/3,5/9,5/12\) 作了独立算术复核。该例严格展示非比例同点合并导致
首项下降。

建议在正文中表述为“同次数的正分量”；不同次数的同点项不能直接用本引理，
而应按 BG.7 先把最小次数项合并、较高次数项转入切向消失背景。

主代理裁决：背景恒零的空余项分支直接引用MX.2；已补W(u)=W(v)及BG.12的精确差因子。角向命题限于已满足H.1的同次权，未将一般Lp角向函数的代数恒等式冒充已证端点类；补锥环带积分验证，并修正角向常数2与空间权2ρ^(τ−n)的区别。
