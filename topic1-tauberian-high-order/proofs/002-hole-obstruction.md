# 多尺度避邻孔机制的限制

以下模型体积不等式可取 n≥1；所有关于真实二阶尺度 tu 的讨论统一限定 n≥2。
本文不提供无约束 Tauberian 常数的上界。

## 2. 有限周期、任意孔形与多尺度的整孔阻碍

### 2.1 精确定义

设 \(\Lambda\subset\mathbb R^n\) 为满秩周期晶格，\(P\) 为有界基本胞，\(|P|=V\)。给定有限个有界可测孔 \(A_1,\ldots,A_N\)，要求全体平移孔

\[
\{A_i+\lambda:1\le i\le N,\ \lambda\in\Lambda\}
\]

两两几乎不交。设 \(m_i=|A_i|>0\)，\(M=\sum_i m_i<V\)。孔可跨越基本胞边界；只要求每个孔本身有界，不要求它位于所选基本胞内。定义

\[
D=\bigcup_{i,\lambda}(A_i+\lambda),\quad E=\mathbb R^n\setminus D.
\]

给定可测的已填部分 \(B_i\subset A_i\)，\(b_i=|B_i|\)。对几乎每个 \(x\in B_i\)，假设存在 \(r_x>0\)，满足

\[
\begin{aligned}
&|Q(x,r_x)\cap(D\setminus A_i)|=0,\\
&|Q(x,r_x)\cap A_i|\ge\theta m_i,\\
&|Q(x,r_x)\cap D|<t(2r_x)^n,
\end{aligned}
\tag{2.1}
\]

其中 \(0<\theta\le1\) 固定。对其它周期副本使用平移后的同一约定。\(\theta=1\) 就是整孔捕获；\(b_i\) 可以为零，即允许完全未填的孔。

令 \(v=(t/\theta)^{1/n}\)。若 \(v<1/2\)，则

\[
\boxed{\quad
1+\frac{\sum_i b_i}{V-\sum_i m_i}\le L_n(t/\theta).
\quad}
\tag{2.2}
\]

特别地，整孔版本的右端就是 \(L_n(t)\)。

### 2.2 密度点上的严格分离

删除每个孔的非密度点，并删除 \(B_i\) 上不具备见证的零测集。此后不改变 \(m_i,b_i\)。由 (2.1)，

\[
r_x>R_i:=\frac{m_i^{1/n}}{2v}.
\tag{2.3}
\]

若 \(y\) 是任意其它孔的密度点，则 \(y\notin Q(x,r_x)^\circ\)：否则 \(y\) 附近一个小立方体包含在 \(Q(x,r_x)\) 内，密度点性质将给出正的交体积，违反第一条。因此

\[
\|x-y\|_\infty\ge r_x>R_i.
\tag{2.4}
\]

该推理只需要每个已填点各自存在一个有限半径见证，不要求见证半径一致、有正的统一裕量、连续或可测。

对两个不同周期孔的已填点 \(x\in B_i+\lambda\)、\(y\in B_j+\mu\)，交换二者角色得

\[
\|x-y\|_\infty>\max\{R_i,R_j\}.
\tag{2.5}
\]

严格不等号对每对点成立；取集合间距离的下确界后只能保留“\(\ge\)”，这里没有把逐点严格性错误传给下确界。

### 2.3 膨胀集合与周期边界

对 \(b_i>0\) 定义开集

\[
S_i=B_i+(-R_i/2,R_i/2)^n.
\]

全体 \(S_i+\lambda\) 两两不交。事实上，一个公共点会给出两个已填点间距离

\[
\|x-y\|_\infty<\frac{R_i+R_j}{2}\le\max\{R_i,R_j\},
\]

与 (2.5) 矛盾。这一步使用开立方体，避免接触边界的无关问题。

若 \(b_j=0\)，把整个未填孔 \(A_j\) 直接保留为占用体积。它与所有 \(S_i+\lambda\) 几乎不交：若其一个密度点落入某 \(S_i+\lambda\)，则其距离对应 \(B_i+\lambda\) 中一点小于 \(R_i/2<R_i\)，违反 (2.4)。未填孔之间本来就几乎不交。

即便 \(S_i\) 或 \(A_i\) 跨越基本胞边界，周期计数仍精确：对任意有界可测集合 \(S\)，

\[
\int_P\sum_{\lambda\in\Lambda}1_{S+\lambda}(z)\,dz=|S|.
\]

所以由上述不交性，

\[
\sum_{b_i>0}|S_i|+\sum_{b_i=0}m_i\le V.
\tag{2.6}
\]

这也说明为什么不能对空的 \(B_i\) 直接套一个正半径的 Brunn–Minkowski 下界；空孔已单独计费。

### 2.4 Brunn–Minkowski 与标量不等式

对每个 \(b_i>0\)，Brunn–Minkowski 不等式给

\[
|S_i|^{1/n}\ge b_i^{1/n}+R_i.
\]

令 \(p_i=b_i/m_i\in(0,1]\)、\(c=1/(2v)>1\)，则

\[
|S_i|\ge m_i(p_i^{1/n}+c)^n=:m_i a(p_i).
\tag{2.7}
\]

这里只用到标准 Brunn–Minkowski 定理
\(|X+Y|^{1/n}\ge|X|^{1/n}+|Y|^{1/n}\)，第二个集合取边长 \(R_i\) 的开立方体；所有参与集合有界可测且 \(b_i>0\)。

对 \(0<q\le1\)，函数

\[
F(q)=\frac{q^n}{(q+c)^n-1}
\]

递增，因为

\[
F'(q)=\frac{nq^{n-1}[c(q+c)^{n-1}-1]}{[(q+c)^n-1]^2}>0.
\]

置

\[
\kappa=F(1)=\frac1{(1+c)^n-1}=L_n(t/\theta)-1.
\]

因而 \(p\le\kappa[a(p)-1]\)。由 (2.6)–(2.7)，

\[
\begin{aligned}
\sum_i b_i
&\le\kappa\sum_{b_i>0}m_i[a(p_i)-1]\\
&\le\kappa\left(V-\sum_{b_i=0}m_i-\sum_{b_i>0}m_i\right)
=\kappa(V-M).
\end{aligned}
\]

证得 (2.2)。证明覆盖任意固定有限的每胞孔数；没有声称它已处理每胞可数无穷孔的所有可测性／求和量词。对有正体积的可数族，若 (2.1) 对所有孔成立且总孔体积有限，同一非负求和与 Tonelli 论证也可进行，但本轮构造与结论只需要有限每胞版本。

### 2.5 对二阶机制的直接约束

若 \(\theta(t)=1-\eta(t)\)、\(\eta(t)\to0\)，则

\[
L_n(t/\theta)=1+2^nt-n2^{n+1}tu+2^nt\eta(t)
+O_n(tu^2+t^2+t\eta(t)^2+tu\eta(t)).
\tag{2.8}
\]

因此 \(\eta=o(u)\) 时，(2.2) 的 \(tu\) 系数仍为 \(-n2^{n+1}\)。这不是对真实算子的断言：它只说明，在仍然避开其它孔的机制里，欲使该上界允许一个更大的二阶系数，捕获损失至少要达到 \(u\) 阶。

式 (2.8) 是**模型上界的展开**，不能作为真实下界证据。

## 3. 同形立方孔：全部半径的点态公式

本节排除最直接的“保留原孔阵列，只改用跨孔半径”路线。

设 \(n\ge2\)，\(0<\delta\le1/(4n)\)，

\[
D_\delta=\bigcup_{k\in\mathbb Z^n}(k+[-\delta,\delta]^n).
\]

对 \(x\in[-\delta,\delta]^n\)，令

\[
s=\|x\|_\infty,\qquad r_0=1-\delta-s.
\]

则

\[
\boxed{\quad
\inf_{r>0}\frac{|Q(x,r)\cap D_\delta|}{(2r)^n}
=\left(\frac\delta{1-\delta-\|x\|_\infty}\right)^n.
\quad}
\tag{3.1}
\]

右端在 \(r=r_0\) 取到（接触邻孔边界不增加体积）。

### 3.1 张量分解

令一维周期缺陷 \(D^{(1)}_\delta=\mathbb Z+[-\delta,\delta]\)，

\[
\ell(a,r)=|[a-r,a+r]\cap D^{(1)}_\delta|.
\]

因为 \(D_\delta=(D^{(1)}_\delta)^n\)，由 Fubini，缺陷比例为

\[
\prod_{j=1}^n\frac{\ell(x_j,r)}{2r}.
\tag{3.2}
\]

### 3.2 小半径与首接触前

若 \(r\le2\delta\)，则每个 \(x_j\in[-\delta,\delta]\) 都满足 \(\ell(x_j,r)\ge r\)，故 (3.2) 至少为 \(2^{-n}\)，大于或等于 \((\delta/r_0)^n\)。后者用 \(r_0\ge1-2\delta\ge2\delta\)。

若 \(2\delta\le r\le r_0\)，本孔完全包含且邻孔尚未进入，缺陷比例恰为 \((\delta/r)^n\ge(\delta/r_0)^n\)。

### 3.3 首接触至半径 1

若 \(r=r_0+z\in[r_0,1]\)，选取 \(|x_{j_*}|=s\) 的坐标。在该坐标上，最近邻孔从 \(r_0\) 开始进入；因

\[
0\le z\le1-r_0=\delta+s\le2\delta,
\]

它贡献至少长度 \(z\)。其它坐标至少贡献本孔长度 \(2\delta\)。于是 (3.2) 至少为

\[
\left(\frac\delta{r_0+z}\right)^n\left(1+\frac z{2\delta}\right).
\]

除以目标值 \((\delta/r_0)^n\)，得到

\[
G(z)=\frac{1+z/(2\delta)}{(1+z/r_0)^n}.
\]

在指定区间内

\[
\frac{G'(z)}{G(z)}=\frac1{2\delta+z}-\frac n{r_0+z}\ge0,
\]

因为

\[
r_0-2n\delta-(n-1)z
\ge n r_0-2n\delta-(n-1)
\ge1-4n\delta\ge0.
\]

故 \(G(z)\ge G(0)=1\)。

### 3.4 所有大半径

先证一维辅助界：当 \(\delta\le1/4\)、\(|a|\le\delta\)、\(r\ge1\) 时，

\[
\frac{\ell(a,r)}{2r}\ge\frac32\delta.
\tag{3.3}
\]

写 \(r=k+q\)，\(k\ge1\) 为整数，\(q\in[0,1)\)。周期性给

\[
\ell(a,k+q)=4k\delta+\ell(a,q).
\]

若 \(q\le2\delta\)，则 \(\ell(a,q)\ge q\)，函数
\((4k\delta+q)/(2(k+q))\) 在该区间递增，故比例至少为 \(2\delta\)。若 \(q\ge2\delta\)，本孔被完全包含，故比例至少

\[
\frac{(2k+1)\delta}{k+q}\ge\frac{(2k+1)\delta}{k+1}\ge\frac32\delta.
\]

由 \(n\ge2\) 与 \(\delta\le1/(4n)\le1/8\)，有 \(r_0\ge3/4\)，所以

\[
\frac32\delta\ge\frac\delta{r_0}.
\]

对各坐标使用 (3.3)，得到全部 \(r\ge1\) 的所需下界。结合前三段及 \(r=r_0\) 的等号，证得 (3.1)。

### 3.5 真正周期 halo 的体积与优化

令 \(E_\delta=\mathbb R^n\setminus D_\delta\)。由 (3.1)，本孔中属于严格水平集的点恰满足

\[
\|x\|_\infty<1-\delta-\delta/u.
\]

因此每胞的填孔体积为

\[
b(\delta)=\left[2\min\{\delta,(1-\delta-\delta/u)_+\}\right]^n,
\tag{3.4}
\]

周期 halo 比值为

\[
1+\frac{b(\delta)}{1-(2\delta)^n}.
\tag{3.5}
\]

这里是实际最大算子的全部半径公式，而不是一个预指定半径代理值。

当 \(u\le1/(4n-2)\)，\(\delta_*=u/(1+2u)\) 落在本节范围内。对 (3.5) 在 \(0<\delta\le1/(4n)\) 优化，最大值为 \(L_n(t)\)：\(\delta\le\delta_*\) 时全孔填入分支递增；在 \(\delta_*<\delta<u/(1+u)\) 上，令 \(K=1+1/u>3\)，部分填孔增益为

\[
g(\delta)=\frac{[2(1-K\delta)]^n}{1-(2\delta)^n}.
\]

其导数符号为 \(-K+2^n\delta^{n-1}<0\)，因为 \(\delta<1/2\) 给 \(2^n\delta^{n-1}\le2<K\)。更大的 \(\delta\) 有 \(b(\delta)=0\)。严格水平只可能移除孔边界等零测集，不改体积。

这证明：在薄的同形晶格立方孔族内，改变孔径、只填孔的一部分、采用跨邻孔半径、采用只截本孔的小半径，都不提高该构造。
