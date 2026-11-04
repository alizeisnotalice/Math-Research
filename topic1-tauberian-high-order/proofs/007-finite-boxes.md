# 任意有限轴平行盒截断的等孔铺砌上界

固定n≥2、u=t^(1/n)≤1/(128n)、0<δ≤1/(64n)。沿用 `005-tiling-holes.md` 的任意单位轴平行立方铺砌及中心等孔。对任意R>0，若E_R具有正测度，则其真实halo比值≤L_n(t)。以下证明不要求R趋于无穷；所用全半径点态公式和b(δ)优化见上述已审附件。

## 4. 新增：任意有限盒截断的统一上界

令

\[
B_R=[-R,R]^n,\quad E_R=E_\delta\cap B_R,
\quad m=(2\delta)^n,\quad
w=\delta/u,\quad a=(1-2t-u)w.
\]

设 `u<=1/(128n)`，因此 `t<2^(-n)` 且 `1-2t-u>0`。

### 4.1 实际捕获点的边界距离

有限 halo 位于 `B_R` 内。令 `x in H_t(E_R)\E_R`，则 `x` 属于某个孔，且有实际严格见证半径 `r`。

当 `r<=2delta`，仅本孔就占见证至少 `2^(-n)>t`，不可能合法。当 `r>2delta`，本孔完整包含，所以

\[
r>\delta/u=w.
\tag{4.1}
\]

记 `d=dist_infty(x,partial B_R)`，即到最近盒面的距离。若 `r>d`，盒外一个薄片占见证比例至少 `(r-d)/(2r)`；严格见证因此迫使

\[
d>(1-2t)r>(1-2t)w.
\tag{4.2}
\]

若 `r<=d`，同一不等式也显然成立。孔中心 `lambda` 到每个盒面的距离从而大于

\[
(1-2t)w-\delta=a.
\tag{4.3}
\]

存在捕获点还迫使 `R>(1-2t)w`。

### 4.2 a>=1/2：完整单位砖分配

由 (4.3)，每个有捕获点的孔，其完整单位砖均包含于 `B_R`。这些单位砖互不交，每块在 `E_R` 内的体积恰为 `1-m`。以 `N_cap` 表示有捕获点的孔数，原稿精确点态公式又给每孔捕获量不超过 `b(delta)`。故

\[
|E_R|\ge N_{\rm cap}(1-m),\qquad
|H_t(E_R)\setminus E_R|\le N_{\rm cap}b(\delta).
\]

当 `N_cap>0`，得

\[
\frac{|H_t(E_R)|}{|E_R|}
\le1+\frac{b(\delta)}{1-m}\le L_n(t).
\]

当 `N_cap=0`，比值为 `1`（忽略密度定理的零测例外），同样成立。

### 4.3 a<1/2：严格低于首项的粗界

此时

\[
w<\frac1{2(1-2t-u)}<1,\qquad \delta<u.
\]

有捕获点的孔中心落入 `B_(R-a)`。其不交单位砖落入 `B_(R-a+1/2)`，所以

\[
N_{\rm cap}\le[2R+1-2a]^n.
\tag{4.4}
\]

所有与 `B_R` 相交的孔中心落入 `B_(R+delta)`，对应单位砖在 `B_(R+delta+1/2)`，故

\[
|E_R|\ge(2R)^n-m[2R+2\delta+1]^n.
\tag{4.5}
\]

以整孔体积粗估捕获量，再用 `R>(1-2t)w` 和 `1-2a>0`，有

\[
\begin{aligned}
\frac{m[2R+1-2a]^n}{(2R)^n}
&=t\left[2w\left(1+\frac{1-2a}{2R}\right)\right]^n\\
&\le t\left[\frac{1+2wu}{1-2t}\right]^n.
\end{aligned}
\tag{4.6}
\]

同时

\[
\frac{m[2R+2\delta+1]^n}{(2R)^n}
\le t\left[2w+\frac{1+2\delta}{1-2t}\right]^n
\le(4u)^n.
\tag{4.7}
\]

最后一步使用 `w<1,delta<u,u<=1/256`，此时括号严格小于 `4`。

所以有限 halo 增益满足

\[
\frac{|H_t(E_R)\setminus E_R|}{|E_R|}
\le t\frac{[(1+2u)/(1-2t)]^n}{1-(4u)^n}.
\tag{4.8}
\]

给定 `u<=1/(128n)`，有

\[
(1+2u)^n\le\frac{64}{63},\quad
(1-2t)^n\ge\frac{16383}{16384},\quad
1-(4u)^n\ge\frac{4095}{4096}.
\]

说明：第二条由 `nt<=n u^2<=1/(16384n)<=1/32768` 及 Bernoulli 不等式；其余同样直接。于是 (4.8) 严格小于

\[
t\frac{64}{63}\frac{16384}{16383}\frac{4096}{4095}<2t.
\]

另一方面，最优孔宽 `delta_*=u/(1+2u)` 在允许范围内，且

\[
L_n(t)-1\ge\frac{2^nt}{(1+2u)^n}
\ge\frac{63}{64}2^nt\ge\frac{63}{16}t>2t.
\]

故小 `a` 分支也严格小于 `L_n(t)`，完成 (F)。这是实际有限输入的结论；不要求 `R` 趋于无穷或大于任何 `u` 的负幂。
