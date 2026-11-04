# 二维半阶 pair 差体与真实体积上界复审

## 结论

对 \(\alpha>1/\sqrt2\)，pair 差体方法严格给出全体输入上界
\[
 C_{1,\mathcal Q_c}(\alpha)
 \le U(\alpha):=1-\alpha^2+(2\alpha^2-1)
 \log\frac{2\alpha^2-1}{\alpha^2}.
\]
该上界只利用任意两点可由同一输入共同严格见证的必要条件，故不声称可达到，也不闭合真实三阶系数。差体、凸化、积分常数和端点展开均正确。

## 1. 两见证到差体

令两严格见证边长为 \(\ell,s\in(\alpha,\alpha^{-1})\)，\(z=(\ell+s)/2\)，两中心坐标差绝对值为 \((u,v)\)。第 17 轮宽度定理给 \(u,v<W=\alpha^{-1}-\alpha<\alpha\)。交叠面积满足
\[
 I\le(z-u)(z-v),
\]
而严格 N=2 容量条件为 \(2\alpha z<1+I\)。因此
\[
 G(z)=z^2-(u+v+2\alpha)z+1+uv>0.
\]
\(G\) 凸且 \(z\in(\alpha,\alpha^{-1})\)，故闭区间某一端点为正。写 \(S=u+v\)，有
\[
 G(\alpha)=\alpha(W-S)+uv,\quad
 G(\alpha^{-1})=(W-S)/\alpha+uv,
\quad G(\alpha)-G(\alpha^{-1})=W(S-W).
\]
若 \(S<W\)，第一式本身为正；若 \(S\ge W\)，第二式若正则第一式不小于它。故必有 \(G(\alpha)>0\)，即
\[
 (\alpha-u)(\alpha-v)>2\alpha^2-1.
\]
反向取两个相同边长 \(z\downarrow\alpha\) 的方体：连续性使 \(G(z)>0\)，个体容量严格可行，N=2 判据便实现任意该不等式内的差向量。于是
\[
 H-H\subset K_\alpha,
\quad K_\alpha=\{d:|d_j|<\alpha,\ (\alpha-|d_1|)(\alpha-|d_2|)>2\alpha^2-1\}
\]
是精确的 pair 可行体；“精确”仅指 pair 层面。

## 2. 凸化和体积

\(-\log(\alpha-|x|)\) 在 \((-\alpha,\alpha)\) 凸，所以 \(K_\alpha\) 为开、对称、凸体。宽度界也使非空光晕有界开；因此
\[
 \operatorname{conv}H-\operatorname{conv}H
 =\operatorname{conv}(H-H)\subset K_\alpha.
\]
对 \(\operatorname{conv}H\) 与其相反集应用 Brunn--Minkowski，得
\(4|H|\le4|\operatorname{conv}H|\le|K_\alpha|\)。正象限面积为
\[
 \int_0^{W}\left(\alpha-
 \frac{2\alpha^2-1}{\alpha-u}\right)du
 =1-\alpha^2+(2\alpha^2-1)
 \log\frac{2\alpha^2-1}{\alpha^2},
\]
且 \(|K_\alpha|/4\) 正是该数。这证明上式 \(U\)。\(K_\alpha/2\) 满足
\((K_\alpha/2)-(K_\alpha/2)=K_\alpha\)，故 BM/pair-only 体积常数本身不能改善；但它不保证一个共同输入同时实现所有差对。

## 3. 边界与展开

严格水平的右连续关系
\(C(\alpha)=\sup_m C(\alpha+1/m)\)（只取 \(\alpha+1/m<1\)）把上界延至
\(\alpha=1/\sqrt2\)，此时 \(U=1/2\)。CAS 展开为
\[
 U(1-t)=2t^2+\frac{10}{3}t^3+\frac{43}{6}t^4+O(t^5).
\]
结合矩形下界
\(R(1-t)=2t^2+2t^3+O(t^4)\)，只可得
\[
 2\le\liminf\frac{C(1-t)-2t^2}{t^3}
 \le\limsup\frac{C(1-t)-2t^2}{t^3}\le\frac{10}{3}.
\]
这不预设三阶极限存在。

## 4. 上包络尾项与高矩应用

置 \(x=1-\alpha^2\in(0,1/2)\)，则
\[
 U=x+(1-2x)\log\frac{1-2x}{1-x}
 =\sum_{j\ge2}\frac{2^j-j-1}{j(j-1)}x^j.
\]
系数均正，故余项至少为下一项；并以
\(b_j\le2^j/[j(j-1)]\) 得
\[
 \sum_{j\ge N+1}b_jx^j
 \le\frac{(2x)^{N+1}}{N(N+1)(1-2x)}.
\]
这是 \(U\) 的任意阶 sharp 截断包络，不能误称为 \(C\) 的展开。

若 \(\mathfrak F(p)=\sup_{|E|=1}\int(M_{1,\mathcal Q_c}\mathbf1_E)^p\)，层析和上界给
\[
 \mathfrak F(p)\le2pB(p,3)+\frac{10}{3}pB(p,4)+O(p^{-4})
 =4p^{-2}+8p^{-3}+O(p^{-4}),
\]
远离顶层部分由既有弱型估计指数小。固定 \(Q_0\) 的
\(2t^2+2t^3+O(t^4)\) 给
\(\mathfrak F(p)\ge4p^{-2}+O(p^{-4})\)：两项的 \(p^{-3}\) 系数抵消。因而
\[
 0\le\liminf p^3(\mathfrak F(p)-4p^{-2})
 \le\limsup p^3(\mathfrak F(p)-4p^{-2})\le8.
\]
上界中只将上确界置于层析积分外；下界固定 \(Q_0\)，未将矩形上确界 \(R\) 不当地积分化。

## 正式稿 06j 与 `mom:pair-bracket` 最终复核

定稿通过复审，无新增缺口。

- `pb:body` 保留了严格尺度窗、严格 N=2 阈值和“反向输入可依差向量变化”的限定；\(K_\alpha\) 的凸性、开性及差体包含均正确。
- `pb:right-continuity` 对每个输入先以严格水平集并集使用测度下连续性，再交换两个上确界，逻辑正确；因此可从右把 \(\mathcal U\) 上界延至 \(\alpha=1/\sqrt2\)，连续值为 \(1/2\)。
- `pb:volume` 中 \(\operatorname{conv}H-\operatorname{conv}H=\operatorname{conv}(H-H)\subset K_\alpha\) 对开凸 \(K_\alpha\) 合法，Brunn--Minkowski 给出的 \(4\) 因子及正象限对数积分均正确。三阶仅给 \([2,10/3]\) 的 liminf/limsup 括号，未越界宣称极限。
- `pb:pair-only` 精确区分“差对逐一可实现”与“同一输入实现整个差体”，范围正确。
- `pb:all-orders` 的正级数、首项尾下界及几何级数尾上界正确，且清楚标注为上包络而非 \(C\) 的展开。
- `mom:pair-bracket` 的上界只把统一 \(\mathcal U\) 放在层析积分外，低层弱型项指数小；下界固定 \(Q_0\)，\(2pB(p,3)\) 的 \(-12p^{-3}\) 与 \(2pB(p,4)\) 的 \(+12p^{-3}\) 正确抵消。因此 \([0,8]\) 的高矩次阶括号成立。
