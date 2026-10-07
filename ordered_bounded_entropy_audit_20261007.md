# 有界连续熵与阈值跳占用：独立解析审计

2026-10-07，gated_radial_energy_audit。只读审计 root 的 `ordered_continuous_bounded_entropy_20261007.md`；不修改该稿，不重跑其数值。本收据审查固定 \(c>0\)、固定物理尺度、固定阈值 \(\lambda>0\) 下的原有序传播。结论：熵恒等式、一次来源预算、下穿越常数均成立。它没有证明原移动尺度/FIRST/geom 交通受控。

## 1. 真实非齐次生成元及来源类别

记 \(B_i\) 为原各轴的非负平移不变生成元，\(G_{d,i}=(I+dB_i)^{-1}\)。各 \(G_{d,i}\) 是保质量的对称正卷积算子。令

\[
F_{z,i}=e^{-z}I+(1-e^{-z})G_{c,i}
=\frac{I+ce^{-z}B_i}{I+cB_i},\qquad
v_z=\prod_iF_{z,i}f.
\]

逐轴的有界算子身份为

\[
F'_{z,i}=(G_{ce^{-z},i}-I)F_{z,i},\qquad
v'_z=J_zv_z,\quad J_z=\sum_i(G_{ce^{-z},i}-I).
\tag{1}
\]

这不是把固定 \(G_c\) 的半群生成元代入；真正跳核在收缩，各轴的时钟率为 1。有限 mask 展开直接给任意 \(f\in L^1\) 的 \(L^1\)-可微性，无须先假设 \(f\in L^\infty\) 或 \(L\log L\)。保质量与正性给 \(v_z\ge0\)、\(\|v_z\|_1=W\)。且 \(\|J_z\|_{1\to1}\le2n\)。有限时间中 \(ce^{-z}>0\)，不需要在参数零端点微分。

## 2. 熵链式微分与对称 Fubini

对 \(t\ge0\)，

\[
\Phi_\lambda(t)=t-\lambda\log(1+t/\lambda),\quad
\Phi'_\lambda(t)=\frac{t}{\lambda+t},\quad
\Phi''_\lambda(t)=\frac{\lambda}{(\lambda+t)^2}.
\]

因此 \(0\le\Phi_\lambda(t)\le t\)，且 \(E_\lambda(g)=\int\Phi_\lambda(g)\) 在非负 \(L^1\) 锥上为 1-Lipschitz。

需准确区分：这并不宣称 \(E_\lambda\) 是整个 \(L^1\) 上一般 Fréchet-\(C^1\) 泛函。沿本题 \(L^1\)-\(C^1\) 正轨道的链式微分仍然合法。设 \(\Delta_h=v_{z+h}-v_z\)，用

\[
\Phi(v_z+\Delta_h)-\Phi(v_z)
=\Delta_h\int_0^1\Phi'(v_z+s\Delta_h)\,ds.
\]

差商 \(\Delta_h/h\to v'_z\) 于 \(L^1\)，斜率有界于 1 且依测度趋于 \(\Phi'(v_z)\)。先用差商的 \(L^1\) 误差，再对固定的可积权 \(|v'_z|\) 用有界收敛于测度，便得

\[
\frac{d}{dz}E_\lambda(v_z)=\int\Phi'_\lambda(v_z)J_zv_z.
\tag{2}
\]

每轴对称化所需的绝对 Fubini 有直接界

\[
\int dx\int G_{d,i}(dh)
\big|\Phi'_\lambda(v(x))[v(x+h)-v(x)]\big|\le2W.
\]

使用 \((x,h)\mapsto(x+h,-h)\) 与核的对称性，(2) 的第 \(i\) 项等于

\[
-\frac12\int dx\int G_{d,i}(dh)
[\Phi'_\lambda(v(x+h))-\Phi'_\lambda(v(x))]
[v(x+h)-v(x)].
\]

有理差身份精确为

\[
[\Phi'_\lambda(b)-\Phi'_\lambda(a)](b-a)
=\frac{\lambda(b-a)^2}{(\lambda+a)(\lambda+b)}\ge0.
\tag{3}
\]

故任意有限 \(Z\ge0\)，

\[
E_\lambda(f)-E_\lambda(v_Z)
=\frac\lambda2\sum_i\int_0^Zdz\int dx\int G_{ce^{-z},i}(dh)
\frac{[v_z(x+h)-v_z(x)]^2}
{[\lambda+v_z(x+h)][\lambda+v_z(x)]}
\le W.
\tag{4}
\]

所有熵都是有限量；没有 \(\infty-\infty\)。有限时间的积分也可先由绝对界控制，再用正 Tonelli。

若采用 bounded-source 逼近 \(f_m\uparrow f\)，不可声称 (3) 的能量密度随 \(m\) 单调：两个端点都增长并不使该比值单调。Fatou 单独足以给上界，但不足以给精确等式。直接的 (2) 避开此问题；或用 \(J_zv_{z,m}\to J_zv_z\) 于 \(L^1\) 和有界斜率，将有限时间的有符号导数积分传极限。

## 3. 无限未来与时钟转换

令 \(v_\infty=G_c^{\otimes n}f\)。乘积望远镜及各因子的 \(L^1\) 压缩性给

\[
\|v_z-v_\infty\|_1\le 2ne^{-z}W.
\tag{5}
\]

由熵的 Lipschitz 性，\(E(v_Z)\to E(v_\infty)\)。因此 (4) 的全未来积分精确等于 \(E(f)-E(v_\infty)\le W\)。此时可以对递增的时间区间用正 Tonelli；不需在 \(z=\infty\) 声称 \(L^1\) 生成元的参数零端点导数。

在 \(r=1-e^{-z}\) 坐标中，测度转换是

\[
dz=\frac{dr}{1-r},\qquad G_{ce^{-z},i}=G_{c(1-r),i}.
\tag{6}
\]

不可漏掉 \((1-r)^{-1}\)，也没有额外乘 \(c\)。这些结论逐固定物理尺度成立；不提供尺度上确界的免费来源预算。

## 4. 下穿越与真正密度加权跳流

给定 \(\delta>0\)，记 \(D_\delta=\{v_z(x)\ge(1+\delta)\lambda,\ v_z(x+h)\le\lambda\}\)，并定义

\[
I_{\lambda,\delta}=\sum_i\int_0^\infty dz\int dx\int G_{ce^{-z},i}(dh)\,1_{D_\delta}.
\]

在 \(U=v_z(x)/\lambda\ge1+\delta\)、\(0\le V=v_z(x+h)/\lambda\le1\) 上，

\[
Q(U,V)=\frac{(U-V)^2}{(1+U)(1+V)}
\ge \frac{\delta^2}{2(2+\delta)}.
\tag{7}
\]

比值随 \(U\) 增、随 \(V\) 减，其最小值在边界 \((1+\delta,1)\)。下穿越及其反向上穿越不交，且具有相同的积分；因此 (4) 的 \(1/2\) 被这两个方向精确抵消，得

\[
\lambda I_{\lambda,\delta}\le
\frac{2(2+\delta)}{\delta^2}W.
\tag{8}
\]

特别 \(\delta=1\) 给 \(\lambda I_{\lambda,1}\le6W\)。含等号的事件没有问题；阈值间仍有严格间隔。

这里 \(I\) 是 Lebesgue 空间边测度与真实跳核的占用，不是由来源密度加权的概率跳流。另有可审查的加权推论：在同一事件上，

\[
\frac{U}{Q(U,V)}
=\frac{U(1+U)(1+V)}{(U-V)^2}
\le \frac{2(1+\delta)(2+\delta)}{\delta^2}.
\tag{9}
\]

先取 \(V=1\)，再用 \(U(U+1)/(U-1)^2\) 在 \(U>1\) 递减即可。结合对称的两个方向，

\[
\sum_i\int_0^\infty dz\int dx\int G_{ce^{-z},i}(dh)
v_z(x)1_{D_\delta}
\le \frac{2(1+\delta)(2+\delta)}{\delta^2}W.
\tag{10}
\]

\(\delta=1\) 的常数为 12。式 (10) 才是本固定传播中、来源密度 \(v_z(x)\) 的向下真实跳流付款。它仍不是原 geom 的 conditional terminal/FIRST 交通恒等式。

## 5. 门的有效范围及最小未付接口

任何在 (4) 的正能量载体上的可测 \(0\le a(z,x,h,i)\le1\) 都使积分下降；无需额外可预测性条件，因为此处是一个确定的正积分。随机门若由归一化条件概率平均得到，也适用。非对称的下穿越门同样可先以 1 支配，再使用 (8) 或 (10)，所以不额外损失常数 2。

但这并不准许把任意来源标签、LCA 或细 history 重复装饰在能量上。若需要保留这些标签，必须给出一个正 lift，使其边测度边际由上述单份能量/跳流支配；否则每个标签都复制能量会重新收取来源质量。亦不能将固定标量 \(\lambda\) 任意换为接收点或来源相关阈值而沿用 (3) 的对称链式证明。

原 actual 交通可能在两端密度相等时仍非零，而 (4) 的能量为零；因此从“任意门可乘正能量”到“原门内交通可支付”仍缺少真实的通量支配或正分解身份。本收据通过原固定有序传播的通用 \(W\) 预算，未闭合原移动 hard 平均与全部 FIRST/history 的空间费用。

## 6. 数值范围

未运行新数值。root 专属标量 Fraction 与 cosine 连续流检查不由本审计重复；此收据的成立性来自 (1)–(10) 的解析推导。离散 mask Jensen 预算亦是另一独立组件，未用于本恒等式。
