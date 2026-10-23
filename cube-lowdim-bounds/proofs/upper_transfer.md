# 参数证书到一般输入的解析传递

令U=R²/(R−1)−epsilon/6。本文引用Lin通用障碍分解定理，未声称重新Lean编译。所有标量前提须由每轮证书独立验证。

## 2. 核与生成元条件

写 \(r=|u|+|v|\)，在原点之外定义
\[
K(u,v)=\frac{(R/r-1)_+}{R-1}-\varepsilon(1-r)_+^2,
\]
原点赋值为 1（零测点，不改变积分）。核偶、非负、紧支撑且可积。在 \(0<r\le1\)，
\[
K-1=\frac{1-r}{r(R-1)}
\bigl[R-\varepsilon(R-1)r(1-r)\bigr]\ge0,
\]
因为 \(r(1-r)\le1/4\)，证书验证 \(R-\varepsilon(R-1)/4>0\)。在 \(1\le r\le R\) 凸起项为零，其他区域显然非负。菱形半径的积分公式给出
\[
\int_{\mathbb R^2}K=4\int_0^R rK(r)\,dr
=\frac{2R^2}{R-1}-\frac{\varepsilon}{3}=2U.
\]

采用作者的未归一化生成元
\[
Ag(z)=\sum_{j=1}^2\int_{\mathbb R}
\frac{g(z+te_j)+g(z-te_j)-2g(z)}{t^2}\,dt.
\]
令 \(\psi=(1-r)_+^2\)、\(H_R=(1-R/r)_+/(R-1)\)，在原点之外
\[
K=\frac{R}{R-1}\frac1r-\frac1{R-1}+H_R-\varepsilon\psi.
\]
作者 `BaseDensity.lean` 的角向下界对每个 \(R>1\) 都成立；`BumpGenerator.lean` 的上界不含 R。在两坐标轴之外，\(0<r<1\) 时它们给出
\[
AH_R\ge\frac4{R(R-1)}J(2r/R),\qquad A\psi\le4T(r),
\]
其中
\[
J(s)=\frac12-\frac2s-\frac4{s^2}\log(1-s/2),
\quad T(r)=2(2-r)(1-\log2+\log r)-2(1-r)\log(1-r).
\]
在 \(1\le r<R\)，分别有 \(AH_R\ge4/[R(R-1)]\) 与 \(A\psi\le4(2-2\log2)\)。因此足够验证
\[
Q(r):=J(2r/R)-\gamma T(r)\ge0\quad(0<r<1),
\quad1-\gamma(2-2\log2)\ge0,
\qquad\gamma=\varepsilon R(R-1).
\]
在 \(r\ge R\)，\(K(z)=0\) 而其余取值非负，生成元非负直接由积分定义得到。轴及分界线为零测集；但这一步尚不等于分布意义的正性，见第 4 节。

## 3. 全区间有理认证，不用有限采样替代证明

\(J(s)=1+\sum_{k\ge0}s^{k+1}/[2^{k+1}(k+3)]\) 在 \((0,2)\) 收敛，所以 \(J\ge1\) 且凸。
在 \((0,1/2]\)，\(T'>0\)，且 \(T(1/2)=3-5\log2<0\)，故 \(Q\ge1\)。
又
\[
T''(r)=-\frac{2(2-r)}{r^2(1-r)}\le-\frac{27}{2}.
\]
后一个估计由 \(r^2(1-r)\le4/27\) 与 \(2-r\ge1\) 得出。证书要求 \(\gamma>4/3\)，从而 \(Q''\ge18\)。任意 \(r_0\in(1/2,1)\) 均满足
\[
Q(r)\ge Q(r_0)+Q'(r_0)(r-r_0)+9(r-r_0)^2
\ge Q(r_0)-\frac{|Q'(r_0)|^2}{36}.
\]
取证书指定的有理 r0。对每个正有理数 q，置 \(z=(q-1)/(q+1)\)，使用
\[
\log q=2\sum_{k=0}^{99}\frac{z^{2k+1}}{2k+1}+E,
\qquad |E|\le\frac{2|z|^{201}}{201(1-z^2)}.
\]
所有运算使用 `Fraction`，随后向外取整到分母 \(10^{35}\)。用这些区间计算 Q、Q'，验证 Q(r0)−|Q′(r0)|²/36>0。

证书使用完全有理运算验证严格正余量、外区条件和核支配条件，见 `results/round1_upper.json`。

## 4. 分布表示与全尺度传递

这一节明确说明新参数为何满足通用传递定理，避免把几乎处处生成元正性直接当作最终结论。

作者 `Delta.lean` 给出 \(G=1/(2\pi r)\) 的基本解恒等式 \(AG=-2\pi\delta_0\)。令 \(g=H_R-\varepsilon\psi\)，则 g 有界、Lipschitz，且在坐标轴和 \(r=R\) 之外分段二阶光滑；\(r=1\) 的凸起导数连续，其二阶导数有界。对紧支撑测试函数，步长积分在折线附近至多出现 \(C(1+|\log\delta|)\) 型局部奇异性，\(\delta\) 为到这些折线的距离，故可积。由二阶差分的平移换元与 Fubini，有 \(\int gA\varphi=\int (Ag)\varphi\)。这正是 `PairingBound.lean` 的证明；其常数只需 R>1、有限 epsilon，并不依赖原来的有理取值。

因此新核在分布意义下满足
\[
AK=\beta-a\delta_0,\qquad
\beta=(Ag)\,dz\ge0,\qquad a=\frac{4\pi^2R}{R-1}>0.
\]
\(Ag=AK\ge0\) 几乎处处由第 2、3 节保证。其局部可积性已说明。用 \(0\le\chi_L\le1\)、\(\chi_L\to1\)、\(\|A\chi_L\|_\infty\le C/L\) 的光滑截断测试分布恒等式，得到
\(\int\chi_L\,d\beta\le a+C\|K\|_1/L\)。Fatou 推出 \(\beta\) 为有限正测度。此处没有遗漏轴上的额外奇异测度，因为对有界部分的绝对可积差分换元已经给出了完整分布表示。

现调用 Lin 的通用障碍分解定理：对有界非负紧支撑 f 及 cap \(\kappa>0\)，存在非负可积 u、\(0\le\sigma\le\kappa\)，使
\[
Au=\sigma-f,\qquad |\{u>0\}|\le\|f\|_1/\kappa.
\]
对每个 \(K_s(z)=s^{-2}K(z/s)\)，尺度齐次性将同一分布表示转为有限正测度减非负点质量。与障碍方程交换卷积（作者通用 `convolution_le_of_eq_zero`）得到，在 u=0 处几乎处处
\[
K_s*f=K_s*\sigma-\beta_s*u+a_su
\le\kappa\int K.
\]
取 \(\kappa=t/\int K\)，同一个例外集对每个尺度有效。先合并正有理尺度的零测集，再用 \(s\mapsto K_s\) 的 L1 连续性及 f 有界性得到所有正尺度。于是核极大算子的弱型常数不超过 \(\int K\)。一般非负 L1 输入由 \(f_n=\min(f,n)1_{[-n,n]^2}\uparrow f\) 单调传递；一般 f 替换为 |f|。阈值始终使用严格超水平集。

单位菱形面积为 2，K 在其上至少为 1，所以菱形平均极大算子的常数至多 \(\tfrac12\int K=U\)。线性映射 \((u,v)\mapsto(u+v,u-v)\) 将菱形映成轴平行方体，行列式为 -2；分子中的面积和分母 L1 质量同时乘 2，弱型常数不变。故结论适用于所有中心轴平行方体和所有尺度。


有限折线分布表示的完整证明见 [distribution_lemma.md](distribution_lemma.md)。
