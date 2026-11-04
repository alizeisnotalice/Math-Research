# 第 20 轮复审：平面定量严格间隙

审查对象为 \(n=2,\beta=1\)、\(\alpha=1-t\)、\(0<t\leq10^{-2}\)。本记录给出一条不调用定量 Brunn--Minkowski 外部定理的自足路线。数值或符号计算只核对候选证书的恒等式；下述间隙证明本身是解析的。

## 1. 支撑函数的等号稳定性

令 \(A\) 为平面紧凸体，\(S=(A-A)/2\)。将 \(A\) 平移后，其支撑函数写为
\[
h_A(\theta)=e(\theta)+o(\theta),\qquad
e(\theta+\pi)=e(\theta),\quad o(\theta+\pi)=-o(\theta),
\]
其中平移消去了 \(o\) 的 \(k=1\) Fourier 模。此时 \(h_S=e\)。采用
\(o=\sum_{k\geq3,\ k\text{ odd}}(a_k\cos k\theta+b_k\sin k\theta)\)，支撑函数面积公式（先对平滑严格凸体成立，再以支撑函数卷积逼近）给出
\[
 |S|-|A|=\frac\pi2\sum_{\substack{k\geq3\\k\ {m odd}}}
 (k^2-1)(a_k^2+b_k^2).
\]
故由 Cauchy--Schwarz 和
\[
 \sum_{k\geq3,\ k\ {m odd}}\frac1{k^2-1}=\frac14
\]
有
\[
 \|o\|_\infty^2\leq\frac{|S|-|A|}{2\pi},
 \qquad d_H(A-b,S)\leq\sqrt{\frac{|S|-|A|}{2\pi}}
\tag{Q1}
\]
（\(b\) 是消去第一谐波所用的平移）。系数 \(\pi/2\)、求和及常数均已复核正确。这里必须明确 Fourier 系数采用上述实三角归一化；否则 \(\pi\) 因子会改变。

## 2. 内接对称体的面积--Hausdorff 估计

令 \(T=\overline K_\alpha/2\)，\(W=\alpha^{-1}-\alpha\)。\(T\) 包含顶点为 \(\pm We_1/2,\pm We_2/2\) 的菱形，并包含于 \([-W/2,W/2]^2\)。所以其内半径可取
\(\rho=W/(2\sqrt2)\)，外半径可取 \(R=W/\sqrt2\)，且 \(R/\rho=2\)。

若中心对称凸体 \(S\subset T\)，对任意单位向量记
\(d=h_T(u)-h_S(u)\)。取 \(T\) 的支撑点和内球上两个横向点，得到一个在 \(T\setminus S\) 内的三角帽，面积至少 \(\rho d^2/h_T(u)\)。取反向帽（除零测边界外不交），故
\[
 |T|-|S|\geq\frac{2\rho}{R}d^2\geq d^2,
 \qquad d_H(S,T)\leq\sqrt{|T|-|S|}.
\tag{Q2}
\]
这里用到了 \(h_S(u)\geq0\)（\(S\) 中心对称且含原点）；即使 \(h_S(u)=0\)，两帽只在零测底边相交。这一帽估计是有效的，不需要光滑性。

## 3. 对一个近极值光晕的合并

对质量一输入的光晕 \(H\)，置 \(A=\overline{\operatorname{conv}}H\)，并令
\(\delta=\mathcal U(\alpha)-|H|\)。两点体结论给出 \(S\subset T\)，而
\[
 (|S|-|A|)+(|T|-|S|)=\mathcal U(\alpha)-|A|\leq\delta,
 \qquad |A\setminus H|\leq\delta.
\]
将 (Q1)--(Q2) 对这两个不同的缺口作 Cauchy--Schwarz 合并，得到
\[
 d_H(A-b,T)
 \leq\sqrt{\left(1+\frac1{2\pi}\right)\delta}
 <\sqrt{\frac43\delta}.
\tag{Q3}
\]
这里不能把两个 \(\sqrt\delta\) 项直接相加后仍声称同一常数；所用的是
\(\sqrt x+\sqrt{y/(2\pi)}\leq\sqrt{(1+1/(2\pi))(x+y)}\)。

## 4. 显式四点证书及传递

取
\[
r=\frac t2+\frac{t^2}{2}+\frac{15t^3}{16}+\frac{9t^4}{4},\quad
q=\alpha+\frac{t^2}{4},\quad \varepsilon=\frac{t^3}{16}.
\]
在每个中心作 \(\ell^\infty\) 半径 \(\varepsilon\) 的独立扰动。对二次式的扰动上界为
\(4\varepsilon(z-2r)+4\varepsilon^2\)，四个小方体并的面积损失上界为
\(16\alpha\varepsilon\)。实际运行 WolframScript 15.0.0（2026-11-05）所得精确余量为
\[
\begin{aligned}
&\frac{t^5}{32}(704+415t+576t^2+648t^3),\\
&\frac{t^4}{16}(23-330t-157t^2-252t^3-324t^4),\\
&\frac{t^3}{16\alpha^2}(24+64t-232t^2+135t^3-66t^4-13t^5+396t^6-324t^7),\\
&\frac{t^3}{64}(304+640t+240t^2+801t^3+1080t^4+1296t^5),
\end{aligned}
\]
分别对应 \(G_r(\alpha)\)、\(-G_r(q)\)、\(-G_r(U)\) 与质量矛盾在扣除上述扰动损失后的余量。这些多项式在 \(0<t\leq10^{-2}\) 上为正；它们也保证 \(r<r_*\)、\(2r<\alpha\)、\(\alpha<q<(\alpha+U)/2\)。

第一余量还意味着每个证书 \(\ell^\infty\) 球包含于 \(T\)：对球内点同时取其负点（落在对角球内），严格的对角两点条件给出 \(2x\in K_\alpha\)。这一步宜在正文中明说，不能只从中心在 \(T\) 内推出整球在 \(T\) 内。

若 \(\delta<\varepsilon^2/3\)，(Q3) 给
\(d_H(A-b,T)<2\varepsilon/3\)。每个证书球中半径 \(\varepsilon/3\) 的欧氏球遂包含于 \(A-b\)：对其中任一点，半径 \(2\varepsilon/3\) 的球仍在 \(T\)，比较支撑函数即可。它的面积是 \(\pi\varepsilon^2/9>\delta\)，故 \(H-b\) 必须与每一个这样的球相交。这给出同一个输入的四个证书邻域光晕点，矛盾。因此
\[
 \boxed{\quad \mathcal U(1-t)-C_{1,\mathcal Q_c}(1-t)
 \geq\frac{\varepsilon^2}{3}=\frac{t^6}{768},\qquad 0<t\leq10^{-2}.\quad}
\]
该六阶下界只是一条可证的非尖锐定量间隙，不蕴含或改变已有三阶区间。

## 结论

路线成立，关键是将 (Q1) 与 (Q2) 的两个面积缺口先分开、再作 Cauchy--Schwarz 合并。正文若采用，应写出内球三角帽、球属于 \(T\) 的对角配对理由，以及 Fourier 归一化/逼近说明；无需援引复杂的定量 Brunn--Minkowski 文献。

## 正式稿终审（`pb:fourier-stability`--`pb:sixth-gap` 与 `mom:sixth-envelope`）

通过，未发现阻塞性漏洞。

1. `pb:fourier-stability` 的 Fourier 归一化、\(\pi/2\) 系数与
\(\sum_{k\ge3,\,k\text{ odd}}(k^2-1)^{-1}=1/4\) 一致，故
\(d_H(A-b,S)^2\leq\Delta/(2\pi)\) 正确。\(h_\eta=1+\eta\cos3\theta\)
确为支撑函数；\(h_\eta+h_\eta''>0\)，\(|S|-|A_\eta|=4\pi\eta^2\)，
且第一谐波平移不能消去第三谐波，足以说明平方根指数在一般凸体类中不能改善。

2. `pb:cap-distance` 的两三角帽互相只可能在零测边界相交，故其系数
\(R/(2\rho)\) 无误。对 \(T=\overline K_\alpha/2\) 的菱形内含和正方形外包
给出 \(R/(2\rho)=1\)。

3. `pb:quant-transfer` 正确地以
\(g+\Delta\leq\delta\) 作 Cauchy--Schwarz 合并；这避免了把两个平方根误加的常数错误。
距离小于 \(2\varepsilon/3\) 时，\(T\) 内半径 \(\varepsilon\) 球的同心
\(\varepsilon/3\) 球确实在 \(A-b\) 内：后者每点仍携带半径
\(2\varepsilon/3\) 的 \(T\)-内球，可逐方向比较支撑函数。孔洞面积与
\(\pi\varepsilon^2/9\) 的比较也严格。

4. 六阶证书的四个扣除扰动后的多项式恒等式被 WolframScript 15.0.0 独立重算；
它们与正文 \(P_1,\ldots,P_4\) 一致，区间下界足以保证严格性。第一余量配对
\(x\) 与 \(-x\) 说明每个证书闭球确在 \(T\) 内，因而可调用传递命题；
于是 \(t^6/768\) 的结论正确且未被表述为尖锐阶。

5. 高矩部分归一化正确：
\[
 \frac p{768}B(p,7)=\frac{15}{16}\prod_{j=1}^6(p+j)^{-1}.
\]
实际 Wolfram 化简返回 `True`。截掉 \([0,99/100]\) 的误差至多
\((99/100)^p/768\)，所以给出的 \(O((99/100)^p)\) 与极限均正确。该项只改善统一上包络，正文没有误称为 \(\mathfrak F\) 的真实六阶系数。

唯一可选的表述加固：`pb:quant-transfer` 中若担忧“非空光晕”却非满维，可把末段“空光晕”改为“零测光晕”。在本文严格光晕为开集的既有可测性结论下，非空光晕自动正测，故当前证明逻辑仍成立。
