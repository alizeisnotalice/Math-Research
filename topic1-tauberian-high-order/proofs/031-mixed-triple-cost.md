# R31：六阶饱和下混合三临界边的额外源成本

状态：完整新引理及证明，供必要独立复核；尚未据此修改论文。只使用已审 R23 的实际源／目标和六阶价格，不引用尚待独审的 R28 七阶轮廓或价格。

## 结论与范围

固定 \(n\ge3\)。沿 R23 同一六层构造的任意输入序列，假设
\[
1+M/V\ge L_n(t)+[K_n-o(1)]av^6,
\qquad K_n=\binom n2n^2/27,quad a=v^n.
\tag{S}
\]
原末层量沿用 R23：实际源 \(F_i\)、完整本目标 \(G_i\)、质量根 \(R_i\)、原／有效／中点盒、完整盈余 \(J_i\)，以及
\[
\sum_i\varepsilon_i r_i^n=O(v^5S),\qquad
\sum_{\rm deep}(\alpha_i+\Delta_i^{\rm end}/r_i)r_i^n=O(v^6S).
\tag{1}
\]

对每条恰有三个临界轴的深处边 \(e=\{i,j\}\)，选 \(R_i\ge R_j\) 的一端为所有者；在这三个轴内排序交长为 \(L_1\le L_2\le L_3\)。令 \(I_e=|X_i\cap X_j|\)。则对每个固定 \(T<\infty\)，有
\[
\boxed{\quad
\sum_{\substack{q(e)=3\\L_2\le T v^3r_i}}I_e=o_T(v^7S).
\quad}
\tag{2}
\]
第三交长只有 R14 的原上界 \(L_3\le C_Avr_i\)。因此 (2) 正覆盖所担心的 \((v^3,v^3,v)\) 混合厚度。它不声称全部三临界边均为 \(o(v^7S)\)，不处理两窄轴超出每个固定 \(Tv^3\) 窗口的尾部，也不证明高维 Slot。

要点是：三临界片不只占用一个已付二轴源臂。其目标真实质量与原二轴全部目标相加后，会把缺源层再抬高。只计算这一抬高部分，就能在不重复支付已付源点的条件下获得 (2)。

## 一、可独立从 R23 得到的两条饱和输入

写实际末层有符号角能量为
\[
H=\sum_{i<j}|X_i\cap X_j|-\sum_iJ_i.
\]
在 R23 最终并集及六层望远镜中，末层暂保留精确的 \(H\)，而不代入其上界。小源质量成本为 \(O(v^7S)\)，环境为 \(O(av^6S)=o(v^6S)\)，早五层总权为 \(O(v)\)，初端为 \(O(v^7S)\)。同一内部图项可丢掉非负剩余，故
\[
V+M\ge A_vM/a-H/(1+O(v))-o(v^6S).
\]
由 (S)、\(V\ge(1-a)S\)，得到 \(\liminf H/(v^6S)\ge K_n\)。R23 已证上界给反向上极限。因此
\[
H/(v^6S)\longrightarrow K_n.
\tag{3}
\]
此处只需六阶精度，不需要七阶账本，也不需要 R28 已称的三角轮廓。

其次，固定 \(0<\sigma<1/4\)，采用 R23 大角支付阈值 \(T_\sigma\)，可将其增大至任意预先给定常数。对小的良好二临界边，用 R23 同源全部目标相加的实际强迫式构造禁层 \(Y_{i,\sigma}\subset(x_i+W_i')\setminus F_i\)。在一组有向轴对的第一排序中，使用统一的高度
\[
h_i(y')=
\left[\frac{(1-\eta_\sigma(v))\sum_jP'_j(Q'_j-y')_+}{nv^2}
       -E_i\right]_+,
\quad \eta_\sigma(v)\to0,
\tag{4}
\]
\[
E_i=C_\sigma\left(v^5+\alpha_i+\Delta_i^{\rm end}/r_i
                                      +\sum_j\varepsilon_j\right).
\]
所有其他坐标离有效盒两端大于 \(H_\sigma v^3r_i\)。第二排序及其他有向轴对按 R23 取并；这些层互不交。高度积分的总误差为 \(O_\sigma(v^7S)\)。常数可取大，使构造和下文共同适用的全部高度窗口均包含在内。

记未用的有效源缺量成本为
\[
Z_\sigma=\frac va\sum_i
 |[(x_i+W_i')\setminus F_i]\setminus Y_{i,\sigma}|.
\]
则存在不依赖 \(T_\sigma,\sigma\) 的 \(c>0\)，使
\[
\boxed{\quad
\limsup_{v\to0}\frac{Z_\sigma}{v^6S}\le C\sigma.
\quad}
\tag{5}
\]

这里的证明不需要形状稳定性。CDF 的精确不等式给
\[
|W_i|-\mu_i\ge |(x_i+W_i')\setminus F_i|
=|Y_{i,\sigma}|+
 |[(x_i+W_i')\setminus F_i]\setminus Y_{i,\sigma}|.
\]
良好源上体积导数为 \((1+o_\sigma(1))v/a\)，全部源上有统一 \(cv/a\) 下界；非良好源取 \(Y=\varnothing\)。故在保留原禁层系数的同时，还可保留余下缺量的固定正系数。R23 用 \(\sigma\sum J_i\) 支付大边，在 \((1-\sigma)\sum J_i\) 中同时留下这一余量，得到
\[
H+c(1-\sigma)Z_\sigma
\le [K_n/(1-\sigma)^2+o_\sigma(1)]v^6S.
\]
三临界及边界旧粗交量只需为 \(o(v^6S)\)，在这里已足够。与 (3) 合用即得 (5)。导数和 CDF 的剩余系数只依赖原固定 \(n,A\)，不依赖角长阈值，因此 (5) 的 \(C\) 可独立于 \(\sigma\) 选择。

## 二、resolved 混合片的真实附加捕获

固定 \(\tau>0,T<\infty\)，先考虑
\[
\tau v^3r_i\le L_1,L_2\le Tv^3r_i,
\qquad \tau vr_i\le L_3\le C_Avr_i.
\tag{6}
\]
丢掉所有者非良好的边，只损失 \(o(v^7S)\)：这些混合片各自至多为 \(C_Tv^7r_i^n\)，相交度有界。还可丢掉
\(\Delta_i^{\rm end}>\tau v^3r_i/100\) 的所有者；由 (1) 其权为 \(O_\tau(v^3S)\)，交量为 \(O_{\tau,T}(v^{10}S)\)。后文仅处理余下边。

反射三个坐标使所有者 \(i\) 在目标 \(j\) 的同一侧。因 \(R_i\ge R_j\)，三个原始面对距离均至多为 \(R_i-L_l\)。在有效源盒中取
\[
d_1=s,\qquad d_2=s+y,
\quad 0\le y\le cL_2,
\qquad cL_3\le d_3\le2cL_3,
\tag{7}
\]
其中 \(c>0\) 只依赖原宽度上下界，取足够小；其他 \(n-3\) 个坐标在离两端 \(H_\sigma v^3r_i\) 的中心区。允许 \(s\le C_{\sigma,T}v^4r_i\)。固定参数后取小 \(v\)，则前两轴各至少捕获固定倍数的 \(L_1,L_2\)。第三轴实际交长为一个真正的一维帽截长，下界为
\[
\min\{w_{j,3}^0,(L_3-d_3-\Delta_i^{\rm end})_+\}
\ge\gamma L_3,
\tag{8}
\]
其中 \(\gamma>0\) 固定。理由是 \(d_3\le L_3/4\)、\(\Delta_i^{\rm end}\le L_3/4\)，而 \(w_{j,3}^0\ge c_wvr_i\)、\(L_3\le C_Avr_i\)。这不是将第三轴强行乘以 \(b_e\)，而是在显式实际源带中逐点保留其真实目标交长。

非临界方向的整个原始目标宽度都进入同一见证。扣目标原始盒缺量后，三临界实际目标给额外捕获下界
\[
c(vr_i)^{n-3}L_1L_2L_3-C\varepsilon_j ar_i^n.
\tag{9}
\]
负的下界允许保留，稍后按总缺量支付。

同时，第 1、2 轴这一有向组内的所有二临界目标，仍被同一个点的同一个原见证按 R23 捕获。三临界目标 \(j\) 不在该列表中，所有这些 \(F\) 两两不交且与 \(G_i\) 不交。因此必须先把 (9) 加到原列表，再用第 1 轴的主面对容量。令 \(p_l=L_l/r_i\)，得到实际源点所需高度
\[
\frac{s}{r_i}\ge
\frac{(1-\eta_\sigma(v))\sum_kP'_k(Q'_k-y')_+}{nv^2}
 +c\frac{p_1p_2p_3}{v^3}-E_i-C\varepsilon_j.
\tag{10}
\]
可取与 (4) 相同的 \(\eta_\sigma(v),E_i\)，将二者常数统一增大即可。目标数和高度窗口均固定有界；容量、截深等 \(O(v)\) 乘法误差纳入同一个 \(v^5\) 项。

重要的是，(10) 在原 12 禁层以上仍然强迫缺源，而不只是再次数到原层。第三深度在 \(\asymp L_3\asymp v r_i\) 的内部带，离该轴两个端点均远大于 \(H_\sigma v^3r_i\)，所以本带不遇任何含第 3 轴的基准二轴层。其他轴也在中心区，因此此处 \(Y_{i,\sigma}\) 只有 (4) 的第 12 组第一排序。

对非负 \(h,b,E,D\)，
\[
(h+b-E-D)_+-(h-E)_+\ge b-E-D.
\tag{11}
\]
把 (10) 与 (4) 比较，积分源带 (7)，即得额外有效缺量成本
\[
\frac va|[(x_i+W_i')\setminus F_i]\setminus Y_{i,\sigma}|
\ge c r_i^n\frac{p_1p_2^2p_3^2}{v^5}
 -C_{\sigma,T}r_i^n v^2(E_i+\varepsilon_j).
\tag{12}
\]
左侧包含整个剩余缺量，只用于该边的一次下界；不宣称不同边的带互不交。其他中心区宽度为固定倍数 \(vr_i\)，在此只需粗正系数。

真实交量满足 \(I_e\le C r_i^np_1p_2p_3\)，而 (6) 给 \(p_2p_3\ge\tau^2v^4\)。因此 (12) 的主项至少为 \(c_\tau I_e/v\)。每源拥有的相交边数有统一上界，故可按这个次数分配左侧，不受源带重叠影响。由 (1) 及有界度，右侧误差求和为 \(O_{\sigma,\tau,T}(v^7S)\)。得到
\[
\boxed{\quad
\sum_{e\text{ 满足 }(6)} I_e
\le C_{\tau,T}v Z_\sigma
        +O_{\sigma,\tau,T}(v^8S)+o(v^7S).
\quad}
\tag{13}
\]
主系数 \(C_{\tau,T}\) 不依赖 \(\sigma\)：\(T_\sigma,H_\sigma\) 只影响高度及中心区开始适用的小 \(v\) 阈值和误差常数，不影响 (8)、(9) 的固定比例。先取 \(v\to0\)，用 (5)，再取 \(\sigma\downarrow0\)，(13) 给 resolved 边总交量 \(o_{\tau,T}(v^7S)\)。

## 三、去掉 resolved 下界

在 \(L_1,L_2\le Tv^3r_i\) 的固定窗口内，如果 \(L_1<\tau v^3r_i\)、\(L_2<\tau v^3r_i\) 或 \(L_3<\tau vr_i\)，则单边
\[
I_e\le C_T\tau v^7r_i^n.
\]
有限相交度使全部退化边总量不超过 \(C_T\tau v^7S\)。先固定 \(\tau\)，用上一节消去 resolved 部分，最后 \(\tau\downarrow0\)，即得 (2)。

## 证据和后续边界

Terra 在本机 Mathematica 15.0.0 实跑 `r031-computation/cap_convolution.wl`，日志 `cap_convolution.log`。实际读取记录：归一化的 (8) 违反集 False，正部差 (11) 违反集 False；完整一维帽卷积三个分支也已算出。这里只用局部带的点态下界，不需要把完整卷积当作强制的全宽系数。计算不替代同一实际源点的加法和额外缺层证明。

本引理把一个特定高维障碍压成严格的额外成本问题：已付二轴缺层之外仍有 \(I_e/v\) 的源成本，六阶饱和不允许它在总量上保持正比例。它不证明全体三临界尾部的七阶可积性，不证明 Slot，也未消除 R26/R28 的真实目标反馈修正。所有输出只在独立 scratch，仓库及他人源码均未修改。
