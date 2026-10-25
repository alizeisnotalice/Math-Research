# 005 — 非顶点谱 ray 的 generic star 系数严格消失（n=5）

## 结论和适用范围

令
\[
 \psi(\xi)=\sum_{\epsilon\in\{-1,1\}^4}
 \left|\xi_1+\epsilon_2\xi_2+\cdots+\epsilon_5\xi_5\right|,
 \qquad P=\{\xi:\psi(\xi)\leq1\}.
\]
在已核的 7 个 body-fan 基锥及其全部符号/排列像中，令一个 ray
为 \(r\)，令 \(q=\psi(r)\)，并写 \(p=r/q\in\partial P\)。若 \(p\)
不是 **\(P\) 的真顶点**，那么在下列 generic 区域中，其简单 \(t\)-极点的
完整 star 系数恒为零：
\[
 S_p(x)=\sum_{C\ni p}\frac{|\det R_C|/\prod_{r\in C}\psi(r)}
 {\prod_{u\in C,\;u\ne p}(\langle x,p\rangle-\langle x,u\rangle)}=0.
\tag{1}
\]
这里和式按出现该归一 ray \(p\) 的全维锥求和；公共 Fourier 因子
\((2\pi)^{-5}\) 未写入。generic 的含义是：所有本 fan 中不同的归一 rays
的 \(\langle x,\cdot\rangle\) 不相等，且 \(x\) 不落在任一 pulling
单纯形两顶点之差的有限个消失超平面上。该区域是非空 Zariski 开集。

这不是“在一个数值点看见零”的推断。下面的函数恒等式先给出无人工
极点的顶点表示，再由有理函数的恒等性推出 (1)。因此可免去每个 giant
symbolic star 的通分。

## 自足的 Fourier--Laplace 论证（不把 Brion 当黑箱）

定义（不含 Fourier 归一化）
\[
 K(t,x)=\int_{\mathbb R^5}e^{-t\psi(\xi)+i\langle x,\xi\rangle}\,d\xi,
 \qquad \Re t>0.
\]
因 \(\{\psi\le s\}=sP\)，Stieltjes 分层和一次分部积分给出
\[
 \begin{aligned}
 K(t,x)
 &=\int_0^\infty e^{-ts}\,d\left[s^5
       \int_Pe^{is\langle x,y\rangle}\,dy\right]\\
 &=t\int_0^\infty e^{-ts}s^5
       \int_Pe^{is\langle x,y\rangle}\,dy\,ds. \tag{2}
 \end{aligned}
\]
边界项为零：\(F(s)=O(s^5)\)，而 \(e^{-ts}F(s)\to0\)；Fubini 由
\(\Re t>0\) 和 \(P\) 紧性成立。

这里可完全不用连续 Brion 定理。任一满维多面体都有只使用其原始顶点的
pulling 单纯剖分：递归选一顶点，并将不含它的边界面已经得到的单纯剖分
与该顶点锥化。令 \(\mathcal T\) 是这样的 5-单纯形剖分。对
\(\Delta=\operatorname{conv}(v_0,\ldots,v_5)\in\mathcal T\)，当
\(\lambda_j=\langle z,v_j\rangle\) 两两不同时，标准单纯形上的一维
divided-difference 递推（等价于逐变量积分）给出
\[
 \int_\Delta e^{\langle z,y\rangle}dy
 =|\det(v_1-v_0,\ldots,v_5-v_0)|
 \sum_{j=0}^{5}{e^{\lambda_j}\over
       \prod_{k\ne j}(\lambda_j-\lambda_k)}. \tag{3}
\]
系数没有额外符号：在一维 \([v_0,v_1]\) 时右边正是
\((e^{\lambda_1}-e^{\lambda_0})/(\lambda_1-\lambda_0)\)；五维式由同一
递推得到。注意 \(|\det|=5!\operatorname{vol}(\Delta)\)，故也没有遗漏
\(5!\) 因子。

将 (3) 对 \(\mathcal T\) 求和并按顶点收集，得到
\[
 E_P(z):=\int_Pe^{\langle z,y\rangle}dy
   =\sum_{v\in\operatorname{Vert}(P)}e^{\langle z,v\rangle}A_v(z),
 \qquad A_v(sz)=s^{-5}A_v(z). \tag{4}
\]
其中 \(A_v\) 是含 \(v\) 的单纯形所给有理系数之和，所以齐次次数恰为
\(-5\)。此式只按**真顶点**求和，而且其证明只用了有限剖分和 (3)。
所需条件在此均满足：\(P\) 是满维紧有理多面体（\(\psi\) 为有限个整
线性型绝对值之和，且这些型张成 \(\mathbb R^5\)），而所选 \(x\) 避开
有限个顶点差 \(v-v'\) 的消失超平面。

将 \(z=isx\) 代入 (4) 并代入 (2)。齐次性恰好消去 \(s^5\)：
\[
 K(t,x)=t\sum_{v\in\operatorname{Vert}(P)}
       \frac{A_v(ix)}{t-i\langle x,v\rangle}. \tag{5}
\]
初始等式在 \(\Re t>0\) 成立；两边都是 \(t\) 的有理函数，故 (5) 作为
亚纯恒等式延拓。特别地，\(K(\cdot,x)\) 的有限 \(t\)-极点只能属于
\(i\langle x,v\rangle\)，其中 \(v\) 是 \(P\) 的真顶点。

另一方面，直接对每个 body-fan 锥 \(C=\operatorname{pos}(r_1,\ldots,r_5)\)
积分，得到
\[
 K(t,x)=\sum_C\frac{|\det R_C|}{\prod_{j=1}^5(tq_j-i\langle x,r_j\rangle)}
 =\sum_C\frac{|\det R_C|/\prod_jq_j}
  {\prod_{j=1}^5(t-i\langle x,p_j\rangle)},\quad p_j=r_j/q_j. \tag{6}
\]
若一锥内五个 \(b_j=\langle x,p_j\rangle\) 两两不同，\(d=5\) 的简单
partial fraction residue 是
\[
 \mathop{\rm Res}_{t=ib_j}\prod_{k=1}^5(t-ib_k)^{-1}
 =\prod_{k\ne j}\frac1{i(b_j-b_k)}
 =\frac1{\prod_{k\ne j}(b_j-b_k)}, \tag{7}
\]
因为 \(i^{-4}=1\)。这正是题述的 \(\alpha/\prod(b_r-b_s)\) 系数，
\(\alpha=|\det R_C|/\prod q_j\)，并非少了符号或 \(i\) 因子。

现在固定一个非顶点 ray \(p\)。generic 性确保 \(\langle x,p\rangle\)
不等于任何顶点投影。由 (5) 它不是 \(K\) 的极点；由 (6)--(7) 取留数
即得 (1)。左边是 \(x\) 的有理函数，在非空开集为零，故在其定义域上
恒为零。重复节点只能在这个 generic 恒等式证明之后取 Hermite/合并
principal-part 极限，绝不能把单锥的分母为零公式直接代入。

## 本实例中哪些 ray 是人工的

`005_vertex_rays_exact.py` 独立枚举 \(2^{16}\) 个符号半空间
\[
 P=\bigcap_{\sigma\in\{-1,1\}^{16}}
 \{y:(\sum_v\sigma_vv)\cdot y\le1\}.
\]
对每个基锥的归一 ray \(r/\psi(r)\)，活跃半空间法向量满秩当且仅当该点
是顶点。11 个符号/排列轨道代表中，仅
\[
 (1,1,0,0,0)/16,\quad(1,1,1,1,0)/24,\quad
 (2,1,1,1,1)/36
\]
是顶点；其余 8 个代表的活跃法向量秩为 \(1,1,1,4,4,3,2,2\)（详细输出
见 JSON）。故本审计所称“人工”恰是那 8 个非顶点轨道，不是凭肉眼
将 star triangulation 的内部 ray 认作非顶点。

## 实际精确检查（仅作一致性检查）

`005_star_residue_specialization.py` 在
\(x=(1,10^4,10^8,10^{12},10^{16})\) 用 `Fraction` 枚举 26,880 个全锥、
134,400 个简单极点贡献。该点分离全部 1,762 个归一 rays：1,482 个
人工 ray 的完整 star 留数都严格为零，280 个真顶点 ray 在该点都非零。
这是对 (6)--(7)、符号和锥计数的精确交叉检查；**它本身不被用作一般
证明**，一般性来自 (2)--(6)。

## 原始文献与边界

* Michel Brion, *Points entiers dans les polyèdres convexes*, Ann. Sci.
  Éc. Norm. Supér. (4) 21 (1988), 653--663,
  DOI [10.24033/asens.1572](https://doi.org/10.24033/asens.1572).
  Theorem 2.2 is the vertex-tangent-cone rational identity; §2.4 explicitly
  permits subdivision of tangent cones. 它是上述顶点局部化观点的背景；本
  报告的实际证明采用 pulling 单纯剖分和公式 (3)，不依赖连续极限版本。
* Jim Lawrence, *Polytope volume computation*, Math. Comp. 57 (195)
  (1991), 259--271, DOI
  [10.1090/S0025-5718-1991-1079024-2](https://doi.org/10.1090/S0025-5718-1991-1079024-2).
  这是顶点局部化/体积路线的原始公开来源；本报告不把它误报为本算子
  上界或热核估计。

该结果只消除了 generic star 的非顶点极点；它不自动给出 \(G\) 的正性、
截断误差、Green 势下界或弱型常数上界。
