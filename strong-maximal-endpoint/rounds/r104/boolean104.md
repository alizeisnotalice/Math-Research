# 二格乘积目录的全覆盖密度斜率

在单位环面的每个坐标等分两格。矩形目录的格集合为布尔立方体的全部子立方体，每个格内允许任意可测二值源。设 $T_D$ 为所有有私有顶点的子立方体族的最大归一化体积和。第101轮的私有格构造给出全覆盖的真实最小密度下确界 $1-\eta T_D$，其中 $0<\eta=1-p\le2^{-D}$。

写 $A_D=\max_{0\le s\le D}\binom Ds 2^{-s}$。存在绝对常数 $c>0$ 使
\[
 A_D\le T_D\le(1+O(e^{-cD}))A_D,\qquad A_D\sim \frac{3}{2\sqrt{\pi D}}(3/2)^D.
\]
下界用所有余维 $s$、指定坐标全部为1的子立方体。对应的恰含这 $s$ 个1的格是私有顶点。漏掉的点可用单格加入，保留每个原有私有点。

对任意不可约族和一个点 $x$，只看包含 $x$ 的盒 $R$。取其私有点 $u_R$，令 $A_R$ 为 $u_R$ 与 $x$ 不同的坐标，$B_R$ 为盒的固定坐标。$A_R\cap B_R=\varnothing$，而 $A_R\cap B_S\ne\varnothing$ 对不同盒成立，因为 $u_R\notin S$ 且 $x\in S$。随机排列全部坐标，事件“$A_R$ 全先于 $B_R$”两两互斥：不同两组的两个交叉交点必须互为前后。故 $\sum_{R\ni x}\binom{|A_R|+|B_R|}{|B_R|}^{-1}\le1$。该计数沿用 Ellis 对 Meshulam 子立方体界的证明。[Ellis，Irredundant Families of Subcubes](https://arxiv.org/abs/1003.2960)

在所有 $x$ 上求和。余维 $s$ 的盒恰贡献
\[
 C_{D,s}=\sum_{l=0}^{D-s}\frac{\binom{D-s}{l}}{\binom{s+l}s}
 =\frac{\sum_{j=0}^{D-s}\binom Dj}{\binom Ds}.
\]
从而体积和不超过 $\max_s 2^{D-s}/C_{D,s}$。在 $s\le5D/12$ 区域，分母中的二项和为 $2^D(1-O(e^{-cD}))$，来自对称二项尾的 Chernoff 界。在 $5D/12<s\le D/2$ 区域，二项和至少 $2^{D-1}$，而 $\binom Ds 2^{-s}=(3/2)^D\Pr(\mathrm{Bin}(D,1/3)=s)$，故这些项相对峰值指数小。在 $s>D/2$ 区域，只保留 $j=D-s$ 项即给上界 $2^{D/2}$。结合 Stirling 在二项众数处的展开得到结论。

该斜率刻画有限目录在极高门槛处的精确响应。连续强极大算子的临界亏损 $d_N/N$ 比 $2^{-N}$ 大，此处的线性区间不能直接代入该门槛。
