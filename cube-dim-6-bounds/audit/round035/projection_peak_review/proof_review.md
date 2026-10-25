# F035：投影峰 trace 公式的独立数学审查

## 1. 商空间的测度与峰值归一化

令 ψ(ξ)=Σ_(v∈body32)|v·ξ|，L 为 r=1,2,3 个独立 body 方向组成的 r×6 行矩阵，W=ker L，d=6−r，D=det(LLᵀ)>0。投影到 W 的实际 Cauchy 概率律的特征函数为 exp(−tψ(w))。选中方向的投影是0，所以它就是删去这些方向后的正交商律。ψ 在 W 的单位球面上有正下界，故 exp(−ψ) 可积。

在 W 的正交 Lebesgue 测度下，Fourier 逆变换给

    q_d(1,0)=(2π)^−d ∫_W exp(−ψ(w))dH^d(w)。

对任意 y∈W，三角不等式给 |q_d(1,y)|≤q_d(1,0)；实际概率律保证 q_d≥0，因此0确为全局峰值。正一阶齐次性给 sup q_d(t,·)=t^−d q_d(1,0)。这一步不涉及空间 Green 核的光滑性。

推前**完整六维**测度 exp(−ψ(ξ))dξ，经满秩线性映射 L 得密度

    F(z)=D^−1/2 ∫_W exp[−ψ(w+Lᵀ(LLᵀ)^−1z)]dH^d(w)。

coarea 因子准确是 D^−1/2。因此 ∫_W exp(−ψ)=√D F(0)，而不是除以 √D。对局部有界 z，ψ 的三角不等式给被积函数≤C exp(−c|w|)，支配收敛证明 F 在整个 R^r 连续，特别包括0。最后应使用

    q_d(1,0)=√D R_L/(2π)^d，κ_L=R_L/(2π)^d，R_L=F(0)。

因此033报告的两种常数测度方向正确，不能将 q_d 与 κ_L 混用。

## 2. fan 完整性与为什么不能在 W 上直接数面

独立程序不调用旧 fan verifier，而用以下充分证据重新证明21个基锥的完整 a.e. 覆盖：

- 每个整数 ray 归一到第一坐标1后位于 section 1≥a1≥…≥a5≥0，六个 rays 线性独立，故 section 是真正的五维 simplex。
- 每个 body form 在一个 simplex 的全部顶点取同一弱符号，因此整锥上 ψ 是线性函数。率 q_j 由32个 body forms 的绝对值准确重算。
- 对每个 simplex 的有理重心，32个 body forms 都非零。21个严格符号字互不相同，故 simplex 内部两两不相交。
- 每个 section simplex 的坐标体积为 |det R|/[5!∏r_(0,j)]。准确总和为1/120，等于 ordered unit section 的体积。因此其并覆盖 section 除零测集以外全部。

将这些锥施加720坐标置换及64独立坐标符号，得到整个 R6 的 a.e. 分割。坐标 ties、坐标零值及锥面都属于六维零测集。在完整六维测度层面，这些面不会造成重复质量。

但是 W 本身可以完全落在这些面上，所以不能直接把闭 fan 锥与 W 的交集体积相加。033 没有作此错误步骤；它先推前六维的正确测度等式，再在所有 cone densities 都连续的 generic ray 上求值，最后只用总 F 在0的连续性取极限。个别锥密度在0不连续无关紧要。

## 3. 一个锥的正交体 fiber 与消失阶

写一个全群锥 ξ=R u，u≥0；ψ=q·u，q_j>0，α=|det R|/∏q_j。代换 y_j=q_j u_j 后，原测度成为 α exp(−Σy_j)1_(y≥0)dy，法向映射是 A=L R diag(q_j^−1)。令 J 为其非零列、m=|J|。rank A=r，故 m≥r；零列变量的指数积分恰为1，不产生额外阶数或权重。

所选 L 每行是 body form，故在此锥上的六个系数不混合正负。逐行反号后 A 的全部列非负；每个 J 列又非零，所以存在行符号向量 s，使 s·a_j>0 对所有 j∈J。对 fiber A_J y=z、y≥0，

    Σ y_j≤(s·z)/min_(j∈J)(s·a_j)。

这给出 fiber 有界性；若右侧为负，fiber 为空。常数依赖锥但锥数有限。这个性质不可省略，否则 m>r 的非负 fiber 可能无界而抵消缩放幂次。

在 z=εv 上，m>r 的 coarea 密度为

    α ε^(m−r)/√det(A_J A_Jᵀ)
      ·∫_{A_J w=v,w≥0} exp(−εΣw_j)dH^(m−r)(w)。

fiber 固定且有界，积分一致有限，因此该锥贡献趋于0。若 fiber 实际维数更低，相应 Hausdorff 积分为0，结论仍成立。

当 m=r 时，M=L R_J 可逆，且 A_J=M diag(q_J^−1)。如果 M^−1v 的坐标严格为正，则沿此 ray 的极限为

    α/|det A_J|=|det R|/[|det M|∏_(j∉J)q_j]；

若至少一个坐标为负则贡献为0。正 q_J 仅改变坐标大小，不改变此符号判断。

## 4. generic ray 上连续性与 density identity 的逐点合法性

不能只凭“a.e.密度相等”便在一条零测 ray 上求值。独立检查如下。

对任意独立 r 列 B，消去其变量，其余 m−r 变量用普通 Lebesgue 积分。near z0=εv 时，前节 fiber 界给所有有效自由变量一个共同紧包络。被积函数包含 solved coordinates A_B^−1(z−A_free w)≥0 的指标。

如果某 solved 边界对 w 是非恒定仿射式，其零集在自由变量空间中测度为0。如果它与 w 无关且在 z0 消失，则对应线性泛函同时湮灭 r−1 个其他 basis 列及全部 free 列，因此 z0 位于某 r−1 个独立非零列张成的超平面。取 v 避开所有这种超平面，就排除了此例外。于是可对紧包络用支配收敛，得到每个 cone density 在 z0 邻域连续。m=r 时只是线性变换的正正交体指标，同样在避开其 facets 时连续。

总 F 连续，各 cone density 在同一 z0 连续，所以推前测度的 a.e. 等式在 z0 逐点成立。锥数有限，可以逐锥取 ε↓0 极限，得到

    R_L=Σ_(m=r,M^−1v>0) |det R|/[|det M|∏_(j∉J)q_j]。

没有额外 1/2，也不能未经论证将 opposite global signs 合成权重2。两侧一般选中相反 normal rays；本次全枚举明确保留所有64个符号。

## 5. 通用方向的准确安全性

所有整数 rays 的 l1 长度≤12。因此 M 的每个元素绝对值≤12。r=2 的 facet normal 可取一个列向量旋转，整数系数绝对值≤12；r=3 的 facet normal 可取两个列的叉积，整数系数绝对值≤2·12²=288。列的 q 缩放不改变张成超平面。

取 v=(1)、(1009,1)、(1009²,1009,1)。对 rank3 的非零整数法向 n，若首坐标非零，|n1|1009²>288(1009+1)；若首坐标为0且次坐标非零，|n2|1009>288；剩余情形显然。因此 n·v≠0。rank2 和 rank1 更简单。

这是覆盖**全部可能 facets**的先验证明，不依赖作者保存的 facet 列表是否穷尽。独立枚举另外检查所有 m=r solved coordinates 均非零，并核对保存的 facet normals 的幅度及非零点积。

## 6. 九类轨道与独立组合计数

rank1 只有一个轨道。rank2 在 body line 的任意代表反号、物理坐标符号和置换下，由 Hamming 距离 k 与6−k 合并，剩 k=1,2,3。相应不同 unordered pairs 数为96、240、160，总496。

rank3 将每个物理列按其首行符号归一，得到四个 projective classes：(1,1,1)、(1,1,−1)、(1,−1,1)、(1,−1,−1)。signed row permutations 准确诱导全 S4（24个不同作用）；坐标置换仅保留四类数目的 multiset。至少三类非空等价于秩3；任意三个上述不同 classes 线性独立。总数6、至少三个正数的四部划分恰为4110、3210、2220、3111、2211五种。

独立程序没有重跑4960 triples，而按组合公式计数：固定第一条 body line（32选），把它归一为全正；第一个物理坐标必为 class0。对某个有序 class 数 n0,...,n3，有 n0≥1 时，其余五坐标分配数为

    5!/[(n0−1)! n1! n2! n3!]。

对给定 multiset 的所有不同排列求和，再乘32并除3!，得到 unordered triples 数480、1920、480、640、1440，总4960，Gram 行列式依次64、96、128、160、192。这样同时核对了分类完整性与后续需要的真实 multiplicity。

在一个固定 normal representative 内，坐标列的重复类数为 n_c，每个不同 coordinate pattern 对应 ∏n_c! 个原720 permutations。我们的递归直接生成“带符号 class 多重集”的所有序列，没有先产生720 permutations再去重；符号是64个独立选项。最后只乘 ∏n_c!，不再乘 row orbit、triple count、global-sign2 或其他稳定子。

这与后续 B008 ordered tuple 权重是不同层次：本审计提供单个几何商的峰，后续 tuple sum 仍须按其自身真实 orbit count 和 ordering 支付。

## 7. 裁决与范围

上述证明、独立完整枚举和260保存项逐一核对共同支持033的九个 R、q_d、κ_L 作为后续严格常数来源。只做 selected-term recheck 本来不能排除遗漏，本次完整一次701568项复核补上了这一点。

该结论不提供任何 Green I、全域 g 或 C6 新界。新峰代入 B008 时必须保留商测度、各次 coarea Jacobian、ordered tuple 数及相应时间幂。不能把 Q3 直接当作 Q5，也不能把一个商轨道常数直接充当全部轨道的和。
