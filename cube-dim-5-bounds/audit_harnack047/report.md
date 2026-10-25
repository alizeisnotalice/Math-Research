# audit_harnack047：043 正混合 Harnack 接口独立审计

## 结论与验证范围

**通过，无阻断 046 的理论问题。** 043 的全空间 log-Hessian sandwich、全空间 log 抛物面，以及保存中心的四维闭盒/椭球证书均通过本轮独立审查。支配微分不需要剩余 Cauchy 向量的矩；无限时间混合归一后是合法概率；方差项及标量最大值 1/4 的方向正确。

本轮没有重新计算中心 H、梯度、Green、对数、指数或 Gamma。中心值与梯度的证明状态是：冻结源码的数学算法通过读审，冻结运行记录与 checkpoint 完整性通过字节核对；区域端点及剩余量另经一次纯 Fraction/SymPy 批独立核验。该批 PID **46260**，0.0735 秒，14 组检查全部通过，退出码 0。它不是第二次生产 043 的 3504 锥数值。

数值证书只在 **x1=1** 仿射平面上，使用四个物理尾坐标梯度。五维 Hessian 定理不意味着缺少的 x1 中心导数已经计算。一个中心的安全区域、或最多 128 个中心的部分覆盖，都不能称为全域 g。

只读来源目录：
`/Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper/cube-dim-5-bounds/upper_harnack043/`。

| 文件 | 实际核对的 SHA-256 |
|---|---|
| report.md | a298ca56ba980fce4fe2542929e03befb97b528afd3d390ed9d15651da68322a |
| certify.py | b463b1f941f47de2063e24ffaf6b4b70e1ffe58654186026fd71d27b9abef82f |
| results.json | 65ae80d993d3b207b54050a6ebf16e701e69234a02e43cb4aed245df25f4727d |

## 1. 正 Cauchy 卷积与积分下微分

固定任意可逆 body 基 N，将 16 个独立 Cauchy 变量中的相应五个作为基变量，其余 11 个在时间 1 的和记为 Y。Y 是有限实向量几乎处处，但无需密度或任何正矩。尺度性质给

\[
 q_t(z)=\frac1{|\det N|}\prod_{j=1}^5\frac{t}{\pi(t^2+(N^{-1}z)_j^2)},\qquad
 p_t(x)=E q_t(x-tY).
\]

同一个时间无关概率空间上的 Y 足够表示所有 p_t；并不声称所选耦合是整个随机过程的真实时间路径。

令 a=N^{-1}u、w=N^{-1}(x-tY)，沿固定空间方向 u 求导。由于 Y 的分布不依赖 x，导数中没有 Y 因子。score 与 log 二阶导数为

\[
 \sigma_u=\sum_j\frac{-2a_jw_j}{t^2+w_j^2},\qquad
 \tau_u=\sum_j\frac{2a_j^2(w_j^2-t^2)}{(t^2+w_j^2)^2}.
\]

对 t>0，\(2|w|/(t^2+w^2)\le1/t\)，可由 \((t^2+w^2)^2-4t^2w^2=(t^2-w^2)^2\ge0\) 直接核对。令 r=w^2/t^2>=0，则

\[
 f(r)=\frac{2(r-1)}{(1+r)^2},\quad
 f(r)+2=\frac{2r(r+3)}{(1+r)^2}\ge0,\quad
 \frac14-f(r)=\frac{(r-3)^2}{4(1+r)^2}\ge0.
\]

所以最小值 -2 在 r=0 取到，最大值 1/4 在 r=3 取到。因而

\[
 |\sigma_u|\le\|a\|_1/t,\qquad
 -2\|a\|_2^2/t^2\le\tau_u\le\|a\|_2^2/(4t^2).
\]

再由 \(\partial_uq_t=q_t\sigma_u\)、\(\partial_u^2q_t=q_t(\tau_u+\sigma_u^2)\) 及
\(0<q_t\le(|\det N|\pi^5t^5)^{-1}\)，分别得到与 x、Y 无关的 t^{-5}、t^{-6}、t^{-7} 支配。二阶混合偏导可极化得到同样的 t^{-7} 控制。对任意 epsilon>0，这些函数在 [epsilon,infinity) 上均可积；其概率空间积分再乘以 1。

逐点导数连续与支配收敛给 H_epsilon 属于 C^2(R5)，并允许两次空间积分下求导。该论证可直接对导数差商使用平均值公式及相同统一支配，避免只凭“积分了逐点导数”断言可微。严格正 q_t 和有限时间积分上界同时给
\(0<H_\epsilon(x)<\infty\) 对所有 x 成立，不受完整 Green 的奇异集影响。

## 2. 无限混合测度、方差及全空间 log 抛物面

原测度 dt dP(Y) 在无限时间区间总质量无穷，但它是 sigma-finite。正可积密度

\[
 d\nu_x=\frac{q_t(x-tY)}{H_\epsilon(x)}\mathbf1_{t\ge\epsilon}\,dt\,dP(Y)
\]

积分恰为 1，故定义概率测度。score 在该测度上有统一界 \(|\sigma_u|\le\|a\|_1/\epsilon\)，tau 也有界。因此以下相减均为有限量：

\[
 u^T\nabla^2\log H_\epsilon(x)u
 =\frac{\partial_u^2H_\epsilon}{H_\epsilon}
  -\left(\frac{\partial_uH_\epsilon}{H_\epsilon}\right)^2
 =E_{\nu_x}\tau_u+\operatorname{Var}_{\nu_x}\sigma_u.
\]

这里是有界 score 的方差，不是 Y 的方差。下界使用 Var>=0；上界使用
\(\operatorname{Var}\sigma_u\le E\sigma_u^2\le\|a\|_1^2/\epsilon^2\)，再用 \(\|a\|_1^2\le5\|a\|_2^2\)。得到准确矩阵意义的 sandwich

\[
 -\frac2{\epsilon^2}N^{-T}N^{-1}
 \preceq\nabla^2\log H_\epsilon(x)
 \preceq\frac{21}{4\epsilon^2}N^{-T}N^{-1}. \tag{1}
\]

方向更细的上界 \((\|a\|_2^2/4+\|a\|_1^2)/\epsilon^2\) 也成立。所有 x 和 u 使用同一个固定 N；没有在路径上切换基而遗漏控制。

沿 x0+s u 积分，二阶余项权重积分为 1/2，因此

\[
 \log H(x_0)+\nabla\log H(x_0)\cdot u
       -\|N^{-1}u\|_2^2/\epsilon^2
 \le\log H(x_0+u)
 \le\log H(x_0)+\nabla\log H(x_0)\cdot u
       +21\|N^{-1}u\|_2^2/(8\epsilon^2). \tag{2}
\]

这对任意有限长度 u 成立；没有必须留在投影互异区、远离 S1/S2 或小邻域的条件。因为 G>=H_epsilon，任何 H_epsilon 的区域下界也是同一区域 G 的下界，无须估计被删去的 (0,epsilon) 时间段。

## 3. 数值固定基与 face 几何

043 登记的五个 body 列给 det N=-48，与保存 N^{-1} 完全一致。限制 u=(0,v) 时

\[
 \|N^{-1}u\|_2^2=v^TQv,\quad
 Q=I_4/4-J_4/36,\quad Q^{-1}=4I_4+(4/5)J_4.
\]

Q 的特征值为 1/4（三重）与 5/36（一重），所以严格正定。对四维符号角点 sigma，

\[
 c_\sigma=\sigma^TQ\sigma
 =1-(\textstyle\sum_j\sigma_j)^2/36\in\{5/9,8/9,1\},
\]

最大值恰为 1。上述矩阵、特征值与全部 16 个符号点均由本轮准确核对；没有搜索新基。

数值来源只计算 x2,x3,x4,x5 的物理导数。式 (2) 在 x1=1 上只需这四个导数，故局部 face 证书完整；要使用数值中心做任意五维方向的 (2)，仍需 x1 导数。固定 epsilon 的 H_epsilon 也不是单独关于 x 的 -4 次齐次函数，不能擅用 Euler 齐次公式补出遗漏导数。

## 4. 中心 H 与四个梯度的算法源码读审

本节接受已独审的 017 锥表作为输入，未重新三角剖分，未重新生成 3504 锥 stencil。043 生产源码的关键数学步骤如下，均通过读审。

**归一化。** 锥射线 r=z/144，determinant 为 d_C/144^5，五个 resolvent 分母各产生一个 144，恰好抵消。中心 x0=n/D，D=104729，物理 face 位移 v 给 k_j=n dot z_j、l_j=D z_{j,2:5}。时间换元 s=144D t 因而给

\[
 H_\epsilon(x_0+(0,v))=
 \frac{D^4}{144(2\pi)^5}\operatorname{Re}\sum_C d_C
 \int_e^\infty\prod_{j=1}^5(s-i k_j-i l_j\cdot v)^{-1}ds,
 \quad e=942561/4.
\]

因 epsilon>0，每个实空间位置上的分母模至少为 s；有限锥和的积分及其一阶空间导数都可用 s^{-5}、常数倍 s^{-6} 控制。这里不把有符号 Fourier 锥项当正项；正性证明来自前面的 Cauchy 卷积。

**时间积分的符号。** 中心节点互异时，部分分式系数为
\(\prod_{l\ne j}(i k_j-i k_l)^{-1}=\prod_{l\ne j}(k_j-k_l)^{-1}\)，因为 i^4=1。无穷端的共同 log R 项因四阶差商消去，留下
\(-\mathrm{DD}_4[\frac12\log(e^2+k^2)]\)。所以源码累加 `-det*value` 的负号正确。平滑 f_e 允许合流极限，但普通 DD 在重合节点不允许直接执行；源码只在已检查全部 280 投影互异的这个中心执行。

**初始一阶量。** 写 e=A/B，则共同常数 -log B 的四阶差商和空间导数都为零。每个初始节点的零阶对数系数为 1/2，物理尾导数为
\(l_jB^2k/(B^2k^2+A^2)\)。源码 `slopes=D*vertices[i][1:]` 正好给物理坐标导数，未遗漏或重复 D。

**DD 一阶递推。** 当 R=P/(d+delta dot v)，有 R0=P0/d、Rj=(Pj-delta_j R0)/d。源码逐层先算 `out[0]` 再算四个一阶分量，且左右节点差与 slopes 差同向。节点按中心 k 排序只选定局部表示，不改变对称差商；中心互异保证存在足够小的邻域用于合法求导。之后区域延伸使用全局 (2)，不要求盒内节点继续互异。

**合并与外包。** 五行分别为一个 H 值和四个 H 导数；每行符号 0 是有理常量，另外五个是对数。以 B^2 k^2+A^2 合并相同对数键（包括相反 k）合法；该键严格正，不与常数键 0 冲突。有理系数用 fmpq 累加完再评估 Arb 对数与 pi，H>0 后区间除法得到四个 log 梯度。`lower()/upper()` 转 fmpq 保存有向 dyadic 端点，没有用展示字符串中点替代原球。

本轮实际核对外置 checkpoint SHA
`1da37f20a866551ac568558f57dffac2fe320e283aca17c7a624c3ffbe3904b6`，保存状态为 complete、next_simplex=3504、五行共 704 个非零系数，fingerprint 与冻结登记/源码/输入摘要一致。未调用其中的生产器或恢复运行。

## 5. 区间梯度、16 角点与整盒

对真实四梯度 g，有已保存的区间 g_j in [g_j,L,g_j,U]。任意 v 的安全线性下界为 \(\sum_j\min(g_{j,L}v_j,g_{j,U}v_j)\)，上界取 max。不能把 g 的中点代入后声称得到全域抛物面。

对每个固定真实 g，下抛物面 g dot v-v^TQv/epsilon^2 为凹函数；若 v 是角点的凸组合，其值至少为角点值相同凸组合，因此闭盒最小值至少是角点最小值。上抛物面为凸函数，类似地最大值不超过角点最大值。随后再在每个角点给 g dot sigma 作有向包围，正是 043 源码的顺序。

043 的固定闭盒下界是 `H_L*exp(min_corner lower_log_remainder)`，上界是 `H_U*exp(max_corner upper_log_remainder)`。本轮只读其完整 dyadic 下端点，准确比较确认其同时大于 1/10000 与 11/100000；未重新调用 exp。

令 Delta_L 为保存的 log(H_L/(1/10000)) 的有向下端点。对每个保存的 b_sigma,L，安全条件为

\[
 P_\sigma(h)=\Delta_L+b_{\sigma,L}h-c_\sigma h^2/\epsilon^2\ge0.
\]

本轮用纯 Fraction 对全部 16 个保存根区间 [r_L,r_U] 逐一核对：0<r_L<=r_U、P_sigma(r_L)>=0、P_sigma(r_U)<=0。因为 Delta_L>0、c_sigma>0，二次式有唯一正根，这些符号准确证明保存根括号。随后核对最小根区间、向下 2^-24 取整及全部有理代回：

\[
 h=53571/8388608,\qquad
 \min_\sigma P_\sigma(h)>0.
\]

最小剩余量与 043 JSON 的完整分数准确相同。另外，本轮从保存的四个梯度完整端点独立重新形成各角点线性下界，得到的最小有理剩余量仍严格为正；不依赖保存的 dot-product 球和分量端点间可能不同的舍入宽度。对原半宽 1/16384 也通过该精确抛物面检查。

准确中心边距为 614/104729，故与单位正 face 相容的半宽是 min(h,margin)=614/104729；源码和保存结果完全一致。最大 h 的意义仍仅限于此固定保守抛物面测试，不是实际超水平集的最大内接立方体。

## 6. 椭球的全部梯度区间及准确余量

Q>0 给对每个真实 g、v：

\[
 g\cdot v\ge-\sqrt{g^TQ^{-1}g}\sqrt{v^TQv},\qquad
 g^TQ^{-1}g=4\sum_jg_j^2+(4/5)(\sum_jg_j)^2.
\]

源码的自然区间计算包住上述二次型，取上端点后再开方取上端点 B_U，方向正确。本轮没有重新开方：二次型关于 g 凸，梯度区间长方体上的最大值可以在其 16 个顶点取到。用完整分量有向端点准确枚举这 16 个顶点，直接验证了

\[
 B_U^2\ \ge\max_{g\ \text{in saved interval rectangle}}g^TQ^{-1}g.
\]

所以保存 B_U 控制整个梯度区间，而非只控制中点。对保存
\(\rho=51857/8388608\)，纯 Fraction 代回得到

\[
 \Delta_L-B_U\rho-\rho^2/\epsilon^2>0,
\]

且与保存的完整分数剩余量完全一致。对所有 v^TQv<=rho^2，函数 Delta_L-B_U r-r^2/epsilon^2 在 r>=0 递减，故该边界检查推得整个闭椭球的严格门槛。椭球本身不必包含于单位 face；若研究域是单位 face，应取交集。

Delta_L 的超越来源仍是冻结 043 中对精确 H_L/target 作 Arb log 并取下端点。本审计检查了该源码路径与完整保存端点，没有从展示的小数再造 Delta，也没有声称独立重算其 log 值。

## 7. 后续接口与交付

046 可对每个新中心保存 H 的正下端点、四个物理 tail log 梯度区间、epsilon、固定 N 与准确区域判据，应用本接口。中心投影互异仅是特定普通 DD 实现的条件；Hessian 定理和区域延伸本身没有这个条件。若实现遇到重合节点，须用合法合流计算或跳过该 anchor，不能把除零结果当区间证书。

部分覆盖只能报告其明确覆盖集合。要交付全域 g，仍须给出整个目标域被安全区域覆盖的准确证明，或对未覆盖部分另立有效下界。本轮通过接口，不代替该覆盖步骤。

复核命令（锁后请用新的输出目录）：

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B \
/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-I/cube-dim-5-bounds/audit_harnack047/verify_endpoints.py \
--output-dir /absolute/path/to/new-audit-output
```

最小白名单为登记、审计器、结果、实际收据、checkpoint、本报告；hash 见 manifest.json。未修改 043、主仓库或旧锁；无活动后台进程，无提交或推送。
