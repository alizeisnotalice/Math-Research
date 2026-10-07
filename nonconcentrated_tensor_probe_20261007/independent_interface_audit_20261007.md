# 独立输入与 reduction 接口审查，2026-10-07

结论：在登记的生产参数 n=16,32,64、L=64、q=(1,2,4) 下，未发现需要停跑或重跑的输入、连续候选、到达聚合、局部坐标、微盒上界或未知分支方向错误。审查通过的范围是 **完整有限原子输入及其浮点 reduction 实现**；不构成区间认证、actualgeom 认证、连续时间误差界或一般能量定理。

本次只读 `probe.py`、`microbox_check.py`、`preregistration.json`、README、summary 及已有验证/完成收据。未执行或 import 实验脚本，未重跑实验，未重构 NPZ profiles，未改动脚本、数据或其他说明。保存 profiles 的重算由 root 独立承担。本文件是本次唯一新增文件。

## 1. 完整 shared-resolution 输入与 source 采样

原输入确为

\[
\mu=\sum_{j=0}^2p_j\left[\frac1{Lq_j}
\sum_{k=0}^{Lq_j-1}\delta_{k/q_j}\right]^{\otimes n},
\qquad (p_0,p_1,p_2)=(1/4,1/2,1/4),\quad L=64.
\]

它有总质量 W=1，三分量是完整张量网格，网格嵌套后的不同原子数为 256^n。`draws` 对每个整个 source 一次抽取分量标签，然后给定标签逐坐标均匀抽取。这正是上式的正混合；不是逐坐标重新混合分辨率。整数 quarter-unit 表示 `Z=k*(4//q)` 对应物理坐标 z=Z/4。边界来源保留，未进行内部条件化。

源和 cone 的有限随机样本只估计外层期望。每个 receiver 的 m(x,R) 使用三分量完整计数乘积，不以这些 source 样本替代 μ。`logaddexp.reduce` 加总分量质量，故同址原子的质量自然相加；同址质量没有被删掉或按分量平均替换。方向先均匀选择 2n 个面之一、其余坐标均匀于 [-1,1]，使用 x=z+rω/2，与归一化 cone 一致。

## 2. 局部坐标与连续 winner 的候选完整性

对任意目标分量 q，令 s=4/q、b=⌊Z/s⌋、e=Z mod s。目标格点指标为 b+o，它相对 source 的精确 quarter-unit 位移是

\[
\ell=\frac{s o-e}{4}=\frac{o}{q}-\frac e4.
\]

代码计算到达侧长 `2*abs(dx-local)`，其中 dx=rω/2，Z 只用于 support 检验 `0<=b+o<Lq`。没有形成大的全局坐标 z+dx 再减去格点，因而避开此前发现的消去错误。

因为 |dx_i|≤1 且 R≤2，相关目标格点必须满足 |ℓ_i|≤2。又 e/s∈[0,1)，故所有可能相关整数 o 均包含于代码的 `-2q,...,2q+1`；其多出的候选不影响计数。对于抽中的源分量，Z 可被 s 整除，o=0 对应 ℓ=0；forced face 的源原子到达值为 2|dx_i|=r，且该恒等式在这里的浮点构造中保留。

对每分量逐坐标排序到达值，以 R=1 的闭区间计数作 baseline。零计数坐标单独记录：第一次到达以零 log increment 配合 first-arrival 标志解除该坐标的零质量状态；以后由 log(k+1)−log(k) 更新。故没有非法使用 log 0，也没有把首个原子计成倍乘增量。所有分量、坐标到达值随后全局排序。

候选包含 R=1、每个在 (1,2] 的坐标到达侧长，以及 R=2。两次到达之间，完整混合的捕获质量固定，m/R^n 随 R 递减，因此这套候选覆盖连续 [1,2] 最大值，无需 R 网格。在同一个 **完全相等的浮点到达值** 上，代码仅给整个事件组更新后的末事件打分；各分量都更新后才合成 μ 质量，所以同址原子和跨坐标 simultaneous arrivals 均包含。R=2 的计数乘积另行直接重算。

`argmax` 的 first-index 规则在计算出的分数完全相等时选择较小 R。任意数学 score ties、几乎相等的 score，以及其他格点到达值的舍入，仍属于浮点限制，未被符号或区间算法认证。局部源 face 恒等式不能推广成对全部数学并列的认证。

## 3. 任意 side-h 微盒的安全 upper 与阈值方向

对任意闭轴盒 B（每边 h），给定候选 Q=Q(x,R)，令 c_ji(Q) 为分量 j 在坐标 i 被 Q 捕获的格点数。长度 h 的闭区间在间距 1/q_j 的格上最多含

\[
K_j=\lfloor hq_j\rfloor+1
\]

点。因此，逐坐标交集给出确定性的分析上界

\[
\mu(B\cap Q)\le G(Q):=
\sum_j\frac{p_j}{(Lq_j)^n}\prod_i\min\{c_{ji}(Q),K_j\}.
\tag{A}
\]

同一个 B 对不同分量未必能同时达到逐分量最大值；把其正上界相加仍安全，可能偏松。有限支持边界只会减少 c_ji，不会破坏 (A)。分母必须是同一个完整输入的 m=μ(Q)，代码确实如此。

`caplogs` 在一般维度中实际使用更松的

\[
\widetilde G(Q)=\sum_j\frac{p_j K_j^n}{(Lq_j)^n}
\mathbf1_{\{c_{ji}(Q)>0\ \forall i\}},
\qquad G\le\widetilde G.
\tag{B}
\]

所以小维验证中它不一定等于 (A)，甚至 \widetilde G/m 可大于 1；这不使上界不安全。本次生产 h=2/n、n≥16、q≤4，有 hq≤1/2<1，故全部 K_j=1。这时 (A) 和 (B) 完全相等，因为每个 \min(c,1) 只记录非空与否。

主比较是 `log(upper/m) < log(eta_n)-tol`，其中

\[
\eta_n=\frac{49}{65536\sqrt n},\qquad
\mathrm{tol}=(n+1)10^{-10}.
\]

方向正确：成功意味着任意 side-h 盒的 captured share 严格小于 η_n（条件于计算出的计数/赢家），即足够判定非集中；失败仅为 unknown，不证明集中。减少 η 使通过条件更严格，符合 η_n 的用途。固定 η0 比较仅为另存的诊断，没有换 λ，也没有替代主阈值。n=32 的平方根及所有 log 比较使用浮点，slack 是数值筛查，不是向外舍入证书。

## 4. Unknown 保留及实际研究集合

`raw_cert_res` 只保留计算 winner 上通过微盒充分条件的点。lower 在接近最大值的有效候选上取最小值，并收紧 band、delay 比较。upper 在这些候选上取最大值，放宽 band、delay 比较，**没有乘 cert**，因此未通过微盒条件的 unknown 全部保留。upper 是完整 hard relaxation 的浮点包络，不能解释为已求出的精确非集中残余能量。

原 raw hardband 是 2λ<u≤4λ，lower endpoint 严格、upper endpoint 闭合；λ 固定为 log λ=−n log64。`Theta_hard` 的权重在激活 d=n log(R/r)≥0 时为 exp(−d)，与 (r/R)^n 一致。`Gamma_v1`、`Theta_hard_v1`、`Gamma_vstar` 的额外 delay 条件与登记说明一致。源 face 的 R=r 被 `rr==r` 保留。near-score、band 和 delay 的浮点包络并不控制所有机器误差。

这里测试的是 Eexists 的补集，其中 Eexists 表示存在某个 side-h 盒捕获超过 η_n m；充分筛查得到该补集的子集。它不是某个 greedy finite paid cover U 的特定残余 E\U：有限 cover 可以把部分精确非集中点也支付掉。因此这组 screened residual 压力不能直接当作某一个 greedy cover 残余的反例。README 已明确该区别。

当前字段没有 FIRST、continuation、CP、GP 或其他实际 geom gates。它们可用作有限完整输入上的 hard relaxation 压力，不能称为 actualgeom 能量的认证；也没有建立原子输入到正宽 L1 输入的下界转移或完整 B41 下界复现。

## 5. 已有验证收据与有限数值限度

本审查读取而未重跑：`validation.json` 的 48 个 n≤4 全原子独立扫描给出最大 log-score 差 1.3322676295501878e−15、最大 winner 差 2.220446049250313e−16；100,000 次 source 标签检查记录了经验 prior。`microbox_validation.json` 的 40 个小维全原子 interval-placement 检查，以及生产维度的 local face 恒等式检查，与上述实现方向一致。它们是浮点实现验证，不是精确有理或区间证明。当前 README 已使用这两个准确误差值。

生产每轮、维度、replicate、A/B 的 seed namespace 不同；同一个 side 在全部 t 节点复用 source/cone 路径。独立 A/B profile 乘积避免直接平方有限 source 均值的偏差，估计的是对应有限节点梯形统计。四对重复样本、近似 t 区间、有限节点 Hoeffding 界以及 coarse/fine 差，都没有认证连续 t 积分、浮点界或未知残余分支。

一个诊断命名需保留口径：`boundary_sources` 统计任意坐标 Z=0 或 Z≥4L−4。对 q=2,4，它包含右端最后一个物理单位厚度内的格点，并非仅仅该分量最外端的单层格点；这只是边界层计数，不参与 source 剔除或任何 gate，因此不影响主接口。

综上，这批数据的完整 μ/source 接口、连续计数 winner 算法及微盒充分条件方向可接受。使用数值结果时仍须同时保留 finite-input、floating-screen、unknown-upper 和 hard-relaxation 四个范围限制，不由三维结果拟合一般 n 阶数。
