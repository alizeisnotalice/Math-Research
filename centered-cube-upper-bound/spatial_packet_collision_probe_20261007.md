# 空间 cell 碰撞：冻结真赢家后处理

本批按新登记，对前轮 15 个完整原子来源/74 个 receiver 完成一次后处理，没有重新最大化或调用旧 oracle。三轮保留原 source、窗口和 seeds：50 个 eligible receiver、14 个 top 分支、10 个 empty 分支。219 项 eligible cell 质量界及 192 项局部 Bernoulli 对数区间比较全部通过，无 unresolved。

这些是有限 local 公式/geometry 守卫；没有核验 6W 全空间 Lebesgue 积分预算。全空间费用需要解析径向积分，不能由这 74 个 receiver 推出。本批 τ=M/q 是逐 row scalar 诊断，不构成同一阈值下的 receiver 集合或概率实验。

## 1. 固定输入与半开 cell

原完整 μ=Σ a_jδ_yj，a_j>0、Σa_j=1；所有重合标签保留并按物理 cell 聚合。沿用保存的闭 cube 全边长 D_j、实际连续 [a,b] 真赢家 R 与 M=m_R/R^n，没有比较新尺度的响应或重新选择 R。三个阶段的 n/N/窗口、完整坐标权重、原 seeds/source hash 均直接复制到新 profiles；输入 manifest 在新 registration 中冻结。

每个实例固定 d=min(a/(8n),(b−a)/4)>0。cell index=(floor(y_i/d))_i，cell=Π_i[d k_i,d(k_i+1))，corner=d k。完整 cell 质量 w_C=μ(C) 包含全部标签；captured cell mass c_C=μ(C∩Q_R)。定义 π_C=c_C/m_R、χ=Σ_Cπ_C²。π 仅在 c_C>0 的 cells 上正，完整 source 分配不变，未捕获 cells 仍完整保存。

74 个 row 的全部 cell 权重和都精确等于 1、标签恰分配一次。负坐标 floor 保持向负无穷取整；261 个原 source 标签含负 cell 坐标，289 个标签至少位于一个 cell 坐标边界，均按固定半开规则处理，未因重合或边界拆分来源质量。

## 2. 空间扩张与 winner 支付条件

记 R+=R+2d。对每个 captured cell，存在 y∈C∩Q_R，因此全部原标签 y'∈C 满足 ‖x−y'‖∞≤R/2+d，cell lower corner 同样满足 r_C=2‖x−corner(C)‖∞≤R+2d。两者只用保存 D 与原坐标核验，无重新最大化。

eligible 定义 R≤b−2d，此时 R+∈[a,b]。赢家保证 m_{R+}≤M(R+)^n，而 m_R=MR^n，因此

    m_R ≥ (R/(R+2d))^n m_{R+} ≥ (R/(R+2d))^n w_C.

所有残差为 Fraction。top 分支 R+>b 时没有把赢家上界外推到窗口外；仍保存 expanded mass 和这些残差的实际值，但将其标为不适用，未列入 eligible 检查通过数。本批 top 分支实际残差未负，不构成窗口外一般推论。

q 固定取 2、4、16，每个非空 row τ=M/q。m_R−τR^n=(1−1/q)m_R>0 的 192 项 strict band 下界比较全部精确通过。没有声称该阈值在不同 receiver 间相同，也没有测试额外上 band。

本批共 240 个 captured-cell incidents（同物理 cell 在不同 row 可重复）。完整 cell 标签在扩张 cube 内以及 corner bound 的 240+240 项全部通过；其中 219 项有 eligible winner 质量界资格。一个 corner 位于原 Q_R 外，但在扩张 cube 内。重要局限：当前固定来源中每个被捕获 cell 的全部标签已被原 Q_R 捕获，partial-capture cell 数为 0；因此这批没有压力到“同 cell 未捕获标签通过扩张新进入”的困难情形。没有为制造该现象调整源或重跑。

## 3. 同一 cell Bernoulli 与真实 collision

对每个完整 cell，令 η_C 为独立公平 Bernoulli，且同 cell 的全部原标签共用该 η_C。ν=2Ση_C μ_C 是完整源的随机正重权。可在整个实例同时定义这些 coins；本批每 row 只计算其精确边缘期望，不声称独立重采样了一个更好的源。

在原 Q_R，m_ν/m_R=2Σπη，V=2Σπη−1。完整未捕获 cells 的 η 对这一查询无影响，因此在边缘 DP 中积分掉；其原质量和标签仍保存。精确 χ=Σπ²，E V=0，E V²=χ，均在保存 DP 上独立重构核验（64 个非空 row 全部通过）。这是 cell 相关、cell 间独立的原份额 collision，没有按源标签人为制造独立 coins。

π_C 为有理数。用全部 π 的共同整数分母 T，将每 cell 写成整数 share u_C/T；subset DP 逐 cell 添加 accept/reject 两个状态，对每个整数 s 保存产生该 subset 总份额的组合数 N_s。其概率 N_s/2^k、V=2s/T−1，故准确积分全部 2^k coin outcomes，而非 Monte Carlo。最大 k=31（2^31 种 outcomes），最大 DP 分母 158、状态数 159，压缩没有近似。所有 DP 质量和为一、subset 均值 1/2、互补对称均精确通过。

## 4. 有理 log 区间与结果

令 z=log q，约定 log_+(0)=0。核验

    E[z−log_+(q(1+V))] ≤ (2+4z)χ.

对 x>1 作 x=2^k t、1≤t<2 的 range reduction。a=(t−1)/(t+1)∈[0,1/3)，使用

    log t = 2 Σ_{j=0}^{K−1} a^(2j+1)/(2j+1) + remainder,
    0≤remainder≤2a^(2K+1)/[(2K+1)(1−a²)], K=100.

log 2 用 a=1/3 的相同严格余项。全部求和和余项均为 Fraction，再向外取整到 192-bit dyadic 格。对 q(1+V)≤1 不取 log，log_+ 精确为零；所有零接受 subset 保留。

使用 expectation 上界（z upper−expected log lower）与 RHS 下界比较；两者都由 Fraction 构造，没有用 Decimal 窄区间冒充证明区间。192 项最终 gap 下界全部严格正，最小 gap 下界约 0.167443280179（精确有理值保存在 saved_dp_audit）。这认证的是本批冻结有限输入上的局部不等式，不是一般 lemma 的独立证明，更不是整体积分支付。

|轮|receiver|eligible/top/empty|captured cells|eligible cell 质量检查|local q 区间通过|χ 范围|最大 cells/DP states|
|---|---:|---|---:|---:|---:|---|---|
|1|18|9/3/6|16|10|36/36|[179/441,1]|3/8|
|2|24|17/3/4|128|121|60/60|[491/12482,1]|31/159|
|3|32|24/8/0|96|88|96/96|[521/5929,1]|14/78|

主后处理正常 exit 0，耗时约 0.967 秒。没有新增输入/阈值搜索，没有 rerun 原 posterior-age 或任何旧 MC。只读 audit 核验 15 新 profile hashes、74 完整 cell 权重和、64 个 π 和/EV/EV²，以及全部 192 local 区间方向。

## 5. 保存和复现

15 个新 gzip profiles 含完整原 case（source、weights、seed、source hash）、原 input profile hash、x/D/M/R、d、全部半开 cell 标签与完整质量、captured labels/π/χ、corner/r、expanded query mass、eligible/top flags、geometry residual、每 row/q τ、整数 subset DP、log/expectation/RHS/gap 的向外有理区间。empty、top 分支未隐藏或重归一化。

新 script 有 register/run 两个模式，audit.py 仅读取本批终态保存记录。按登记执行的文件及所有输入/输出 profile hashes 在 final_receipt；复现应使用一个空的新目录，以免 exclusive-create 输出拒绝覆盖冻结数据。任何后续 partial-cell 构造或全空间积分实验需要另登记，不属于本批。
