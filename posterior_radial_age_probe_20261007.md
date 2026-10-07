# Posterior 径向年龄：连续窗口精确有限守卫

本批已按独立登记完成三轮：15 个完整正原子输入、74 个固定 receiver（64 非空、10 零响应），467 个正质量候选 L。1,955 项 Fraction 质量比较全部通过；467 项 σ=1/2 的有理包围比较全部通过。年龄均值的 160 位小数诊断也全部通过。这些点不是 Lebesgue 积分样本，不提供 geom 资格、总体概率或维数阶数证据。

## 1. 原接口与实际连续赢家

完整来源 μ=Σ a_jδ_yj，所有权重严格正、总质量 W=1。重合标签未删除或单独归一化，查询质量是全部标签的总和。全边长 D_j=2‖x−y_j‖∞；闭 cube 在 D_j=L 时捕获该标签。每个实例共用窗口 [a,b]，M=sup_{a≤R≤b}m_R/R^n。

程序使用 Fraction 坐标和权重，只在 a、b 与 (a,b] 中全部不同 D_j 上比较响应。同距离标签先全部加入质量，再比较；两个 arrivals 之间 m_R 不变而 R^n 增大，因此这是真连续窗口赢家而非半径网格。多个响应精确相等时选择最小 R。本批实际跨尺度赢家 tie 数为 0；但有 130 个多标签同距离组，92 个 D_j=a 和 88 个 D_j=b 的闭边到达被保留。每个输入都含正质量重合标签。

对任意候选 L∈[a,b] 且 m_L>0，令 π_L(j)=a_j1_{D_j≤L}/m_L，c=exp(g)=M L^n/m_L≥1，以及 A_j=n log[L/max(a,D_j)]（仅 posterior 正质量标签）。年龄支持在 [0,n log(L/a)]。

本批主核验不依赖浮点对数：对 s=a、s=L 及 [a,L] 内每个真实 arrival，准确比较

    m_s ≤ M s^n,
    m_s/m_L ≤ c (s/L)^n.

两个式子由 Fraction 分别核对。a 下的全部源质量已并入 a 初值；同距离组 simultaneous 更新。因此对 0≤t≤n log(L/a)，设 s=L exp(−t/n)，有 Pπ(A≥t)=m_s/m_L≤exp(g−t)；更大 t 的尾为零，t=0 的尾为一。实 s 位于两个 arrivals 间时质量不变、右端单调增加，离散登记不漏掉连续参数。由尾积分得到 EA≤g+1 和 Eexp(σA)≤exp(σg)/(1−σ)，σ∈(0,1)。这是接口推导；本批点核验不能替代任意输入的证明。

## 2. 登记与输入

所有数组、weights、source hash、receiver 坐标在响应 oracle 执行前已写入 registration。三个阶段无自适应输入调整。第 3 轮在同一 n 独立变化 N，且 n=64,N=32 与 n=256,N=32 各保持同一完整 μ 改变窗口。R3 同 source 的 hash 相同；receiver 构造随窗口改变，未声称 receiver 独立比较。

|轮|来源构造|n|N|窗口|固定 seed|非空 receiver|正质量 L|精确比较数|
|---|---|---:|---:|---|---:|---:|---:|---:|
|1|small_random|1|4|[1,2]|852001|4|10|20|
|1|small_random|2|7|[1,2]|852002|4|12|26|
|1|small_random|3|10|[1,2]|852003|4|17|37|
|2|correlated_ray|4|12|[1,2]|8530412|5|42|151|
|2|nested_sign_layers|8|16|[1,2]|8530816|5|33|119|
|2|separated_clusters|8|24|[1,2]|8530824|5|22|64|
|2|finite_lattice|4|32|[1,2]|8530432|5|22|60|
|3|large_block_layers|16|8|[1,2]|8541608|4|26|86|
|3|large_block_layers|16|32|[1,2]|8541632|4|55|311|
|3|large_block_layers|64|8|[1,2]|8546408|4|23|71|
|3|large_block_layers|64|32|[1,2]|8546432|4|41|192|
|3|large_block_layers|64|32|[1/2,3/2]|8546432|4|43|201|
|3|large_block_layers|256|8|[1,2]|85425608|4|18|45|
|3|large_block_layers|256|32|[1,2]|85425632|4|41|184|
|3|large_block_layers|256|32|[1,3]|85425632|4|62|388|

来源构造只是明确有限输入压力：随机 dyadic 点；强相关 ray；符号嵌套层；分离簇；有限格点；高维两块相关中心加独立符号微偏移。没有把这些名称当作已有一般下界模型的资格。具体坐标和每个权重全部保存，N 是原标签数，physical distinct count 在 results 中另列。receiver seed 为上表 seed+500000。

固定停止规则是完成全部 15 实例/74 receiver，实际一次执行正常 exit 0，耗时约 0.315 秒；没有旧 oracle、旧 MC 或旧守卫调用，也没有额外搜索。

## 3. 数值范围与证据边界

|轮|receiver 非空/总数|正质量 L|Fraction 比较通过|σ=.5 包围通过|EA 数值诊断通过|最大 EA（约）|真实赢家最大 EA（约）|
|---|---|---:|---:|---:|---:|---:|---:|
|1|12/18|39|83|39|39|2.079442|0.095894|
|2|20/24|119|394|119|119|5.545177|0.389566|
|3|32/32|309|1478|309|309|228.830449|0.211065|

实际赢家的 g=0 在 64 个非空 receiver 上精确成立，故其 EA≤1；非赢家候选可以有很大的 EA，但相应 g 同样增长，不能把 EA≤1 强套到任意 L。本批全部候选的 g+1−EA 为正小数诊断。EA、g 和 nearwinner logarithmic gap 使用 Decimal precision=160，没有对 log 实施外包围，故这部分不是区间证书。

σ=.5 时 exp(A/2)=sqrt((L/max(a,D_j))^n)，bound=2sqrt(c)。对每个正有理 z，以 floor(sqrt(z)·2^128)/2^128 为下界、至多增一格为上界，并精确核验 lo²≤z≤hi²。posterior 加权 moment 的上界均≤bound 的下界；467 项无 unresolved。这认证的是这些有理有限实例的矩比较，不认证一般定理或任意 σ。

nearwinner 不使用容差改定义。第二高节点响应（保留所有其他 R）与 M 的相对 gap 是精确 Fraction；三轮最小 gap 分别为 23/135、41/16848，以及 45665213578082641706564353/1104427674243920646305299201。无数值 nearwinner 判定参与真实赢家选择。全 gap、M、R 与响应表在 profiles 中。

1,955 是按每个 L 的实际 s 比较累计数，包含跨 L 重复的相同 cap；64 个非空 receiver 的端点/真实 arrival 上有 286 个不同全局 cap 记录，不将重复比较当独立样本。零响应 receiver 不定义 posterior、g 或年龄矩，显式记录而不丢弃。

## 4. 保留的 finite-grid 假帽反例

单源 μ=δ0，n=1、x=3/4、[a,b]=[1,2]。真实 arrival 为 3/2，实际连续最大 M=2/3、赢家 R=3/2。只用粗网格 J={1,2} 得 M_J=1/2；遗漏 arrival 上 m_{3/2}=1>M_J·3/2=3/4，精确违规量 1/4。

取 L=2、t=log(4/3)，posterior 尾为 1，而使用假帽得到 exp(g_J−t)=3/4。真实 M 的 cap 在 arrival 上等号通过。这只是“有限网格最大不能替代连续窗口最大”的反例；不是有限 J 自身 maximal 的反例，不是 weak/geom 反例。

## 5. 保存与后处理

15 个 gzip profiles 每个含完整输入 μ、每个 receiver x、所有标签 D_j、全部 endpoint/arrival 的同时捕获 mass/response/cap residual、实际 R/M，以及每个 L 的完整 posterior（包含未捕获标签的零权）、g 与年龄参考 max(a,D_j)。posterior 不重整或删除完整 source。也保存 L≥R 标记，便于后续独立登记的 periodic-coordinate/global-test 精确后处理；本批没有执行该工作。

只读重构检查全部 15 profile/source hashes、467 posterior 归一化、10,838 标签 capture 关系、286 cap residual 和 467 包围方向，全部通过；未新调用 oracle。

复现命令使用 bundled Python：

    python posterior_radial_age_probe_20261007.py register
    python posterior_radial_age_probe_20261007.py run

输出采用 exclusive-create；若要从头独立复现，应复制脚本到空的新目录，不能覆盖本次冻结登记与终态数据。registration 已保存完整数组，不需要凭 RNG 重建当前输入；numpy 仅用于首次登记，实际 run 核验为标准库 Fraction/Decimal。

主输入/执行/结果 hash 及每个 profile hash 在 final_receipt 中。不存在概率置信区间，因为这里是有限输入精确守卫，未对随机 receiver 的 Lebesgue 分布作推断。一般 posterior-age 接口的可用性不等于原 geom 门或 source-once 预算已支付。
