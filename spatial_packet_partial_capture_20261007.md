# 同 cell 部分捕获：新构造的精确三轮守卫

本批新登记 9 个完整正原子输入，三轮 n=1、16、256，每轮三个窗口。没有修改或重跑之前的 posterior-age/cell-collision 输入与程序。一次主执行正常 exit 0：9 个 row 均为真实 partial-capture，9 个 lower corner 均在原 winner cube 外而在 expanded cube 内；135 项赢家/几何检查、27 项真实候选 cap 和 27 项逐 row threshold 比较全部通过，empty/top/failure 数都为零。

这是补齐前批没有实际 partial cell 的局部实现压力，不是全空间 6W 积分、近优输入、weak 反例或原 geom 资格认证。来源是原子测度，不宣称本批已做 positive-width L1 迁移。

## 1. 原输入与解析构造

每个 n 使用窗口 [1,2]、[1,1+1/(16n)]、[2,3]；a,b 固定后

    d=min(a/(8n),(b−a)/4), ε=nd/(10a)>0,
    C=[0,d)^n,
    y_in=(9d/10,...,9d/10), y_out=(d/10,...,d/10),
    μ=δ_y_in+εδ_y_out, W=1+ε,
    x=(a/2+d/2,d/2,...,d/2).

完整来源没有归一化为质量一，也没有删除未捕获 out 标签。两原子确实都在固定半开 cell C 内。输入全部为显式有理数，没有随机 seed 或自适应调参；完整数组在 registration 与 results 双重保存，各 source hash 冻结。

|n|窗口|d|ε|D_in|D_out|R+2d|corner r|
|---:|---|---|---|---|---|---|---|
|1|[1,2]|1/8|1/80|9/10|11/10|5/4|9/8|
|1|[1,17/16]|1/64|1/640|79/80|81/80|33/32|65/64|
|1|[2,3]|1/4|1/80|9/5|11/5|5/2|9/4|
|16|[1,2]|1/128|1/80|159/160|161/160|65/64|129/128|
|16|[1,257/256]|1/1024|1/640|1279/1280|1281/1280|513/512|1025/1024|
|16|[2,3]|1/64|1/80|159/80|161/80|65/32|129/64|
|256|[1,2]|1/2048|1/80|2559/2560|2561/2560|1025/1024|2049/2048|
|256|[1,4097/4096]|1/16384|1/640|20479/20480|20481/20480|8193/8192|16385/16384|
|256|[2,3]|1/1024|1/80|2559/1280|2561/1280|1025/512|2049/1024|

## 2. 实际 continuous winner 与闭边 arrival

实际全边长 D_j=2‖x−y_j‖∞。第一坐标给 D_in=a−4d/5、D_out=a+4d/5；n>1 的其余坐标全距离为 4d/5。由于 d≤a/(8n)，a−4d/5≥4d/5，故这些确实是 max 距离，不是把第一坐标误当完整 l∞ 距离。

在 R=a，只有 in 被捕获，m_a=1。out 在其真实 arrival D_out 上按闭 cube 等号捕获。完整实际候选集合是 {a,D_out,b}，不存在被遗漏的窗口内到达；区间之间质量不变、响应随 R 增大严格下降。代码用全部真实 Fraction 距离生成这些候选并比较完整质量响应，而不是预设 a 为赢家。

另核 Bernoulli 解析残差

    (1+4d/(5a))^n ≥ 1+4nd/(5a) > 1+ε.

因此 (1+ε)/(a+4d/5)^n < 1/a^n，实际唯一赢家 R=a，M=1/a^n。n=1 的首个 ≥ 是等号，第二个 > 仍严格；所有维数均保留两者的 exact residual。b 的响应更小。程序最小 R tie 规则保留，实际九例 tie count 都为 1。

所有 n 的原查询 m_R=1、captured labels=[0]，而完整 cell mass w_C=1+ε、labels=[0,1]。未捕获部分正质量 ε 未在 cell weight 中扣除。这是真实 partial cell，不是同位标签拆分或 posterior 重归一化造成的假象。由于仅有一个被捕获 cell，其 π_C=1、χ=1；本批不额外重复前批的 Bernoulli log 实验。

## 3. 扩张与 corner 支付链

expanded=R+2d=a+2d，2d≤(b−a)/2，故严格留在窗口内，top=false。两标签满足 D_j≤expanded，完整 cell 标签全部进入 expanded cube。原 half-open lower corner 是 0，其真实全边长

    r=2‖x‖∞=a+d>R,
    r≤R+2d.

所以本批代表 corner 确实未被原 Q_R 捕获；仍可通过扩张预算处理，不假定 representative 就是捕获 atom。

保存并精确核验

    M(R+2d)^n−m_(R+2d) ≥0,
    m_R−[R/(R+2d)]^n w_C ≥0.

这里 m_(R+2d)=w_C=1+ε，分母、质量都是完整 μ 的真实值。27 个 actual-candidate cap residual 也全部非负，赢家上 a 的残差等于零。

阈值诊断固定 q=2,4,16、τ=M/q；27 个残差

    m_R−τR^n=1−1/q>0

全部精确通过。τ 每个 row 依赖 M，不把九点当成同阈值 Lebesgue 水平集。

|轮|n|records|partial cell|corner 在原 cube 外|赢家/几何检查|候选 cap|strict τ 检查|empty/top/failure|
|---|---:|---:|---:|---:|---:|---:|---:|---|
|1|1|3|3|3|45|9|9|0/0/0|
|2|16|3|3|3|45|9|9|0/0/0|
|3|256|3|3|3|45|9|9|0/0/0|

这些检查全为 Fraction 算术，无浮点、log 区间或 receiver 数值积分；不存在 Monte Carlo 置信区间。实际二高 response 的精确相对 gap 也保存在每 row 中，全部严格正。不要把重复检验次数当独立 receiver 证据或维数阶拟合。

## 4. 冻结与复现

同前缀 .py 是 standalone 新构造和 oracle，不导入旧脚本。registration 在响应执行前冻结全部九输入、checks、预算和脚本 hash；results 保存完整 μ/x、每坐标距离、全部真实候选/捕获标签/响应、source hashes、partial/expanded/corner/residual、每 q τ。代码输出 exclusive-create；复现要复制到空的新目录，不能覆盖冻结旧批。

最终收据包含脚本、登记、结果、执行收据和本文的 SHA256。原先 collision/age 批没有被读取、改动或重跑，新增输入不回填其已结束的压力结果。
