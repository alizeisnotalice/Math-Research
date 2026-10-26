**123 Sol：6237个新 EPS128 jet 与局部闭覆盖独立审计（2026-10-26）**

两批有限审计全部通过：6237个新H与24948个物理导数，共31185个分量，由直接EPS128纯整数有理部分分式路径严格内含生产球；6237个fresh Delta、99792个闭角点、18713个日志事件及最终覆盖树重建均通过。新增572个通过叶、5665个几何未证父盒的闭二分。最终98579 passed、57549 splits、11330 pending，全部52360根保持精确闭Kraft覆盖。**11330个待证叶仍在，完整g110尚未成立，全局S=0，正式C_cert≤C5≤32不变。**

本审计仅处理原123两批4096+2141个冻结原前沿任务，不消费新128 anchor-dependent曲率；所有角点罚项仍为1024。没有重跑旧H或已保存成功前缀，没有派发任何本批新geometry child，没有修改生产、共享原源或旧审计。未提交或推送。只写独立clone的audit_g123_sol/及外置scratch/audit_g123_sol；旧Astra静态材料只读，没有新调用或续派Astra；原三、四维项目保持只读暂停，没有调用共享MATLAB。

主全文审阅准备代码后，授权回执SHA为 `1079d0769b231f40b71808433c0ad154f265444ea321ec3233be96d2e5fac0b4`。准备manifest SHA为 `065ac9038dbb5b207877b543f00cc35873593ebc4268b8171f5201aae3ffb24c`；source lock SHA为 `7a580c2d68cd55e2ac03faa7cbab3de2325a4d6f2d8cadf0cd1dbe2d546e57f6`。授权之后九份已审数值、绑定、准备、故障及结构代码保持字节不变。最终锁检查确认全部4234条冻结输入锁不变。

|有限批|生产manifest SHA|审计登记SHA|新点/分量|新pass/split|
|---|---|---|---|---|
|1|fa43cf5be513eca42d452eaae8062d5331994a18f9acf589bbf2f0033a619a7f|293de5bb03a69f307334ca95fc3cf4d9c2ef2a42f344ae60928f02cf2b7d3b07|4096 / 20480|307 / 3789|
|2|cc4f1f8cd6c186cff073f289e8725037fd56f97b067c27558648aa191c032ff8|b9e0d9e18b2c8d70808e51ddd25ed5691b61aa32250e70696913a2f2b5c0d8df|2141 / 10705|265 / 1876|

旧source的消费边界见准备报告与preflight.json：seed的98007个passed逐对象等于105的96025、107的850、118的1132；51884 splits逐对象等于098的48638加107的3246；184796个H64与4290个H128完整raw cache逐对象匹配旧源。accepted mapper的leaf SHA、record SHA、epsilon与完整五维point逐条绑定。116 manifest `111840a3858a2a8b616d6fe731c16a28f8c04a18515dc9ac8fa39e879a484f2d`、122 manifest `a20ee932eedb8d83e2644bf7fdf1e74ec0aa8f93d310f25ab8ee7e48368d5e40`为冻结依赖；本轮没有重新求这些旧值、旧Delta或旧corner资格，也不转移主对其全文和总体证明的最终责任。

全部6237个新数学对象使用epsilon=1/128、target=11/100000、同一五维alpha1/16个未单位化body方向/Fourier(2π)^−5模型。typed key由epsilon和完整有理五维point组成；所有新key彼此互异，与旧4290个H128 key、095/103冻结包含EPS128的scalar/planned点集无交。故6237个任务全部为新fulljet，零cache hit；逻辑typed cache总量189086+6237=195323，不能将同点H64对象混作H128命中。840个同点H64 guard只在新旧完整point相等时使用，变锚点guard为null。

独立数值只用Python整数与Fraction，参考已审081的部分分式及舍入设计，不导入生产DD递推、Arb、flint、mpmath或作者求值函数。对每个锥的不同整数投影k_j，直接求有理pole权重及其物理导数。令d_j=∏_{m≠j}(k_j−k_m)、w_j=A²+B²k_j²、s_j=Dv_{j,tail}，其中A/B=144D/128；H项为−det·log(w_j)/(2d_j)。导数为

\[
\frac{\det\,T_j\log w_j}{2d_j^2}
-\frac{\det\,B^2s_jk_j}{d_jw_j},\qquad
T_j=\sum_{m\ne j}(s_j-s_m)\frac{d_j}{k_j-k_m}.
\]

全部3504个锥实际求和，归一化D⁴/(4608π⁵)。log采用768位向外整数atanh级数及显式余项，π采用Machin公式和216项交错级数包围；结果保存为512位向外dyadic球，并在保存后再次检查各球严格内含生产球。直接EPS128计算不调用116的EPS64(2x)缩放数值。单pole物理导数公式另经实际SymPy符号恒等式验证，残差0。

共享依赖明确为017真实fan与其模型/归一化资格、081 JSON reader及算术公式设计、075/099全空间log-Harnack理论，以及已锁旧source值资格。不同的数值路径提供独立计算证据；本轮没有独立重证几何剖分、模型或所有旧source，不把这种共同依赖说成完全独立证明。

每个实际raw checkpoint均核complete、partial=None，raw/returned/point_manifest/record逐字典绑定，同epsilon/target/model/kernel、精度、attempts、calls、PID与原任务一致，五个生产球均满足已审宽度条件，H下端严格高于目标。没有H128 obstruction或数值诊断退出。840个有效guard逐条核同point跨epsilon关系；guard未代替新H或梯度。

每批先只派一个pilot，全部五分量严格内含并fsync后才派剩余块，每块最多128点。batch1共33块（含pilot）、batch2共18块（含pilot）；每点持久保存后coordinator才ACK。全部正常完成，无fault、无重启或成功前缀重算。无wall cap；实际worker全部nice5。准备阶段坏pilot故障探针另实际验证idle/ACK worker可释放，详见fault_probe.json；该探针无H调用，不计入6237新点。

|实际入口|master / worker PID|秒|退出|
|---|---|---|---|
|batch1 numeric|36981 / 36990–36993|315.550228625|master和全部worker exit0|
|batch1 structure|38350|5.651627250|exit0|
|batch2 numeric|38809 / 38818–38821|178.467414208|master和全部worker exit0|
|batch2 structure|39483|4.764212583|exit0|

structure门只在相应新值全部严格内含后运行。fresh Delta用独立整数log核Delta≤log(HL/target)；同epsilon四个物理导数/H区间用全部端点交叉乘法核对。真实闭盒由root坐标和每四位交织path位直接生成，Delta、梯度区间、box−anchor位移统一到整数分母，穷举同一个anchor的全部16个闭角点。075理论允许盒外anchor和负Delta，只要求HL>0；本批新锚点恰在盒内且HL>target。没有按角点switch锚点，几何未证没有当数值精度故障处理。

batch1的4096次dispatch/received/merge与finish共12289个事件，batch2的2141次三类事件与finish共6424个事件，canonical hash链全部独立重放；重建overlay与生产state逐字段一致。每次dispatch取冻结原任务前缀，合起来恰为6237个原seed pending集合。失败父盒的axis=depth mod4、两个闭child准确分区，全部child保留在queue，没有被本计划派发。

最终98579 passed与11330 pending组成109909个活叶，57549 splits恰为全部严格祖先。逐根dyadic整数区间准确铺满单位区间，prefix-free且所有52360根Kraft=1，最大深度16；不是仅核计数恒等式。独立结构结果SHA分别为 `f7f26e1ea07526175f0825759428808573c7034515a835c25da24d3ba323d8a0`、`0f7ea604631d652f835299d7b2d9979660bacc9e83fc50bd5fb60e5ac326ff6c`。

最终实际ps检查确认两个numeric master、八个worker与两个structure PID均已不存在；没有残留pool。完整结果见final_check.json，SHA `f046695b0398106e38dcd294d56ffee14513d1b574ec2f89e4710edce2275858`。

最小归档只保留冻结数学代码、来源/授权与两个有限登记、零H准备回执、两个结构证明回执、本报告和全部独立值证据。外置原始数值53个JSON合计57227099 bytes继续保留；归档把各分块重复的input_locks抽成每批一份，全部proof内容及原文件SHA/bytes完整保留在 `evidence/numeric_batches.json.gz`（6030635 bytes，SHA `14a644142c7bd585cad6959b73e638970ca51b527d2da308db005135cb612be9`）。实际verify_archive.py逐文件重建标准JSON字节后确认53个原SHA和bytes完全相等；没有重新计算任何H/log。该新增归档校验器不属于九份已冻结的数学程序。

复现原执行使用准备报告给定环境与两个registration；不得重跑已经成功的6237点。只读核验归档可运行 `venv/bin/python -B audit_g123_sol/verify_archive.py`，实际路径以本目录和外置venv为准。manifest.json提供精确最小白名单及文件SHA，source_lock和两个registration列出必要外部输入。主仍保留最终证明接纳责任；本轮仅给原1024曲率下的完整有限值审与局部叶资格，没有全域g110、I或C5改善结论。
