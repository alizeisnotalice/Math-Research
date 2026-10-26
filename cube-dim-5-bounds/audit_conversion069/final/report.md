# audit_conversion069/final：066 全17块及最终物化独立审计

**冻结的全17块与最终 full state 通过本轮独立审计。** 确认34038个测试源父，其中15143个新父严格通过、18895个未证；最终 passed94644/splits56048/queue13764。终点没有未测试的 eligible 旧源父，但 queue 仍非空，所以这不是全域 g 或 C5 结论。主任务保留终审，本轮没有派发新核。

## 实际程序、总wall限制与既有首块结果

本轮另登记在 `audit_conversion069/final/`，首块069与旧062各四件运行前后均保持原锁。实际 PID **63038**，exit 0，运行后 OS 核验 PID 不存在。单进程单线程，总wall **22.037174833秒**，峰值RSS **1268908032 bytes**。

遵守更新后的**每程序单次总wall≤180秒**限制，含输入加载、审计与落盘：正常在160秒前保存可恢复前缀，174秒有信号兜底；每块另限2048父/120秒。本次没有触及任一限额，后16块证书循环最大 **0.616972792秒/块**。进度记录每128个测试及每块完成后原子保存，恢复须显式给进度JSON及其SHA，并核代码/登记/输入锁一致。

按新指令，首块2048个81点结果**不重算**，由已锁首块069结果 `31a9967b74785c61ddadda63c6975259dab515e2642066dbe711e7233c5b4a8b` 复用，逐行核node/源索引/准确minimum/判定/81点计数与当前冻结首块证书相同；只重放其结构delta以连接后块。本次实际新算 **31990行、2591190点**；连同首块既有收据共核 **34038行、2757078点**。

新H、gradient、kernel、Green、log、exp、Gamma、Arb、MATLAB调用均为**0**。没有导入或调用作者模块、box、partition、certificate、reconstruct函数；没有执行作者重建命令。只使用Python标准库与本审计者自己的整数/JSON/dyadic实现。

## 不可变输入锁

- 唯一source061：`/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover061/batch0001/state.json`，SHA `d0d62d99d57b3f15d3acb1f4a635b188e196b0de6e155b461460ccbd64a72d31`。
- 最终state：`/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover066/final/state.json`，SHA `42c67156078cb0eb134910a5248821eb3c826c96981afb3d1d5b2fe3242a871b`。
- 066冻结manifest：`75d5754fde6813285fe59d36b68680fd31675c375ae46cfcbb261f2c6649e8a6`，包含17块全部registration/certificates/retired/checkpoint/summary文件锁。
- 最后checkpoint：`791f4a18c5e6506b724a4580ec52580cae9ec3b9712349e4b522bd06e1f463db`。
- run registration：`b8cc3a17538649edaf9472175f863c5335b51c44859bf98e89f91c1cd25b1c34`；作者代码：`e4061fa845d60f8606119c79d6a997924fd7907cc09a3daabb58443714990b46`。

先核文件字节SHA再解析，运行结束再次核所有输入及保护文件锁不变。完整绝对路径和所有逐文件SHA均在result.json的inputs中，未复制full state、sourcecache或全部块账本。

## 选择、跳过、退休与逻辑链

从source061活queue的每条node真前缀独立算初始后代数，用 `(path长度,−初始queue后代数,root,path)` 排序源eligible split。每块由当前活source索引和tested集合重新构造准确候选前缀，并额外确认候选仍有当前剩余queue后代。与注册的源索引/node/初始queue后代数逐条相等。

按注册前缀顺序逐个重放：仅当该候选是同块较早通过父的后代时才允许跳过；其余必须准确对应一条冻结测试行。所有tested源索引全链唯一，测试时仍为活split；全部accepted父在全链前缀无关，不同块也不重叠。

对当块每个accepted父，直接从当前活source数组索引计算严格后代叶及自身/后代split，退休账本每条 `source_index/node/unique_new_parent` 必须完全一致。新父真祖先继续活跃；其他逻辑记录继续以原source索引引用。各块前沿独立转换成每root共同深度D的dyadic整数区间，排序后要求从0准确连续接到2^D；全部 **52360 root** 在源及17个新前沿均prefix-free且Kraft和准确为1，活split集合同时准确等于叶path的全部真前缀。

每块独立重算逻辑对象 `{source_sha256,active_source_indices,accepted_certificates,tested_source_split_indices}` 的规范哈希，核parent逻辑前沿、parent checkpoint文件SHA、所有产物链接及post逻辑前沿。最终逻辑前沿SHA为
`c303bfd053ff1a2fae9cafbc05669e121b63cde87c93ddabac747fe12d81ab5d`。

|块|测试|通过|未证|退休queue|剩余queue|
|---|---:|---:|---:|---:|---:|
|1|2048|428|1620|3424|58054|
|2|2048|228|1820|1824|56230|
|3|2048|922|1126|5354|50876|
|4|2048|584|1464|2336|48540|
|5|2048|908|1140|2508|46032|
|6|2048|1252|796|5008|41024|
|7|2048|833|1215|3332|37692|
|8|2048|1065|983|4260|33432|
|9|2048|349|1699|1396|32036|
|10|2048|1233|815|3590|28446|
|11|2048|359|1689|718|27728|
|12|2048|410|1638|820|26908|
|13|2048|1863|185|3726|23182|
|14|2048|1358|690|2716|20466|
|15|2048|1663|385|3326|17140|
|16|2048|1110|938|2220|14920|
|17|1270|578|692|1156|13764|
|合计|34038|15143|18895|47714|13764|

累计退休source passed **3043**、splits **35614**、queue **47714**；所以最终
`passed=82544−3043+15143=94644`，`splits=91662−35614=56048`，`queue=61478−47714=13764`。
最后前沿叶数 **108408=52360+56048**。剩余source eligible父均已测试或已退休；未证记录不是H反例。

## 81点准确有理认证

固定epsilon=1/64、target=1/10000、alpha1和x1=1。对每个新测试源父，物理盒由root坐标和path每隔四位的二进制子串独立构造，准确比较保存box/anchor。旧failed_certificate完整SHA、源索引、Delta、四个**物理**log-gradient端点、anchor均与source061对应记录一致；H区间严格为正，导数/H端点交叉不等式准确核，accepted行另要求H下端大于target。

对`{lo−anchor,0,hi−anchor}⁴`逐点清分母后用纯整数计算

\[
\Delta_L+\sum_j\min(g_{j,L}v_j,g_{j,U}v_j)
 -4096\left(\sum_jv_j^2/4-(\sum_j|v_j|)^2/36\right).
\]

每条保存minimum和worst_offset均准确一致，保存passed布尔值与准确最小余量的符号一致。全部accepted父均严格通过；全链accepted最小准确余量为

\[
\frac{17501515405697254342509087167603042911706741169832190893155570471227181738815845}
{17484020030031561638275157161411670159781661551469255691266311687650138566846876483584}
>\frac1{1000000}>0.
\]

首块margin来自其冻结独立收据，后16块由本进程新算。闭象限凹性及81点控制整闭盒的分析定理沿用056/060，不重新证明；H/物理gradient球和Delta directed-log资格沿冻结源链继承，不重新生产超越函数区间。18895条未证只表示此充分证书未通过。

## 最终full state逐对象物化检查

完成逻辑链之后，另同时流式读取source061与最终full state。直接比较解析对象，**不只比较计数或整表哈希**：

- **178722**条cache的key及完整记录逐项相等。
- **79501**条保留旧passed、**56048**条保留splits、**13764**条保留queue的所有字段及相对顺序逐项相等。
- **15143**条新增passed与链累计accepted证书逐项相等。
- 除passed/splits/queue外，源其他**42个旧顶层值**全部相等；只新增`task066/status066/fingerprint066/conversion066`四键。

缓存规范SHA独立为 `f3b174e2bab3fcaed2c291090b378ddac7414393a6fb7a7bceffedca0ffaccac`，与每块cache_ref及最终conversion均一致。source中的所有旧状态、资格、计数及前轮061字段保持。

最终conversion源路径/源SHA/源fingerprint、run登记SHA、最后块号/checkpoint、最终逻辑前沿、tested/accepted、cache SHA和false新核标志全部与独立重放对象逐键一致；另独立重算 `fingerprint066=digest({registration,last_checkpoint})`。最终run receipt的计数、源/新state hash、chain终点、zero-kernel和无全域结论标志也一致。

## 精简文件、锁与复现

本final阶段仅四件：registration.json、check.py、result.json（逐块准确摘要及全部输入锁）、本报告。全部34038行完整输入证书已在不可变066账本，结果只保存逐块独立结果摘要和规范结果SHA，避免再复制近百MB数据。源和首069/旧062四件均保持原锁；无提交或推送。

|文件|SHA-256|
|---|---|
|registration.json|`22cba37ce20e15a3d1985ee281c3d3b0feb6e0de95bd97facaa25939314b8114`|
|check.py|`e1db3d4302afa71128d8d445c2f17f77864a4cef6158af742f50b2dcae4e53ca`|
|result.json|`cda5506b45ed89976547325b885efca0243aa01bd1f65c682b963da15fd4ca6e`|

```sh
cd /Users/zhengzhihao/Desktop/T/Math-Research-cube5-lower-audit-cont
python3 -B cube-dim-5-bounds/audit_conversion069/final/check.py --run
```

若触及总wall预算，结果状态只能是`verified_prefix_only`。恢复用`--run --resume-json /absolute/prefix.json --resume-sha256 <its-frozen-SHA>`，不会把未完成前缀写成全链通过。主任务仍负责最终接纳，final queue13764尚待以后授权研究处理。
