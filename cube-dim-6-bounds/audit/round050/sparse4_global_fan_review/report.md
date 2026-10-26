# F050 sparse4 全fan独立审计：完整数学入口 PASS

接受beta=1/80真实152列frame的完整fan：21个parent、478条primitive rays、4266个full cells、13229个非退化full simplicial cones。首parent继承已独审F046，其余20个parent的全部保存active ledger及全局拼接通过。

这是完整数学fan的裁决，**不是数值引擎放行**。没有热核/Green计算、旧tube或新S2继承、全I/g或Cn上界更新。

## 覆盖证明的关键验收

独立reader逐parent从真实frame重建crossing forms和H，避免仅核保存数据内部一致。对全部783994记录采用消去u6后的5×5整数det，区别作者6×6矩阵检查；只核保存witness，不求新解。完整lexicographic组合序列、S/O/I分类、first-appearance ID、顶点见证、压缩与未压缩hash全部吻合。

|分类|独立条数|
|---|---:|
|Singular|247176|
|非奇异且在父simplex外|319333|
|非奇异且在父simplex内（未去重）|217485|
|总计|783994|

Closedcell的全部顶点由完整pool准确过滤。所有rank4 cut facets都强制有单bit邻居，其交集恰为对应facet；保存图连通且对所有强制邻居闭合。真实全邻接图通过避开codim2交集的generic内部路径证明连通，故保存图不是遗漏部分区域的BFS子图。低维零/tie由有限closed full cells的稠密闭包保留。

全部20888个新pulling节点的facet、apex和递归自底向上核验。覆盖来自完整半空间partition和facet/pyramid归纳，不由volume或BFS旗标推断。每child真实152forms同号由对应fullcell的完整152form检查和逐child顶点包含严格推出；未重复超过预算的冗余逐child/form循环。

全局所有primitive rays的q由真实152列重求；所有新child的signed primitive/slice determinants和scale恒等式通过。全部28个内部旧facet的诱导triangulation匹配，70个外facet均有准确chamber墙见证，126个parent facets全覆盖。Global ray order保留首parent的原顺序；其局部构造未重跑。

完整证明见proof_review.md，其中独立说明BFS连通、closedcell边界覆盖、pulling归纳、正比例截面转换、hyperoctahedral复制和global负号权重。

## 归一化与计算接口边界

全局xi1=1截面准确体积为1/120；每parent的相对体积为1。两者都是证明外的冗余一致性检查，不能代替coverage。

群大小46080，global负号factor2未额外插入。未来完整群逻辑项为609592320；这仅是结构成本，不是实测kernel时间或Cn改善。正时间Fourier锥积分数学恒等式成立，但rational q/int64域、zero/repeated Fourier节点、精度/cancellation及time0积分交换仍未审。旧21cone数值引擎不得直接调用本fan。

## 实际独立工作与停止状态

PID51829，一次独立批耗时11.102822秒，源锁至数学核验完成204.307071秒。783994个det/witness记录；0新solve、0新vertex、0新cell搜索、0首parent重枚举。matrix/结构操作总计1695425，低于2000000；详细分类在execution_receipt.json与results.json。

计数口径为每ledger准确det及保存witness闭合一条；cache-miss rank、其它det、ray/form dot、cell/form sign与列明结构对象单独计数。没有把普通整数乘加伪报为独立matrix，也未以batch为单位隐去大规模form检查。

0数学失败/重试，0作者writer/checker执行、0kernel/native/MATLAB/grid/新时间或空间点/扩招/push。开始时只读查找不存在theorem.md返回缺文件，随后读取实际method.md；此准备期文件查找已登记，不涉及计算重试。所有冻结源前后hash不变。

本独立目录仅保存reader、证明、逐parent和操作收据，未复制或改写作者783994原始ledger。它们作为锁hash的必要输入保留在B050包中。未用完的时间和操作额度不授权后续数值工作。

## 准确来源

- B050 manifest: 891160dc9e98571218bcb56717b374a6c4e81f05db58c5fe340c387d3544ae51。
- global_fan.json.gz: 2405ccad9b2b0c344167539e769e3c64e3275e6e4ead208fbedc3167be47c948。
- 解压内容: 522308c63bf203f4b5d6170ada14a81f3e824d0b014fe1d50aa78bb5d22795d6。
- F046继承manifest: 6b347c5b9aa900f11893afbd0a7e706fda066c7891496d119eed7226bdd59edd。

本目录由manifest.json封存；父代理保留最终证明与集成责任。本批完整PASS后停止。
