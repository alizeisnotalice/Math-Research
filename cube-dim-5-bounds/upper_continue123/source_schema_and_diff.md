# upper_continue123：只读种子、完整6237原前沿与dispatcher差异

**当前0新H/gradient/log，尚未启动生产pool。** 主授权的123来源/8worker/两批dispatcher变更现已具体化，待主快速审后才派发。所有旧118/107/105/098文件保持字节锁。

## 真实只读source schema

外置 `seed.json` 为761623221bytes，chmod444，SHA `52fe7f9f4ecd8d4c16c1c4d48423f78bc09821770fb8e3c5d8899a9dfa4f691b`。它是完整状态，不盲套旧source105：

|字段|内容|
|---|---|
|`passed`|98007个原对象：105的96025+107的850+118的1132，逐字段/哈希保持；包含原194个EPS128已110通过叶和118盒外锚点叶|
|`splits`|51884个原对象：098的48638+107新增3246，无118新增split|
|`roots`|52360个ordered root，原根列表不变|
|`queue`|准确118的6237个原对象顺序：877原未测+5360剩child|
|`cache_by_epsilon['1/64']`|原184796个完整raw对象，含历史固定H64障碍|
|`cache_by_epsilon['1/128']`|原194+1074096=4290个完整jet，逐对象保持|
|`inflight/resume_jobs`|空|
|`source_locks`与mapper/proof refs|完整源hash、98007 leaf来源/record SHA映射、118 proof library及4290已存target110 Delta/gradient资格|
|`status/globalg110_qualified`|partial candidate/false，不能由当前source预计S|

seed活叶104244=98007pass+6237queue=52360root+51884split。准备过程核整树闭Kraft/反链；作者另逐对象核typed cache、原leaf和全部几何，避免只凭数量作覆盖结论。seed保留原118已通过的盒外anchor证书，不要求它们回到leaf内。

## 完整point/time注册与有限输入

先保留原877顺序；5360children按118真实“最佳单anchor整盒最小余量”升序、准确physical volume降序、root/path排序。该max只作优先级，不是把不同anchor角点拼成证书。`frontier_order_mapper.json`保留每个新序位置对应的原6237 queue index。

全部6237在真实闭盒内部按046有限midpoint/三种eta perturbation规则预检280个非零互异投影。与095完整13406点和103完整28427点在准确time1/128/五坐标上查交：111个103冲突、11个095冲突，共122，均选下一合法候选2。最终候选1有6109、候选2有122、候选0有6，6237互不重复；0个EPS128 fulljet hit、0个scalar/planned重复、0个已失败同region同EPS128同anchor重复。

877个旧107未执行H128 reservation逐点转接到123，不重复已运行点。全局4290个EPS128 fulljet全部可查；preflight若遇到完整hit会直接consume（当前实际hit0），而scalar-only/planned H128不能假装具备梯度。112 manifest声明 `new_H_jet_log_exp_asinh=0`，它新增的scalar mean/cell不产生任何新H point keys。

冻结两个输入包：第一4096、第二2141，总6237，全部是原pending。第一包按旧877优先后再取3219 child；第二2141也是固定原child。每个任务明确 `ordinal`、`batch`、`batch_ordinal`、seed queue index、node/physicalbox、EPS128/target110、candidate tests、typed cache key和checkpoint路径。原877只有精确同点才带H64 guard；新的child没有“原H64候选”，用−1/null，禁止混用父点guard。

## 生产core与准确dispatcher diff

`worker.py`直接逐字复制107v2，SHA `62375ec07bf0a030c397c4a925b0809ab1f625edf75a19f434e4eb7f8a492e8b`。binding/cache/certificate三个core文件也与107逐字相同；target11/100000、EPS128、3504完整锥、physical四梯度、256/384/512数值梯级、complete时partial=None均不改。底层098 kernel仍只读128d98…。

完整差异见 `dispatcher_vs107v2.diff`。dispatcher的实质变化：

1. 来源改成上述唯一只读seed，accepted98007，完整旧typedcache189086；overlay只追加新对象，不覆盖旧source。
2. worker_count由4改8，仍nice+5、no wallcap，先登记实际master/pool PID再H。
3. CLI `--batch`只允许1或2；count从固定4096改为冻结plan的4096/2141，选定节点prefix/hash逐项核对。没有无限controller。
4. 第二批必须读第一批真实CURRENT/state/journal SHA，且上一批正常complete、无inflight/fault/obstruction。第二批重置本batch local ordinal/journal，保留累计cache/通过/split/queue和parent checkpoint。
5. 新几何child始终排在未消费原输入之后，两个plan只选原6237；第二批prefix必须准确等于冻结2141原child。新child只pending，不能混进第二批派发。
6. 加明确的complete fulljet reuse路径：旧raw对象不改，只为job建立显式EPS/target/model view；0新H、不重复插入cache，target110 Delta/gradient从冻结saved proof逐字段消费。其独立16corner helper允许global075盒外anchor，仅用于g，不用于积分mean。已知同region同anchor失败在preflight排除。

每点checkpoint fresh-only，任意旧batch/point路径存在就拒绝盲重启。每journal事件fsync，state每128事件及finish/fault保存。严格H128 upper<110只固定cutoff障碍，保留inflight并停止；真实fault同样保留成功prefix/诊断。每个有限程序结束close/join并核完整根覆盖/cache恒定，代理经作者复核后才接已授权下一原prefix。

## 实际零H核验与主审边界

准备PID22890，22.955秒；独立作者source/preflight PID26295，18.745秒。作者以source逐对象SHA比对98007个leaf、189086 typedcache、51884split与完整6237 queue，核所有6237真实box/280projection/候选拒绝原因、全部point/time去重、877 reservation转接及52360根闭Kraft。只读seed保持444。

`review_checks.py` PID27201，0.125秒：核心字节相同、两有限prefix限定；synthetic legacy fulljet消费只形成job view而rawcache不变；synthetic outside-anchor单证书16corner通过但错target拒绝；模拟第一批几何split后仍优先第二原prefix、第二批重置local ordinal且保留累计cache/children；全部6237任务均通过不变v2的EPS/target/guard验证。未构造/执行真实数值H来测dispatcher。

116已锁4096新值独审，122已锁118复用独审；最终prepared包另列这两receipt SHA，不改在较早时刻物化的seed或旧source。主仍保留其接纳责任；新的123值未来也需指定独审。

当前0新H/梯度/核/log，生产文件夹batch0001/0002均不存在。准备包只是让主可审具体版本；审阅后实际入口是 `run.py --batch 1 --target 11/100000 --reviewed-manifest-sha <准备包SHA>`，第一批无fault且作者闭树/cache核验后才调用同一版本 `--batch 2`。不需另向人类确认，仍不得派任何本两批新生child。

目标是主指导的“完整g110+112同模型I”正S资格；112/121 scalar/plan和完整I图、g110全域、准确幂门须全部经主/独审后才可入S。当前partial前沿的局部通过不算全局S收益。本目录未计算/宣布I/C5结论，未提交、MATLAB或改旧source。
