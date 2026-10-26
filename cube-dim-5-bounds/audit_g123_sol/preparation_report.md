**123 Sol 独立审计准备（2026-10-26）**

状态为“来源与程序准备完毕，待主阅读源码后授权 actualaudit”。本阶段新增 H、梯度值、数值 log 均为0。只写独立 clone 的 `cube-dim-5-bounds/audit_g123_sol/` 和外置 `cube-finite-dim-scratch/audit_g123_sol/`，未提交、推送或改生产。旧 Astra 生成的冻结材料仅只读消费，没有调用或续派 Astra。原三、四维项目保持只读暂停；没有调用共享 MATLAB。

来源锁定生产 prepared manifest `1edbc7912626737020fec146c80247f3f5876603732d97b0d6ecdd45aa46165f`、worker `62375ec07bf0a030c397c4a925b0809ab1f625edf75a19f434e4eb7f8a492e8b`、run `0bd96824021c94334aa9c68ff46e8b81f5de83609431954b53db2086b2ca5cf9`。全部具体绝对路径与字节 SHA 在 `source_lock.json`，最终准备包 SHA 在 `preparation_manifest.json`。

**旧来源消费与零 H 检查。** 实际 Python PID34528，exit0，15.902275292秒，峰 RSS755908608 bytes。只读流式逐对象核 seed：其98007通过叶等于105的96025、107的850、118的1132逐对象拼接；51884 splits等于098的48638加107的3246；184796个 H64和4290个 H128完整 raw cache对象与各旧源完全一致。52360个有序根等于明确组合索引集合。每个 accepted mapper 的 leaf SHA、epsilon/point pair key及record SHA逐条绑定，不以数量代替来源。没有重跑105、107、116、118、122已成功的旧数值或旧 Delta。

把98007 passed与6237 pending一起重新构造104244个活叶的严格前缀祖先集；全部52360根的dyadic整数区间逐根准确铺满单位区间，Kraft=1，且全部51884 splits恰等于这些严格祖先，最大深度16。6237 pending仍全部待证，不能从闭覆盖推出它们通过。

6237个有限原前沿任务对应相同seed queue集合，前877精确保留原顺序。全5维有理point、EPS128 typed pair key、target11/100000、真实root/path闭盒及280个整数投影generic资格逐项核对。所有point严格位于各自盒内，6237新点彼此互异，与4290个旧H128 key及冻结095/103所有包含EPS128的scalar/planned点集无交。本批0旧fulljet命中，全部6237为新H128对象；旧H64同point对象有840个有效guard，变锚点不携旧点guard。完整6237原任务分4096和2141，新增geometry children只入queue，不属于本计划派发集。

固定输入延续116 manifest `111840a3858a2a8b616d6fe731c16a28f8c04a18515dc9ac8fa39e879a484f2d` 和122 manifest `a20ee932eedb8d83e2644bf7fdf1e74ec0aa8f93d310f25ab8ee7e48368d5e40`。其旧资格作为依赖消费；本准备不代替主对116/122全文的最终数学接纳，也不把作者qualification标签当新值独审。

已锁生产batch1 raw point manifest `40342fbe1d0bf596dc3a31aac741237ad926b17539c8f96a5988f1e3e331b656`、state `75aaf1ef83f62b2e72c2b524a18af9749ab1ff50f5c276a0610be2ba6834e4a9`，以及4096个实际complete checkpoint，逐point匹配原任务、完整record和partial=None。生产作者所报307新pass、3789新split尚不是本审计的结论；须等待独立新值及Delta/corner/journal门。

**独立新值程序。** `pf128.py`只使用Python整数和Fraction；参考081已审的有理部分分式和向外舍入设计，直接令a=144D/128，而非116的EPS64(2x)缩放调用。对每个锥的5个不同投影，设d_j=∏_{m≠j}(k_j−k_m)、w_j=A²+B²k_j²、s_j=Dv_{j,tail}，直接求−det·log(w_j)/(2d_j)，四个物理导数用显式权重导数与log导数求和。没有导入生产DD递推、Arb/flint或任何生产求值函数。共同依赖明确为真实017 fan、同一模型和归一化D⁴/(4608π⁵)，以及075全空间log-Harnack理论；这是不同数值算法的独立核验，不是模型或几何的全独立重证。

所有log以768位向外整数atanh级数求值：2N+1的显式几何余项，N=282；π以Machin公式和216项交错级数包围。每个H及四个物理导数的独立有理球均须严格内含生产球，512位向外dyadic保存后再次检查严格内含。实际值计算尚未运行。SymPy已实际验证单pole项的物理导数公式，残差严格为0；这里只是符号恒等式，没有新H或数值log。详见 `symbolic_probe.json`。

**实际批次方案与恢复。** 提议每原batch最多4个nice5 worker，先1个pilot的全部5分量成功，再派固定最多128点块；整个两批最多6237×5=31185新分量，每点fsync保存后coordinator才ACK。恢复核同一code/input锁，跳过每个已持久保存成功点，不重复成功前缀。无wall cap；finite count、STOP/fault、错epsilon/target、partial未清、production球不包含、源SHA变化、worker死亡、异常返回等均停止并保留前缀。旧H128目标障碍只属固定组件，不声称排除G或C5。

新coordinator显式避免旧116潜在故障挂住：fault先保存前缀、停止派发，再给每个tasks队列sentinel、每个ACK通道False；有界join故障worker后终止仍存活者，最终join登记exitcode。真实零H故障探针PID35169令pilot在raw header处失败，另外3个worker空闲，四个PID35179–35182均正常exit0，0.096344084秒释放。没有调用P.values；详见 `fault_probe.json`。

**值后结构门。** `structure.py`须看到每个新point全部5分量独审成功；再独立整数log验证fresh Delta≤log(HL/target)、四梯度/H全部端点交叉乘法，按root网格和每4位交织path直接重建物理闭盒，以共同整数分母穷举同一个anchor的全部16角点。允许075一般情形的盒外anchor和负Delta，只要求HL>0；本有限新锚点恰在盒内。绝不按角点switch不同anchor。

随后独立重放dispatch/received/merge/finish canonical hash链，重建overlay全部字段，并重构全部根closed Kraft树；保留未证及新closed children、pending/inflight边界，不能由pending存在时声称global。两批完成仍可能只有局部叶进展。当前正式C_cert≤C5≤32、完整已审g105与I112/121 cap资格不由本准备改变，全局S记0。主保留最终证明责任。

复现准备与主阅读顺序：`prepare.py`→`pf128.py`→`numeric.py`（特别release）→`bind.py/raw_gate`→`structure.py`。先读最终manifest和source lock；`bind.py`可零H绑定一个实际完整batch，再由主明确授权 `numeric.py --main-authorization actualaudit123-approved-after-source-review` 并传入实际registration SHA。没有该授权不得启动新值核。`structure.py`只在该batch的全部新值通过后运行。
