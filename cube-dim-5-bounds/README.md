# 五维中心方体弱型常数

固定Lebesgue、中心轴平行、共同半径、一般L1、全部正尺度。

有效界为 C_cert≤C5≤32，其中 C_cert 是约1.9943164615684859的完整代数构造，严格大于19943/10000。准确表达式、来源和公开/内部比较见 results/bounds_ledger.json；不从展示小数重造证书。接续门槛冻结在 round005/continuation_reference.json。

第013轮：同一alpha1模型的全域g=21/200000与完整I≤54702517/10^12已主审及独审；配对保证成本仍约38.47983413，尚无新的C5界。固定模型成本区间跨32，不能据上包失败排除模型。后续候选与唯一计算入口以STATE最新条为准。上界优先，下界低资源；程序运行不预设wall上限，按有限工作量保存可恢复检查点；总研究算力不预设上限。

阅读：STATE.md → results/bounds_ledger.json → rounds/013.md → 各证书及证明。当前完整100998叶可携入口为 upper_portable096/check_compact_packet_v3.py，数据在 upper_portable096/run098_v3/，不依赖外置大state。新194点由 audit_mixed100/ 不同整数路径独审，主证明与接纳见 round013/primary_globalg105_review.md 和 primary_globalg105_acceptance.json。完整I账本可携入口为 upper_adaptive095/verify_portable.py；它核几何、凸弦及全空间预算，H球仍依赖明确锁定的017数值公式和保存来源。

只有新下界严格高于L_ref且新上界严格低于U_ref并通过资格核查，才记录双侧时点；此后再执行三轮完整同步迭代，再终审及编译中文论文PDF。原cube-finite-dim-bounds及原资料包保持只读，暂停状态未改；内部结果不称公开发表成果。
