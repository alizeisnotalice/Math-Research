# 五维中心方体弱型常数

固定Lebesgue、中心轴平行、共同半径、一般L1、全部正尺度。

有效界为 C_cert≤C5≤32，其中 C_cert 是约1.9943164615684859的完整代数构造，严格大于19943/10000。准确表达式、来源和公开/内部比较见 results/bounds_ledger.json；不从展示小数重造证书。接续门槛冻结在 round005/continuation_reference.json。

第012轮归档：同一alpha1模型的全域g=1/10000组件已完整认证；全积分及均值预算继续收紧，尚无新的C5界。后续候选与唯一计算入口以STATE最新条为准。上界优先，下界低资源；程序运行不预设wall上限，按有限工作量保存可恢复检查点；总研究算力不预设上限。

阅读：STATE.md → results/bounds_ledger.json → rounds/012.md → 各证书及证明。可携2048盒入口为 round011/verify_portable.py；完整97960叶的可携核验入口为 upper_portable084/check_portable_tree_v2.py，26块数据无须外部680MB状态。旧H区间沿用历史生产资格，新增1318点用不同整数算术独审；完整资格见 audit_kernel081/ 和 audit_globalg083/，主接门见 round012/globalg_acceptance.json。

只有新下界严格高于L_ref且新上界严格低于U_ref并通过资格核查，才记录双侧时点；此后再执行三轮完整同步迭代，再终审及编译中文论文PDF。原cube-finite-dim-bounds及原资料包保持只读，暂停状态未改；内部结果不称公开发表成果。
