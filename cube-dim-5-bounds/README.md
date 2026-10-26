# 五维中心方体弱型常数

固定Lebesgue、中心轴平行、共同半径、一般L1、全部正尺度。

有效界为 C_cert≤C5≤32，其中 C_cert 是约1.9943164615684859的完整代数构造，严格大于19943/10000。准确表达式、来源和公开/内部比较见 results/bounds_ledger.json；不从展示小数重造证书。接续门槛冻结在 round005/continuation_reference.json。

第011轮归档：符号基闭盒证书、零核树转换、完整S2管预算及正时间分段积分组件；没有新的C5界。后续候选与唯一计算入口以STATE最新条为准。上界优先，下界低资源；程序运行不预设wall上限，按有限工作量保存可恢复检查点；总研究算力不预设上限。

阅读：STATE.md → results/bounds_ledger.json → rounds/011.md → 各证书及证明。可携2048盒入口为 round011/verify_portable.py；完整15143新父的数据在 upper_portable071/。底层H区间共享生产输入，树转换独审见 audit_conversion069/final/。

只有新上下界均严格超过接续门槛并通过资格核查，才记录双侧时点；此后再执行三轮完整同步迭代，再终审及编译中文论文PDF。原cube-finite-dim-bounds及原资料包保持只读，暂停状态未改；内部结果不称公开发表成果。
