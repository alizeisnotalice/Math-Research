# 独立覆盖核对：90秒有限预算

仅在export标记完整时运行一次run_verify_once.py；90秒外硬包含导入/读写，现存证据拒绝重复。上限1500叶/24000条准确有理行。不调用Arb、AS_model、APsi或LP，不重新计算生成元。

对每个父格用全部叶的x端点切成有理slab；每slab内的y区间须从0至1准确相接，无缝隙或内部重叠。另核每parent面积恰1及local/physical坐标一致，三格物理面积合计3/1024。该方法不同于导出器四叉树分账，不靠面积相等单独推覆盖。

核固定baseline文件hash、alpha/R/h；所有accepted模型hash/parent/local/mode/rows绑定。用模型保存的有理Bern端点及非负余项，重新组装c/s两列，再核每叶16个不同Bern索引的准确c1*c+s与保存余量相等且严格正。const与Taylor按整叶固定mode选择，不逐系数混拼。pending/unresolved/in_flight必须空。

这个独立检查覆盖几何与有理行算术，不独立证明其Arb解析包围，也不把局部覆盖当全域证书或新C2界。失败保存状态，不修数据或自动重跑。
