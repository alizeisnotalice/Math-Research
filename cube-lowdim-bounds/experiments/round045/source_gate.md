# 第45轮wrapper独立静态gate

冻结`run_stage.py` SHA `7ff9ee59e16df6e8363fdaa204e571f235fcd28a4103c7e0afeca1e66b477454`，结论GO，可按PLAN分别唯一启动first/second。本审读未运行数学。

- 全文审wrapper，并静态核两实际CLI/parser。first修正cert_run+hlcert各pin，参数Q8/Adepth7/Kdepth6/6procs/MKA；second原hl_check+hl_core各pin，claim3.59087、maxdepth8/kdepth10/6procs/完整c1,c2,c3,c4，无cells或tight缩域。
- 两driver均显式将自身源码目录加入sys.path，故`python -I`与同目录依赖导入一致，不会因孤立模式漏导入正确hlcert/hl_core。工作目录为相应checker目录。
- 固定候选/structured/准确成本/独立候选保存审核/solve运行与源码逐文件pin；plain严格六字段并与structured一致，beta=alpha=47/40，成本证书绑定structured。准确upper<B1lower且upper<359087/100000后才启动，未把显示小数当门槛。
- 每阶段拒绝已有meta/log/report，只有一次Popen。独立session拥有进程组，180秒从准备开始计入launch及完整checker/报告写入；超时kill并wait，finally仅清理此自有组。终态hash记录所需短收尾不是额外数学运行，不能把elapsed展示伪称恰等于计算时限。
- source/input/wrapper运行前后hash、候选、报告、日志及PID/退出/超时记录齐全。源码或输入变化使正常终态失效。wrapper的terminal只是进程状态，尚须后审完整report；first M输出质量球，本次准确质量门槛由独立成本证书及second C3支撑。
- 完整成功仍须first K/A零失败与质量一致、second C1/C2/C3/C4及worker provenance全部通过；超时或负充分下界均不能当实际负点。两实现共享Arb及解析族，不称全部数值原语独立。
- 一般L1/全部尺度/奇异分布比较桥已审并见`plan_bridge_audit.md`。本候选沿用Tomar2.2、Lemma3.1与Lin外部存在性依赖；旧alpha=1有限源测度证明不适用于当前47/40核。只有全域条件通过才可声称新界，17/5目标仍另列。
