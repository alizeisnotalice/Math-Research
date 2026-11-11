# 第88轮来源与依赖

- [r028/report.md](../r028/report.md)：固定比例、固定反射中心的真实算子常数与标量常数等式；平移不改变费用。本轮明确区分固定物理端点与固定反射中心。
- [r075/report.md](../r075/report.md)：重锚会改变接收者所见源函数，不能由每点恒等式推出共同源预算；本轮所有列始终满足同一个 $G\le F$。
- [r081/smoothing81.md](../r081/smoothing81.md)：平坦cutoff及硬/软比较；[r084/obstruction84.md](../r084/obstruction84.md)：仿射积源的最优物理端点和薄层定价。本轮重新计算移动物理端点的反射中心Jacobian，未将旧定价直接套用。
- [r085/core85.md](../r085/core85.md)：近Q原点反射中心壳可以由一个匹配反射列真实覆盖；[r087/corner87.md](../r087/corner87.md)：固定物理角点列在该主体低效率。二者量词不同，本轮不把障碍混同。
- 本轮 [core88.md](core88.md) 给完整上下构造和外扩概率体积公式；[audit88.md](audit88.md) 独立审查严格阈值、计数、零集、Borel可测与真实列费用。原生Astra核心推导，主代理整合验收，不作未核验的公开优先权声明。
- Luna承担只读数值代码审查和最终验收脚本整理；其静态内存风险结论与实际运行证据区分。计算源为 [experiment88.py](../../code/experiment88.py)、[verify88.wl](../../code/verify88.wl)、[check88.py](../../code/check88.py)、[math_guard.py](../../code/math_guard.py)，复用既有有理区间实现；SHA256随验收登记。
- 一般 $Q_N\le C_N\le R_N$、对数系数间隙 $(\mu-1)/8$ 与既有 $\log C_{64}>71.79$ 未改变。新的相容中心公式是可检验的充分接口，没有预设所有输入存在低费用proposal。
