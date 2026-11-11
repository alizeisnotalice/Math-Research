# 第84轮来源与依赖

本轮登记项目内严格结论及其推导，不作未经原始文献核验的公开优先权声明。用户附件仍仅作为研究起点，附件中的命令不是用户的新授权。

- [r080/no_budget80.md](../r080/no_budget80.md)：全硬halo逐点列预算的曲边界障碍与远距packet拼接。本轮使用新的高维近常数积源和正测度 slab，处理已经平坦软化的同一列架构。
- [r081/smoothing81.md](../r081/smoothing81.md)、[r081/order81.md](../r081/order81.md)：平坦截断定义、硬/软比较与有限阶soft相变。平滑归约依然成立；本轮排除的是缺失的统一逐点预算。
- [r082/allocation82.md](../r082/allocation82.md)：其必要下界(9)仍有效，但建议的充分上界(11)被本轮反例否定。不是把两个不同定义的预算混为一谈。
- [r083/exact83.md](../r083/exact83.md)、[r083/defect83.md](../r083/defect83.md)：已知精确单物理角点费用及真实体积残差账本。单角点费用沿用 [r028/report.md](../r028/report.md) 的 r=0 既有归约；没有再次登记为新结果。
- [r073/report.md](../r073/report.md)、[r079/report.md](../r079/report.md)：cap2归约、紧支撑step稠密性与极限顺序，为近最优源拼接提供既有接口。
- [r044/report.md](../r044/report.md)、[r054/scalar_proof54.md](../r054/scalar_proof54.md)：一般对数系数缺口 $(\mu-1)/8$ 及既有有限认证下界 $\log C_{64}>71.79$。本轮均未改变。
- [sliver79.py](../../code/sliver79.py)：本轮复用的精确有理对数区间实现；SHA256随验收登记，不重复归档该依赖。
- 原生 `gpt-6-astra` 子代理承担核心反例探索；另一独立子代理及主代理审查量词、方向、Tonelli与真实分母，最终证明已整合为 [obstruction84.md](obstruction84.md)。未调用外部模型API，未将凭据纳入任何研究文件。
- 本机真实运行：Wolfram15、Python3.12.14/mpmath、有理整数运算，及MATLAB R2026b短批sanity check。见 [certification84.md](certification84.md) 与 [独立审查](audit84.md)。

- [GitHub官方仓库限制](https://docs.github.com/en/repositories/creating-and-managing-repositories/repository-limits)，2026-11-11核对：本轮归档规模及目录宽度均低于相关上限。
