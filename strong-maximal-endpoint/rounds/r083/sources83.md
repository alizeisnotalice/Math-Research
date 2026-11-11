# 第83轮来源与依赖

本轮不登记未经原始文献核验的公开优先权。物理固定角点的精确重排不是本轮新定理。

- [第28轮报告§2—4](../r028/report.md)：比例r∈[0,1]^N明确包含r=0；该边界正是signed物理固定角点Hardy。其capA弱型常数等于纯Gamma标量最优值，覆盖任意非负物理输入。第83轮直接接续该已证明归约。
- [第30轮KKT与次梯度](../r030/report.md)：高度1处次梯度[1,N]及封顶Fenchel分支的既有背景；本轮求解单Gamma约束的激活根与实际费用。
- [第82轮径向费用](../r082/corner82.md)、[其常数阶展开](../r082/asymptotic82.md)、[独立证明类审查](../r082/kernel_review82.md)：有效的计费类最优上界，本轮量化其与真实算子常数的差，不把两者混同。第82轮报告和审查已补第28轮边界归属说明。
- [第82轮软预算必要下界](../r082/allocation82.md)、[第81轮平滑归约](../r081/smoothing81.md)、[第73轮cap2归约](../r073/report.md)：必要预算的cap-safe链。第81轮过强目标仍被排除；第82轮校正后的log阶目标尚未证明可达到。
- [第54轮显式有限cap2源](../r054/scalar_proof54.md)：维数64、真实有限L=10^8输入的log score>71.79，供本轮有限预算必要下界使用。
- [certify82.py](../../code/certify82.py)与[第82轮算术契约](../r082/certification82.md)：复用定向基础算术、正Gamma有限和、对数/指数正级数。只维护一份依赖，SHA256随本轮验收登记。
- 本轮 [exact83.md](exact83.md)、[asymptotic83.md](asymptotic83.md)、[defect83.md](defect83.md)、[certification83.md](certification83.md) 支持所有新数学和认证结论。更大维数的失败扩展不作为证明证据。
- [Python Decimal官方舍入文档](https://docs.python.org/3/library/decimal.html#rounding-modes)与[GitHub官方仓库限制](https://docs.github.com/en/repositories/creating-and-managing-repositories/repository-limits)，2026-11-08核对。归档规模远低于限制；没有新增大文件、缓存或重复依赖。
