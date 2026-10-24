# 第37轮外区常数型基线准入：独立静态审查

**最终 gate：GO。** 已全文审`export_cover.py`、PLAN、`run_once.py`及`verify_saved.py`，实际读取并核对冻结文件hash。没有运行模型、核、积分、LP或核验器。许可仅是主审批准的单父格const-only有限实验，不预判能完成覆盖。

| 文件 | SHA256 |
|---|---|
| export_cover.py | `d1cec9425986e2434bf7a67fcf973ad75b68126e83d5d3d198baa1279d11ff68` |
| PLAN.md | `3b5f72fbe49376903ed72db0521cf59f4acedc04da97dd7148128979262a6898` |
| run_once.py | `4177c92064366ca45564d9cf8738b12743485319df69bfc9828fdb8650b0eae9` |
| verify_saved.py | `99a903292addd762d407d9c1ae2465ce0db4d290968a38e0b00677b1e7423808` |

## 实际数学接口

固定B1及独立闭式源码分别pin SHA；旧独立全域报告也固定hash并核spec SHA/c4成功。父(52,1)局部域[0,1]²转物理域`[13/8,53/32]×[1/32,1/16]`，每个生成盒准确断言正坐标且`u0+v0>R`。因此常数API只走外区单调分支，不受u跨坐标值R的整盒分支限制。

已确认调用为`APsi_lower_const(u0,v0,u1,v1)`，准确区分API顺序和保存box顺序`(u0,u1,v0,v1)`；没有调用内区Taylor、外区Taylor或隐式整单元checker。`AS_model(p,q,x0,x1,y0,y1)`正确使用同一knot内局部坐标。已有phi_table的端点chord处理保持；不跨内部节点假定高阶光滑。

所有模块导入后设192bit，之后重建MC3/MC3T并构造Spec/AKChecker。AS整体模型保存其内部抵消；Psi列不乘c，AS列不乘s。余项要求有限且非负；Bernstein各端点经有限精确Arb端点`man_exp()`转dyadic Fraction。每盒16行准确计算`c1*Lpsi+BernS_lower-E`，全严格正才接受，未引入数值容差或半余量。本阶段只证明基线两列安全，未来方向和LP不在此执行。

## 计数、分割与停止

根域严格为公共GRID18×18=324闭盒；分割仅准确中点四分。几何键和组件键分别拒绝重复。每开始盒最多1次AS及1次常数调用，boxes≤1500、AS≤1500、const≤1500、0Taylor/T/LP，depth≤4。首个depth4未接纳盒完成模型/16行/未决checkpoint后立即停止；不继续别的根，不自动换模型或加深。

盒弹出后先保存in_flight；组件开始、返回、失败、中断分别记录。返回后时钟终止仍只把该组件计completed，盒则中断，不伪称整个盒已完成。部分模型按阶段更新并留hash，最终node模型另存且拒绝覆盖。异常、非finite、负余项及预算退出保留未决、pending、in_flight；剩余队列不能称已检查区域。

外层300秒在启动前计时，以独占session进程组执行；包括导入、模型计算和写盘，超时杀该组，终态清理该组。内部guard不能抢占单个模型，外层硬限承担此作用。已有结果/cover/model目录拒绝复跑。

export自身完成路径检查接受面积恰为1，未独立重验每个split；完整覆盖还应以主侧标准库核验为依据。`verify_saved.py`仅读保存模型：核模型hash，按16索引重组准确行、准确中点子盒、parent/depth及接受判据；partial状态只报告实际保存模型/接受面积，complete时才另做横向slab无缝覆盖并要求总面积1。它不独立重算超越函数或AS模型，不能称第二份数值公式证书。

## 结论边界

本轮不读取或重算三旧父格基线/T模型。准入失败只说明特定常数界/多项式误差未能在既定预算内接纳已证基线，不是AK真实负证据；首次depth4失败早停时尤其不能把pending解释成一组已定位难点。完成也仅覆盖这一父格，仍不替代全域认证。下一阶段T建模或优化需另行冻结范围；目标C2<3.4尚未实现。
