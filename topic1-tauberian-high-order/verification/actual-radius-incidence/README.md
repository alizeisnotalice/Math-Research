# R27 四胞 slot 的最小组合核

`four_cells.wl` 使用 Mathematica 15.0.0 实际运行，日志为 `four_cells.log`。全部变量取实数。

第一项假设 ax≤bx、cx≤dx、ay≤cy、by≤dy、dx<ax、dy<ay。对 cx≥bx 或 by≥cy 的违反集，Reduce返回False；全称结论 cx<bx 且 by<cy 返回True。它表达两条对角中一条严格相交时，另一条在两轴上严格分离。

第二项假设 0<ε<1/10、ε≤x≤1+ε、−ε≤y≤1−ε，并有 max(|x|,|y|)≥1、max(|x−1|,|y−1|)≥1。偏离 1≤x≤1+ε、−ε≤y≤0 的违反集为False；全称蕴含为True。

本目录仅检查题述端点与坐标约束的代数后果，不证明这些条件从铺砌、前极限或slot唯一性中产生。脚本无自动时间限制。
