# 显式矩阵算术（不等于捕获结论）

令 A 为 4×4 全 1 矩阵，A=vv*、v=(1,1,1,1)，故 A 半正定。对 x=(1,1,1,1)，x*Ax=16，而只取对角项得到 4。该有限例说明交叉项不可丢弃；它不声称 A 可由私有方向捕获几何实现。

已实际核验，见 `audit_current/g/case_results_20261007.json`（G03-positive-explicit-PSD-quadratic）。
