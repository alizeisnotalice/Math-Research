# R20 源角代数核验

这些计算只核验所列实数不等式、积分和幂次；不证明全部输入的局部几何或小角误差预算。`cost_checks.wl` 以Mathematica 15.0.0实际运行，完整结果见 `cost_checks.log`。

对x,y≥0，使用交换与符号对称归一化为0≤x≤y、0≤c≤y。若c≤x，两个截断项平均减去xy(x+y)恰为3c²(x+y)≥0；若c≥x，写y=x+s,c=x+t，余式展开为非负系数多项式。两分支的负反例集均False。越过Clip左端或右端时相应一侧两个正部同时为零；相应违例集均False，包含零坐标与等号边界。

两个目标的真实积分按Q序与两种P序分别切分，两项和均为½Σc_j P_j Q_j(P_j+Q_j)。这是一点态单目标积分和有限线性叠加的核验，不假设目标角不重复。两端点归一化为2r^(n−3)/(n v³)，再乘积分中的½得到 r^(n−3)/(n v³)，均返回True。

在n≥2整数、v>0、r>0及0<q≤p下，H=min(q,pq/v²)给出的规范化恒等式为True；H≤q、v²H≤pq的违例集为False。额外若v<1、p≤C₀v、X≥c₀H且c₀,C₀>0，则按p≥v²及p≤v²两支分别核得 pq≤max(C₀,1)vX/c₀，两个违例集都为False。

`rational_boundaries.py/json/log` 对144个有理Clip输入（24个等号例、严格左右截断各30例）、3组双目标排序与30组维数/大角分支做补充Fraction检查，均通过。有限样本不代替上述全域符号论证。小角误差项涉及外部源预算，未在本目录作独立核验。所有脚本无自动时限。

复现：从本目录运行 `wolframscript -file cost_checks.wl` 与 `python3 rational_boundaries.py`，输出分别打印/写入相应文件。脚本未改，主代理不重复既有成功批次。几何现由 proofs/020-sharp-source-cost.md 与 reviews/020-source-cost.md 独立审定；仍不以有限算例证明几何。原稿的重复六层符号检查不另复制，第21轮保留更完整的六层准确分母核验。
