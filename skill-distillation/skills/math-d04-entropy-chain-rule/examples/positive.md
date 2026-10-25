# 适用案例

非零KL核验例：P=[[1/10,3/10],[1/5,2/5]]，Q=[[1/4,1/4],[1/8,3/8]]；两者均为满支撑2×2联合分布，边缘和条件项均非零，直接核验KL(P||Q)=KL(P_X||Q_X)+Σ_x P_X(x)KL(P_{Y|x}||Q_{Y|x})。

本例已由本轮独立执行；实际输入、条件、结果及当前 Skill 哈希见 [本地证据解析说明](../references/evidence-guide.md)；bundle-relative artifact：`audit_current/cd/cases.json`（由说明中的resolver命令解析）。该实例检查不等于一般定理证明。
