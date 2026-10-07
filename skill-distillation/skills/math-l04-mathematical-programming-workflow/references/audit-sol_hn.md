# L04 蒸馏与局部审计

审计者：GPT-6.1 Sol。本记录读取真实证据卡并执行所述有限演练；full_read只记录源阅读与覆盖，不是独立证明认证。来源行号指转换原文，准确卡/页码由provenance回查。

## 原文结论与适用量词

P-0af5d66d0492bded LP limitedprecisionoracle lines376–525、547–674、721–937、1229–1465是rational 细化/reconstruction与最终exact verification；P-666b3fde09c43a4d exactconic对偶基于real algebra。MPAX/SDPNAL/HDSDP/CCOpt/LCQPow给数值算法、KKT/convergence各范围。

## 证明骨架与已检查推导

complete模型→判convexity/integer→solvercandidate→independent原问题/对偶/rationalcheck；锥规划弱不可行性需要精确链/face信息，浮点KKT不是exact证书。

本次实际检查：max x+y≤1,x,y≥0；(1,0)原问题1，对偶 multiplier1upper1，opt1。非凸驻点不能以status排除更好点。

## 正反例与失败处理

正例：最大化 x+y，约束 x+y<=1、x,y>=0。原问题点 (1,0) 的目标值为 1；对约束 x+y<=1 取对偶乘子 1，也给出上界 1。

条件缺失例：非凸模型的求解器返回一个局部值，却在没有全局证书时被报告成全局最优；应将其标成候选值。

## 中心立方体迁移推断

用优化软件探索有限参数族，再以显式数学检查验证抽取出的候选。

迁移前逐项给输入映射、工具前提为何成立、界的方向、n/scale常数和连接引理。公开原文没有自动提供这些接口。

## 待验证命题

弱不可行、退化半定规划、互补约束局部平稳性各不同；容差结束条件不会自动变globalproof。

求解器只处理输入的模型。约束放松、方向写反、容差或局部最优都可能改变数学结论。

行为结果只认证记录中的有限输入；部分演练通过不表示高级定理全部证明或一般cube主张验收。旧文逐段审计保持局部scope，不认证同段其余结论。
