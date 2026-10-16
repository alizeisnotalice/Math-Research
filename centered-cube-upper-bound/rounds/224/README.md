# 第224轮：定量共同重叠与终端带J正交通付款

一般 n^(1/2+o(1)) 一致弱型上界尚未闭合。本轮有两项经过独立审查的推进。

1. 对全部实际尺度、独立软度和任意来源位置，精确min归一共同核满足 D_m²≤1−κ(1−ω_m)，κ=25/2^68。证明保共同硬窗口，覆盖hard-soft端点及第三J腿。固定参数张量后保score阈值，完整clock含上端薄层得到Oκ(log n)；它不是自适应接收积分的免费包络。
2. 以a_J=√J/(√K1+√K2+√J)分配原正交通，保终端带g exp(−d/n)≤t≤g。同点FIRST、高阈值、完整后验与不交接收域给P[a_J1_Td]≤εX+128Λ²/(27ε)d²e^d Wν。取d=(log n)/2得到εX+32Λ²/(27ε)√n log²n Wν，且Wν≤||f||1。全部实际L/R/s与R=b均保留，hard帽只用一次，不消费负信用。

当前浮动补域还允许在wt/wg≤√n时从good中吸收J，只剩两同soft腿。该吸收与J分配付款是可选处理，不重复领取。仍未付：两条移动目标核的自适应来源列、更深clock及一般联合几何回代。未把固定核定理或一个已付分配当作整体闭合。

## 实际验证

- 359行详细LaTeX、全部几何分支与原量回代、完整独审保留，不编译PDF。
- 本机Python/mpmath三阶段220个配置，6组40dps、2组60dps复核；最小观察gap比约0.024678145213，仅数值探路，非全域下界。
- 自足Fraction脚本：12个共同来源模型（n至10000）、3个receiver并集模型、24组Young、4个吸收例、4225软系数及所有有理常数检查通过。有限模型不认证FULL几何。
- 两脚本均经独审独立实跑；最终checkpoint迁移复现耗时约3.870秒与1.094秒，两份完整JSON分别256094与356922字节一致，12项manifest前后通过。

## 精简归档与复现

[research.zip](research.zip)：13成员、85069字节，SHA256 `2c552ceaf9bff19d483d195dab927bfb08c8f5797db1dda680f94b2454552d2c`。新增2文件，另更新项目README。来源快照哈希与最终文件绑定区别已在报告说明。

展开备份：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round224_20261016_final_checkpoint`；独立ZIP：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round224_20261016_research.zip`。工作目录、展开备份与ZIP逐成员相同后才移除本轮work。

复现：`python3 reproduce.py --verify-only`；完整运行 `python3 reproduce.py --output-dir <新的包外目录>`。需Python及mpmath1.3.0；可用`--mpmath-path`指定已有安装路径。最终回执：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round224_20261016_driver_replay/receipt.json`。
