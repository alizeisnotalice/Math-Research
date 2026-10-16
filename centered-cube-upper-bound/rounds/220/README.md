# 第220轮：一般尺度共同重叠与软度选择费用

一般 n^(1/2+o(1)) 一致弱型上界仍未闭合。本轮把R219同尺度几何扩展到全部实际尺度：b>0、L,R∈[b/2,b]、A∈[L,b]、cj≥min(L,R)/√LR，任意独立双软度s,σ及任意来源y,v,z下，归一共同重叠D≤1。三条实际目标腿全部覆盖，不要求硬来源z属于终端Zg。

张量化给固定参数真实接收积分 ∫√(P_s,L(x−y)P_σ,A(x−v)) h_R(x−z)dx ≤(cj/g)^n√w_g(y,v)。原硬门仍精确保留为两个坐标乘积之差，量纲与来源除数完整。证明结合饱和、两种尺度顺序、远位移与中间位移凹性；不是从有限积分样例猜一般参数界。

最后标量V(t)≤(1+t)^−2采用计算辅助严格证明：完整解析公式与单调性，加64份r右矩形、40个t闭区间及纯Fraction指数包围。每格的严格归一上界均≤23/25<1，覆盖全部t∈[0,1]。独审完整核源码、全部40格、解析接口及实际重跑；所有证据可复现。

另证公共软度网格：sup_s H_s≤4Σ节点H，实际积分内选择费用≤32+16√n log(n−1)（n≥2）；节点独立来源位置。它仅支付固定其他参数时的soft选择，不能把Q、capture、history或连续物理尺度藏入核。

## 剩余缺口

有限费用桥明确列出原Q下界、全硬来源质量、clock端点质量比及所有j/尺度菜单的真实价格。一般输入下这些合计尚未获根阶预算，连续尺度选择也尚未支付。因此本轮登记一般尺度核估计和soft包络，**不登记新的原正交通根阶付款，不宣称一般弱型闭合**。

## 实际核验

- 677行详细LaTeX完整独审通过；保留文本，不编译PDF。主代理核56个标签、55处引用及环境配对，修正β任意时的sup量词与一处符号歧义。
- 普通探路：120个真实尺度配置、240次多精度积分，最大D约0.65594；9个泛仿射配置、18次多精度网格检查。两组均不是一般性证明。
- 严格证书：纯整数/Fraction，40格完整有理覆盖；与前两组普通求积明确区分。
- 主代理从最终checkpoint启动driver，再迁移脚本至新的临时目录：15项manifest哈希前后通过，三份完整JSON逐字节一致。运行时间分别4.091、0.135、0.179秒。
- 外部回执：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round220_20261016_driver_replay/receipt.json`。

## 精简归档

[research.zip](research.zip)：16成员、98077字节，SHA256 `8adb0d0a96a3598b3e4f5756544bb8d31fb2a92243842d979e75c9c2c5c5f69d`。包含完整证明、解析接口、独审、三套脚本/JSON、数值笔记、报告与manifest。

展开备份：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round220_20261016_final_checkpoint`；独立ZIP：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round220_20261016_research.zip`。所有工作文件与两份备份逐成员字节核对后才移除本轮工作目录。仓库只新增本README和ZIP，另更新项目README；未改历史或其他项目。

解压后运行`python reproduce.py --verify-only`；完整复现用`python reproduce.py --output-dir <新的包外目录>`，普通积分需mpmath1.3.0，严格证书只需Python标准库。提交和远端接收以Git及推送确认记录为准，归档存在不代表推送完成。
