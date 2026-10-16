# 第242轮：尺度支付的窗口几何与来源权指数

日期2026-10-17。一般n^(1/2+o(1))一致弱型上界尚未闭合。

保留真实比例节点(M,c²M)的尺度差，得到T=min(Σmin(|dy_i|/b,1),εn)与I0≤BM^-n exp(jh−T)，无需原统一窗口的(b−a)/2无衰减区。更精确的M^-n归一化允许合并J=max(E,Gsharp,T)，并在一般同源Holder中将尺度代价精确抵消为exp(−4θjh−rX)。两尺度费用分别1/θ与1/r²，未免费删去。

取p_n=max(1/2,1−1/log(n+2))<1，单一来源门J≥log+ρD/p_n只领一次来源预算，得到(21/2)Λe⁵max(2,log(n+2))[1+log(e ceilπ√n)]³||f||₁。这是O(log⁴n)的明确子支付款，包含旧门；未声称所有输入上价格更小或新增FULL交通有正比例。不能取p=1。

补域为原全部条件与E<t、Gsharp<t、T<t的交。仅在t=log+ρD/p_n≤εn时，严格推出||dy||₁<2bt；阈值饱和时不作该推断。这一联合近来源容量及其它原补项仍未付。

## 实际核验与归档

[research.zip](research.zip) 含11个必要成员：完整LaTeX、尺度/插值附件、独审、合同、数值脚本/结果/说明、报告与重放驱动。不编译PDF。

- ZIP 50,008 bytes；SHA256 `98aeeefde98f302c8510bfbc836307367d8dd315979b80ac6f08d22e5af2c85f`。
- 本机47组正式配置：22物理核/产品、19共同Holder数组、6含cap的完整级数；另136项内嵌检查通过。
- 包外完整重放通过，119,700-byte JSON逐字节一致，SHA256 `467c3ce6d67d7fb997646b929235626d7d392ac0369eaf72ad91593ac79e01b1`；10项manifest前后通过。
- 真实积分30/50位复核，产品直接用raw scalar计算log；普通高精不是区间证书，合法物理工具节点不认证原FULL历史或一般阶数。

解压后执行 `python3 reproduce.py --output-dir <包外不存在的新目录>`，依赖路径见numeric_notes.md。完整checkpoint：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round242_20261017_final_checkpoint`；回执：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round242_20261017_final_replay/receipt.json`。

本轮仅新增README和ZIP两个tracked文件，另更新总览。本轮work经过checkpoint、独立ZIP备份及逐成员核验后移除。下一瓶颈保留为低J交域的真实receiver/FIRST联合容量，不能用距离收缩替代付款。
