# 第197轮：有限位移截断核的反集中几何

一般 n^(1/2+o(1)) 一致弱型上界尚未闭合。本轮证明实际产品接收核的有限位移反集中界：固定任意软度和L/b∈[1/2,1]，若m个来源—目标坐标位移至少3b/4，则两块固定接收质量与Sperner反链给r_m≤min(1,sqrt(128/(m+1)))；真实截断平方核积分至多B r_m/(1-exp(-1/16))。证明不要求位移趋无穷，不假设来源坐标独立。

该界保留cap和原overlap费用。实际g若随接收selector变化，采用固定Wg=Cprod并支付wt/Wg。固定参数的平方根因子尚不能免费通过自适应选择、nu/e来源归一和全部j/time。一般来源预算与少分离坐标补集仍未付。辅助稿另给正负矩、显式全selector包络、远尾熵-clock和保留high的有限逆核费用，均按准确范围使用。

最终driver在本机核验13个文件哈希，完整重跑21个物理核矩（25/40位）及9组Fraction反链证书，两份JSON逐字节一致。解析截尾上界不等于浮点区间求积认证；有限证据不证明一般维度阶数。

[research.zip](research.zip)：41,189字节、14成员，SHA256 `ef3ef8cb4cda4bbd672233100d5d9be805829da07857227d14afa2abdf87727d`。含详细LaTeX文本、独审及必要复现材料；未编译PDF。解压后 `python3 reproduce.py --verify-only`；已有mpmath 1.3.0环境中 `python3 reproduce.py --output-dir <新目录>` 完整重放。

只新增本说明和归档，更新项目索引；展开文件已仓库外备份并逐字节核验。复用1 Astra high、2 Sol high、1 Sol medium；无可获得的逐模型实际扣额数据。
