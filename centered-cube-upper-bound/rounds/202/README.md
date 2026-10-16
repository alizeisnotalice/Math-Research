# 第202轮：原捕获归一化与来源重叠估计

一般 n^(1/2+o(1)) 一致弱型上界仍未闭合。本轮从真实ΠJ改写加权后验，避免额外ν(V)/e放大；条件方差与改变中心的偏差分别保留，原α及实际target系数不重选。

证明材料偏差由真实selected/lost卷积重叠o控制的一阶/平方显式界，并在原2ρ·αJ/a≤1/2下给归一化分量价格。该界不需全来源位置矩，保留B*、β、j/time与窗口费用；o≈1的一般中间区域仍有线性一阶费用。原capture、变帽权和selector未付。真实时间熵及删除差已记录，但没有熵到材料Dirichlet耗散的已证桥。

最终本机driver实际核验11个文件哈希并完整重跑：3个Fraction证书、3个物理25/40位积分、8个物理MC案例及其同样本材料/overlap扩展；完整JSON逐字节一致。平移与材料得分分开，经验o代入只是诊断；求积、插值和MC均不冒充区间或统计认证。详细LaTeX747行，全文独审和静态检查通过，未编译PDF。

[research.zip](research.zip)：69,499字节、12成员，SHA256 `2494df09b9f7e2e05d748c68f31de8e88ea74cf413baa084efe67cd58a01bd3e`。解压后 `python3 reproduce.py --verify-only`；完整复现 `python3 reproduce.py --output-dir <新目录> --numeric-python <具numpy和mpmath的Python>`。本机依赖NumPy2.3.5、mpmath1.3.0，具体解释器/PYTHONPATH见numerical_review.md。未使用SciPy。跨环境字节差异须检查实际原因。

展开检查点已在仓库外备份并逐字节核验，归档成员亦核对一致。仅新增本说明与归档，更新项目索引。复用1 Astra high、2 Sol high、1 Sol medium；缺少逐模型实际扣额数据。下一步集中攻o≈1区域的材料偏差与真实capture/selector联合费用，不以小重叠分支代替一般闭合。
