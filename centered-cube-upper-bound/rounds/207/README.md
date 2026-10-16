# 第207轮：近线性帽、来源时间乘积与clamp面价

一般 n^(1/2+o(1)) 一致弱型上界仍未闭合。本轮的几何进展是降低来源帽幂次、保留时间与尺度关联，并估计此前单列的clamp额外面速度。

1. 原非线性硬面核满足 N_B≤65B√n min{√[B(1−s)],√log(e+B(1−s))}。共同分数指数优化将B^(3/2)费用降为B乘平方根对数，保留真实面测度与完整尺度归一。
2. 原t≤g与current比例带给合法逐尺度帽B_L。其相对尺度积分按q(ξ)=∏max(Ci,ξ)精确分段；m=#{Ci≤1/2}时unclamped路径获得至多1/m因子，m0仍保对数。原γ=1+1/n的relative clamp平台另保logγ。
3. 保原时间权后，线性帽乘积精确化为B0 I(D;x)/∏Di，Di=max(Ci,ξt)。小于时间阈值的个体倒数在这个线性因子中取消；非线性对数因子及原e_t归一仍有来源费用。
4. clamp额外速度b/2的正核路径价≤B*，B*=B0/∏Ci，没有额外维数因子。证明用G≤BP、真实一维边缘归一及n(γ−1)=1；不将它当作来源预算已付。

尚未支付：真实来源归一和对数尾、receiver-adaptive尺度/soft选择、共同α与capture入口、储备补充及其它端态和门交换。所有βΠ²、j/time/window费用保留。核界与原n阶能量分开作Cauchy仍可能只有n^(3/4)，一般根阶目标未完成。

本机有限验证包括9个θ优化案例、11个精确分段产品案例、2个非法域守卫、20个Fraction clamp区间；mpmath25/50位部分为非区间诊断。n128低阈值例积分为(1−2^−128)/128，n1024全Ci1例仍为log2；没有据有限数据猜阶数。

复现：解压research.zip后 python3 reproduce.py --verify-only；完整运行 python3 reproduce.py --output-dir <新目录> --numeric-python <有mpmath的Python>。使用已有Python Fraction与mpmath 1.3.0。保留详细LaTeX，不编译PDF。

主代理最终driver实际exit0：11文件哈希前后均通过，完整JSON逐字节一致。634行详细LaTeX全文终审和结构检查通过，未编译PDF。

[research.zip](research.zip)：58,855字节、12成员，SHA256 `ac828c1435fa04c11728f3a12fa295392594cef8c1c9f72c2a544dbf1c068f3b`。仓库外展开备份与压缩包均逐字节核验后才移除本轮自身work。仅新增说明和归档、更新项目索引，共3个仓库文件。复用1 Astra high、2 Sol high、1 Sol medium，使用B03逐点支配纪律；无逐模型实际扣额数据。
