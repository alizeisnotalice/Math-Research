# General-R流式试验的数学接口

本稿是第38轮接口在R4/8/16/32的显式推广。两路预运行独审先核对程序，冻结PLAN后才进行一次四组试验。普通FP的完整向外误差未认证；下列实数算术证明与保存数值的证据资格分开。

## 来源与实际prior

K255、B4、v8、w4、n6704800、共同D与38相同。表b的plus格微相位是(label/4+(b+.5)/R) mod1，minus格是(label/4+(.5−b)/R) mod1。使用整数odd numerator/(2R)，R32允许箱峰边界恰碰粗格边，不跨出。每格峰质量相同，D对所有表相同，所以真实类别质量η共同，但不是1/2。

三个以上分段的整数切点按名义nη+σΦ⁻¹(b/R)四舍五入后固定。这些是固定整数，不是精确二项分位点。对共同实际Bin(n,η)段概率分别保留BE及η参数误差上下包π_b^-、π_b^+，不把它们归一到1，也不直接取1/R。gated源在第b段选第b表，周期质量准确为Σπ_b=1。不同R相位网格不嵌套；R4原子兼容38，但本轮R4切点不同于38的±1sigma切点，不是38的重放。

## 条件矩和PB

给定有序query bins，隐藏位置独立；目标对坐标置换对称，因此count histogram可用规范有序代表。对每个表b、尺度j的真实类核均值m_±，令a=m_++m_-，p=m_+/a。准确条件均值是Σ_b A_b T_b，其中A_b=∏a_i，T_b是该表PB计数落入自己的第b段的概率。gated项不再乘π_b。prior才用Σπ_bA_b。

输入仅以counts N_q表示后，四个充分统计是ΣN_q log a_q、ΣN_q p_q、ΣN_q p_q(1−p_q)、ΣN_q log(1+r_q)，r_q为该表两个class的相对方差帽最大值。所有维数都由这些有限和准确表达，不用n项PB递推。

两class均值各有绝对误差e，概率参数误差≤e/a；PB均值与方差误差均≤ne/min a，产品均值的log误差≤n[−log(1−2e/min a)]。宏观真实参数在[.3,.7]，独立非同分布BE常数.5606按38轮已核对的文献给CDF误差≤1/√n。整数半格CDF差的两个端点各付一次；首尾支持端点准确为0/1。每表使用自己的相邻两个切点，未将一个表的PB段概率套给其他表。

39/38已证的正多项式分层Minkowski给绝对方差上包(Σ_bu_b√r_b)²，u_b为该表gated或prior项的均值上包，r_b=exp(ΣN_qlog(1+r_q))−1。固定表同样可用class最大帽。均值区间与这个方差帽进入Cantelli：各尺度下包取max，上包按并集求和截到1并结合Markov。这里始终是共同j先组合全部表，之后选择尺度。

## Proposal及两侧帽

冻结N256。每个R独立seed，先以1/2概率抽均匀counts；否则在全部RJ表尺度对中均匀选(b,j)，按名义a_q=exp(log a_q)显式归一后抽Multinomial(n,p)。记Z_bj=Q⁻¹Σa_q，A_bj(N)=∏a_q^{N_q}，理想混合RN为
\[
 D(N)=\frac12+\frac1{2RJ}\sum_{b,j}A_{bj}(N)/Z_{bj}^n.
\]
用的是混合后的RN，不是被抽中component的RN；multinomial系数消去。Z没有强制设为1，实际nominal向量和normalizer定义在采样与stream中相同。

以确定source计算C_Z=max(1,max Z_bj^n)，C_E=exp(max logmean inflation)。各模型每j均值上包≤C_EΣ_b A_bj，D≥Σ_bj A_bj/(2RJ C_Z)。条件尾下包B为正时λB≤均值；上包U有union/Markov控制，故对两个方向都有
\[
 0\le\lambda B/D,\lambda U/D\le
 \min(2\lambda,2RJ C_ZC_E).
\]
C_Z/C_E在抽样后stream才汇总，但其输入、计算函数及使用规则在PLAN已冻结，完全不依counts或观测尾；没有用样本最大贡献设帽。

## 统计及实现边界

四个R、三个模型、21个固定h、两个方向，共504个单侧Bernoulli-KL区间，总名义α=.05。区间只覆盖这些冻结阈值与有限233尺度；不成为全部阈值或连续尺度的上包。额外L=1只为weak平凡锚1，不属于原233点网格。

保存逐样本lower/upper贡献、logD、整个component和的log、proposal标签。可以独立重算均值、ESS、每项CI、同h paired描述SE及同时CI端点差。没有另对paired统计声称未分配α的CI。计数与全部原始核统计只保留hash，不能从小payload完全重构每个RN分量；完整同环境重放才可重建，它不是独立确认。

程序只分配N×Q counts及4×N×R×J统计，不建立4×R×J×Q特征。int64 counts在hash后转float64，每项≤n<2⁵³，整数转换精确。源的三角函数、CDF、BLAS、随机提议及所有FP尚未完整向外认证，经验guard不作证书。四次pilot全程未改冻结文件、未按收益增样本。

周期来源转真正L¹输入仍采用已证的整数周期截断与内域因子，固定维数再令截断尺度趋无穷；本试验没有支付匹配渐近源资格、全部阈值或连续尺度费用。因此未产生新C_n阶数结论。
