# 第043轮：准确分数矩及受限定的列截断障碍

本轮没有新的Cn上界。已完成并独立ACCEPT的是指定独立坐标乘积、统一列截断及其moment尾部支付结构的范围限制。当前已接受原算子界仍为R042；R045的新逆列权重方案正在独审，不能在此升级为新定理。

固定r,ρ,g>0且g≤min(r/2,ρ/2)，自由源坐标X均匀于[0,ρ]，ell=1−(X−ρ+g)+/(2r)。对全部实p>1，精确矩为1−g/ρ+2r[1−(1−g/(2r))^(p+1)]/[ρ(p+1)]。均值m=1−g²/(4rρ)，U=ell/m满足EU=1，EUlogU≥g³/(64r²ρ)。N=n−1自由坐标归一乘积Zhat=∏U，其logmoment的凸性给log E Zhat^p/(p−1)≥Ng³/(64r²ρ)，无需把固定p展开的不明误差除以p−1。

因此，若只用E[Y1{Y>HEY}]≤EY·EZhat^p/H^(p−1)来支付固定θ比例的尾能量，这一种证书必须满足logH≥Ng³/(64r²ρ)+log(1/θ)/(p−1)。p趋1不能消除三次项；第二项可在确定性Y的例子中纯属Markov证书代价，不能混为真实几何必要性。

另对冻结slice的Y-product envelope作energy tilt dQ=Zhat dP，令A=Ng³/(r²ρ)，则EQlogZhat≥A/64、VarQlogZhat≤A，Chebyshev给Q{Zhat≤exp(A/128)}≤16384/A。该结论只诊断Y自身的能量律。R042只证明真实列C_A≤Y，没有相等或正反向下界；因此绝不能由Y的截断必要性推出真实C_A能量的必要性。作者初稿84d版本在这一桥上表述过宽，本轮审稿指出后已作最小scope修正，最终d8fe版本的21个display公式字节完全未改，独审确认ACCEPT。旧草稿仅在仓库外保留。

在g~d,r~1,ρ~n/L且证书缺额形式≤Cnd/H的明确假设下，上述moment支付模板只能认证O(nL^-1/3)阶。它不限制真实缺额、Cn、任意kernel、mixed-radius law、actualμ或row-dependent weighting。boundedA固定Y-energy比例η还要求H≥η，补齐相同阶论证。R045改变并重新计算行权重预算，正属本轮未排除的路线。

全解析机制及证明见fractional_moment_cut/PROOF.md；原算子一般L1/全尺度定理没有在本轮重证或更改，且本轮模板命题也不声称其自身是算子定理。作者实际10项Fraction矩核对在scope修正前运行，只验证未改公式；独审另实际执行实p符号导数/归一矩展开及新512bit Arb熵点核验，未重跑作者检查。实际运行各有receipt/source hash，不把旧操作日志当作新运行。

必要来源及复核代码保存在本轮MANIFEST与audit/round043/fractional_moment，数值点不替代全参数证明。公开结果比较不变，本轮没有内部新Cn界，也不称内部成果公开发表。

同步说明：正常fetch/rebase/push已恢复，R031–042与24页论文分别保留独立提交，远端已读到38786ecb。rebase后提交ID会改变，原记录中的旧ID保留实际溯源；可移植恢复原13页稿应使用稳定Git blob：TeX a2867a0957281e932d260de6e4cea40319b619f1，PDF cbfa926d3e77e744d1c2d9a212e37621f68d0c17，执行git show BLOB即可。主clone出现外部生成的未跟踪“ 2”同名兄弟文件，未把它们当作研究源或暂存；读取诊断遇到阻塞已停止，不改删这些来源未明文件。
