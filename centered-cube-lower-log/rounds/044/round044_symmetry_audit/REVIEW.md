# 第44轮对称来源数值脚本独审

本目录 `checks.py` 已实际执行，`CHECKS.json` 为PASS；一般自含证明见 `PROOF.md`。没有修改主脚本/PLAN/原结果，没有提交Git或启动子代理。

## 核与标签

读取 `round044_symmetry/probe.py`。背景按θ、θ+.5切成minus/plus段，峰按符号赋相应class。b2负峰中心0的跨周期split覆盖[0,δ/2]和[1−δ/2,1]，两段均在该table的minus内。其他峰也完全位于自己的class内部；broken control负峰加宽到2δ仍不跨class边界，峰质量由反比高度保持，故仍各class质量1/2。保持class质量并不保证mirror identity，control正是在检验这个区别。

核的周期shift −3..3足以覆盖当前所有query∈[0,1)、L≤2.9；整数L分支每个class平均=.5正确。独立Fraction直接积分的336个kernel（含跨0峰、整数窗、长短窗）与主值最大差3.34e−15。主root顺序为[minus,plus]，独立checker原顺序为[plus,minus]，比较时已交换，避免标签反序造成伪通过或伪失败。

DP的第k项是plus count k：minus核不增加次数、plus核向后移一格。fixed通过两类核之和逐坐标相乘，DP系数和应等于fixed；主partition error检查方向正确。table-specific labels使每表的count函数不同，但输入质量仍由每表Bin(n,.5)质量及逐列权重和1确定。旧43 helper只接收已计算好的feature、fixed与配置，复用的fit/value是通用逐count映射训练策略，不会重新使用旧source的class定义。π用Bin(n,.5)质量加总，正确；type数量占比不能替代它。

## 扫描与高维资格

主SCAN含4 table×7 length×7 tilt×mirror/control=392项。score倾斜以dx网格为基准乘a^(1+t)，与Q=a dx一致。H=2E pp^T是原η=(.5,.5)归一化Gram；μ_t=.5时leading eigenvalue1。SCAN的mirror μ误差约4.20e−14、γ误差约3.87e−13属浮点求和量级；独立100个16384节点倾斜检查μ误差≤2.34e−15，另100个对称Fraction有限查询律精确μ=.5。任意t的结论来自证明，而不是有限扫描。

broken control的最大γ约1.605出现在短窗L=.6（q0），不受43的q≥1 κ帽控制，也不与本轮长窗9/8帽矛盾。mirror短窗τ≈.998同样说明不能把长窗帽移到q0。

高维训练冻结同一count selector后才生成独立test seed；输入不是逐点winner，允许多个table同时开启，质量证明已经重新检查。浮点guard不是向外误差证书；没有正候选时按探索统计报告即可，不能据其小数升级为严格C_n有限下界。不同table的R_b、source score与query律仍须联合付款，单表mirror结构不能证明多表无增益。

## 新增固定共同长尺度的联合界

root在真实低维grid实际探路后提出的joint目标独立审查成立，已加入 `PROOF.md` 第7节。关键是先对actual F/G作逐点Cauchy，才把各表积分分别解释为不同Q_b^n期望；没有在原共同query中逐项换元。跨表共同结构在原来源端Bin(n,.5)的Krawtchouk基，非query空间。同一count列的Σ_bw_b=1精确支付每个degree的系数平方和，保留全部score后得到∫F²/G≤9/8、χ²(F||G)≤1/8（q1）。π=0项先删除，避免0/0。

读审 `joint_probe.py`：F=旧通用value返回的真实门控平均，G=真实π加权fixed平均；separated中的component²/(π_b fixed_b)恰为A_bT_b²/π_b。主有限grid最大约1.02277、分离积分约1.07899、点wise CS无正违反，均与证明方向相符。F_mass grid偏差最高约6.74e−4说明该grid仍是普通求积，不能声称原探路有向外误差认证。其证明目标已由独立一般推导闭合，数值本身继续按探索记账。

新增 `JOINT_CHECK.py/JSON` 的全部精确小例实际通过。该联合定理不需mirror，但只适用一个共同长L和t=0原Q_b；不取消短窗或连续尺度sup的未闭合费用，也不推断C_n无log。

## finite K共同算子补充

新增条件性frame接口已独立核对并写入 `PROOF.md` 第8节：Σ_b(w_b⊗w_b)/π_b≤I，若posterior Gram存在共同G0，则联合积分≤trG0。共同单坐标H0=uuT+D经正张量与对称压缩给相应G0，费用≤∏(1−d_j)^−1。原42逐表trace帽允许明确取D=β_qP_u⊥，所以无需假设新的低trace D存在，即有有限K费用(1−β_q)^−(K−1)。该费用指数依赖K，不是K_n增长时的维数无关结论。

读审 `joint_k4_probe.py` 的共同原source-type与score Cauchy量，和此抽象接口相符；普通粗网格下小于1的F²/G数值不能当作χ²为负，仍按未向外认证的求积误差理解。没有重放该网格或声称其为严格数值积分证书。

## 稀疏长尺度推论

最后核对 `scale_budget_probe.py/SCALE_BUDGET.json`：Q/[4(Q+1)]的telescoping和、单位区间密集J项的J/8均正确。固定prior-low事件上的见证尺度满足error²/prior>λ(1−c)²/c，故逐尺度积分预算推出所提weak比较。第9节已补自含推导，包括可数集合的合法单调截断（固定无限prior事件，而非同时增加两个finite max）。连续L可测性成立，但dense预算发散，不得交换积分与sup。这只约束特殊门控族，不是原问题的全f强L1/全HL统一结论。
