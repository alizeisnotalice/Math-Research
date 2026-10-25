# 11H：固定三维tau1/3的连续积分上界

2026-10-25。cb11完整表示作者及主阅读/MATLAB通过，最终manifest7b22d25442f8d13837ece005849d2d142a18fd0e33d57543663e244de05701fa；公式独审pending。030奇性帽预算独审通过。固定tau=1/3，目标I=∫_[0,1]^2 G(1,y,z)^(3/2)≤57/100000；真实g由continuous_minimum单独处理，不重复。

复用原upper_continuous_bridge，独立bridge clone只upper/continuous_bridge/cb11_integral_*。主锁定snapshot为011_minimum_snapshot和030；作者另登记本批精确source及参数。首先选择可验证连续积分机制：建议正time-cutoff G_epsilon的无人工pole积分表示，付030/11A给出的off-line短时差与帽误差；或与minimum代理共享已验证confluent evaluator，不重复推导。有限批初始≤4096空间/时间分区诊断，固定tau不再优化；允许一次根据明确失败原因的表示修订。达到I目标优先锁证书与最小checker送原audit，未达则保存准确未闭合区域/误差分解，不把细化趋势当通过。

帽误差、短time、time大尾、空间cell包围、求积和算术误差需逐项支付。不得只用中心值或浮点Gauss阶差；不得删负cone或zero/coalescence层。主提供已证明全face帽上界（n3每face8个交点，全face总48）；single positive-face quadrant须按覆盖关系明确换算，不能重复漏算。Python/Arb实际，MATLAB主串行；共享库/表示明确，不当独审。

这是新上界关键路径，可提出更有效表示/批次，但不自行无限扩大。原agent不能扩招、提交推送、动其他目录/总账；主承担一般L1/allr/准确成本整合。当前U3=8不变，I目标本身不等于新界。
