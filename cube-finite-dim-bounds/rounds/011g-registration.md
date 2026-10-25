# 11G：四维纯体对角真实Green成本pilot

2026-10-25，承接u11稳定接口SHA fdfa422f2105710c84e7e43fdff1a4dcbef20dcc4e8a69d5a42116db3010b99d；作者准确锥、12非退化点与e1退化符号例通过，完整独审待续。主只跑n4纯body，不套用三维成本。

64固定scrambled Sobol内点（seed20261025，坐标映到[1e-5,1−1e-5]）及一个≤192eval局部min，体积Gauss6x7x8和10x11x12两阶，成本4^(4/3)I/g_sample^(4/3)。genericguard抛错时准确Fraction allnodes coefficient fallback，frame线∞不当0；无法求值则整阶标未完成。浮点I和min都非全域认证。

有限一批，串行避免重载；完成后无自动加密。若有<16迹象转认证，否则保留成本诊断再选新机制。主src/round011_n4_cost.py/results/round011_n4_cost.json；不可变源放scratch，完整源将随代理manifest归档；MATLAB代理固定例交主串行。无新上界预设。

首批结果14.31825触发准确单点诊断（不扩搜索）：仅用allnodes接口核数值极小趋向的x=(1,1,1,0)，保存C/pi³准确系数及Arb包围。这是g的上界证据，不是全域下界。另仅核拟定目标g_low=21/25000、I_up=19/100000的准确充分不等式，两个目标仍待全域证明。

数值best点接近重投影平面，另预登记仅该已保存点一次exact-Fraction/Arb校准，检查generic浮点取消；不重新min搜索。差异只校准该点，不认证全域。
