# 第607轮：真实混合律的全行动gap核查

一般 n^(1/2+o(1)) 一致弱型上界仍未闭合。第606轮两条固定混合律的全行动gap下界占能量约0.99824、0.77114，均严格超过1/16。因此不再把它们作为low-full-gap分支的证据；原有限混合收益仍有效。

一个另行核验的n128实际方体实例说明：八个可达phase不提供八份廉价容量，必须保留真实g⁻¹成本。该例W项容易付款，不是完整目标反例，没有新的统一预算。

[研究包](research.zip)含详细LaTeX、最终独审、两份必要脚本和结果、SHA256SUMS。不重复存放第606轮代码，不编译PDF。

复现：先从[第606轮研究包](../606/research.zip)提取numeric.py，运行 `python gap.py replay.json /path/to/606/numeric.py`；另运行 `python phase_probe.py phase_replay.json`。需python-flint。实际使用本机Python3.12.14/python-flint0.9.0，gap用Arb512，phase用Arb256及Fraction。根代理从历史ZIP提取后重放gap，全部数学字段一致；两计算均已执行，独审核验量词、方向、真实几何及源付款边界。

ZIP SHA-256：`6dc27790a77af3d611af94e0c08b9ede05bde17d53ccf495c422ae20fc55c8b8`。

剩余缺口仍是一般输入的完整联合几何预算，包括单phase静态核心、真实半径费用和全部交叉项；有限策略的高gap不排除同输入的其他策略。
