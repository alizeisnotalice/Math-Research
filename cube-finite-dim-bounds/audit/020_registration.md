# 第8轮接口预审：原点消元、有理化与双层超周期

2026-10-25。已读主rounds/008-registration.md及results/round008_input_snapshot.json。目标门槛保持015完整准确上端、原R4完整上端、U3=8/U4=16。本轮只新audit/020_*起及必要literature；不旧证据/总账/提交推送/扩招，非唯一工作者，MATLAB由主串行。

先做有限解析/准确Python预审，不等待作者搜索：从B的准确原点值和导数推n3/n4二元轨道原点消元关系与Hessian；明确全胞元ψ非负/原点等式的有理化不能靠缩放修复。对固定若干含负端点整数块及两接收层，独立枚举parity split、总质量2、质量斜率和共同第四根余项的适用条件，验证符号/跨超胞元表示。尚不重算旧R4体积，不跑任何LP/相位/自由参数优化。

预算：一批独立Fraction/SymPy有限恒等与块枚举（最多一次实现修复）；作者cb8/vg8稳定源到达后另登记精确hash、只补读新公式/接口和必要反例测试，不重复旧认证。实际新界优先；没有稳定候选时只给严格接口/必要条件及未闭合项，绝不把D>8/16误记Cn下界。

## vg8实际稳定源到达后的补核（计算前）

锁vg8_superperiod.md SHA 82e9a5628938bb50352b06dfad845ec1411ad62b3f759455c844b10aef44a857、vg8_parity.py SHA 42ac3828f335f75ec46c8c0ef7c6e5708895e304007b96284c95b79fb5a9ded8、vg8_blocks.json SHA 05b2f7804dfdf8d16749a844fe7ed0361379e1ac02ecfd1f714763f398ffb551。独立读取固定12输入块，另用逐原子奇偶分类与floor/ceil公式核16键、391实例、3质量及层平移；准确N81负斜率q4095/4096例独立核八raw接收端点余项。不导作者parity/guard/checker，不重算R4联合。已读锁源全文；MATLAB仅主实际回执。
