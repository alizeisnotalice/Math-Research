# R070：物理非线性截帽与有限递归小核

[完整证明](../../proofs/070-physical-feedback.md)和[独立重建](../../reviews/070-physical-feedback.md)证明：给定可计算解析κ，二维指定纯宏观canonical构造族有每个预选有限阶的兼容递归及严格有限输入恢复。这里的实际有限计算核验代数，不承担几何证明或全输入最优性。

复现：`python3 physical_curve_check.py`（仅标准库Fraction，order=6、depth=4，写同目录JSON）。Wolfram 15中分别运行 `wolframscript -file physical_branches.wl`、`high_check.wl`、`physical_j1.wl`、`scalar_rebuild.wl`。各wl均独立，不依赖R069运行文件。

作者Fraction核保留完整参数曲线prefix，沿κ=1、ρ=(1+2u)^−1路径核至u⁶；Mathematica独立恢复首层六个系数、Φ/Ψ互逆、两支接点及有效域内非负性。reviewer另核严格预算根和曲率增益16/729。它们实际成功运行的原日志均保留；未纳入第一次因SymPy不可用而失败的脚本。日志中既有编码显示瑕疵原样保留。

这些检查是短小代数核，未声称15核持续满载或机器区间有限恢复已执行。任意阶/截断误差由解析证明给出，不由有限jet推断。`manifest.json`绑定全部公开字节。
