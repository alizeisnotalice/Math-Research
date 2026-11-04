# 角向合并的精确算例

正文paper/angular_grouping.tex给出锥概率测度的极坐标公式，以及同次数角向权的Lp范数表示和Minkowski等号判据。端点命题限于(H.1)的正双侧可比权，不凭一般Lp恒等式扩大端点适用类。锥测度不能替换成任意曲面的均匀面积测度。

`python3 angular_certificate.py`仅用标准库：n1–3、δ1/2与3/4共6个有理矩配置，各3组正系数组合，精确核对平方系数及严格合并下降；Decimal显示不承担证明。`wolframscript -file angular_verify.wl`为Mathematica15的6个精确算例复核。两脚本均实际从归档路径运行通过。一般公式及等号判据由正文证明，不由这18个系数组合推断。
