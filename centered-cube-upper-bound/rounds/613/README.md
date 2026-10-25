# 第613轮：共同尺度降维的实际费用核验

一般 n^(1/2+o(1)) 一致弱型上界仍未证明，本轮没有完整上界改善。

低捕获层的精确坐标展开产生随尺度漂移的单侧区间，不能免费使用固定中心的低维弱常数。即使假设每层已有√k费用，逐层独立最大权之和也严格至少为3⌊n/2⌋/20；这是指定收费架构的损失，不是实际算子的下界。停止这条独立层付款安排。

共同尺度分块给出有完整常数的未闭合递推：n=2m≥64，d=⌈3√n⌉，O_n≤16(2d+√(πn))+αD_m，α=(256/225)/(1−16e⁻⁹)<√2。O_n是单octave常数；D_m是低层块与完整块使用同一t的同步张量弱常数，仍作用于2m维输入。尚未证明D_m≤O_m+A√m，不能把递推当作已降维或一般上界。原低层接收集及所有一般输入费用仍保留。

实际验证：Python3.12.14、SymPy1.14.0、python-flint0.9.0；五维标量费用的精确有理/Arb证书，129个层由主代理另用320位有理端点核验；n=64/128/256/1024四组超几何尾另用交换抽样的组合式独立重算。Taylor有理余项认证8000<e⁹<8192及α<5120/4491<√2。两个脚本完整重放稳定字段一致。有限计算不认证算子预算；一般抽样关系依详细解析证明。

[研究包](research.zip)含详细LaTeX、独审、两份必要标量脚本及结果、独立复核脚本和SHA256SUMS。解压后以用户指定数学环境运行：

```sh
/Users/zhengzhihao/.local/share/math-research-tools/venv/bin/python numeric.py replay.json
/Users/zhengzhihao/.local/share/math-research-tools/venv/bin/python split_certificate.py split_replay.json
/Users/zhengzhihao/.local/share/math-research-tools/venv/bin/python verify.py .
```

重放比较仅排除command、runtime_seconds、peak_rss_bytes_macos。不编译PDF，不归档探索草稿。团队复用1个Astra及3个GPT-6.1 Sol，任务完成后停止。

ZIP SHA-256：`428e01f254c16b3348e7981f8ce0f4ad4ec9c22c015e682f99b85ea2abaaa74c`。
