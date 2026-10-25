# 第609轮：双来源坐标轮廓与共同半径约束

一般 n^(1/2+o(1)) 一致弱型上界仍未闭合。证明了任意可测半径图的双点乘积包络；特定对角位移满足 2^n J/√n→√π。但完整位移积分严格等于1+2n ln(b/a)，不能把单点对的平方根尺度外推为一般预算。

有限L¹均匀输入、符合原bulk条件的中心族中，实际常半径核行质量为1，逐点对包络却为1+2n ln2；其总量至少线性于n乘W，实际核总量不超过W。损失在于放松同一中心的共同半径和冻结标签；未反证一般目标，也未排除保留这些约束的其他双来源方法。一般维数费用没有改善。

[研究包](research.zip)含详细LaTeX、最终独审、必要脚本/结果及SHA256SUMS，不编译PDF。未归档重复草稿或额外反例。

本机Python3.12.14、SymPy1.14.0、python-flint0.9.0实际执行18行有理数/Arb256轮廓核验与3行有限输入公式核验；主代理重放后数学和版本字段完全一致。渐近与全部维数的行积分由解析证明确定，有限表只核实现。没有采样n维数据或优化一般行动场。

解压后复现：`/Users/zhengzhihao/.local/share/math-research-tools/venv/bin/python profile.py replay.json`。比较时只排除command、runtime_seconds、peak_rss_bytes_macos。全部输入、严格区间端点及版本在JSON。

ZIP SHA-256：`5eeb8ad59c6beb0e173669084f8afcc1df6c32d7ac88cb94e70acd5a7ca07cdf`。下一步须保留同一父/中心的共同半径、真实质量窗与跨层源预算；原完整母账和历史费用不变。
