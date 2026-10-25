# 第614轮：同步张量的完整难度与可放宽的临界比较

一般 n^(1/2+o(1)) 一致弱型上界仍未证明，没有完整一般上界改善。

固定第一坐标块的高层乘第二块完整方体，其共同尺度强L¹费用h_m≤4√m−1。保留真实接收测度和任意相关来源，精确有D_m≤O_(2m)≤D_m/θ+h_m/(1−θ)。因此D_m仍携带完整2m维问题，不能当作已降维的m维算子。

固定分块已足以给这个比较，无需继续R613的随机分块归约。真正待证的联合比较可放宽至D_m≤√2 O_m+A√m：若成立，θ_j=j/(j+1)的完整调度推出O_n≲√n(1+log n)²，再沿用既存全尺度定位。该联合比较仍未证明；本轮不将条件递推或已付补集冒作闭合。现有独立maximal只给O_m(1+m log2)，没有改善完整费用。

实际验证：指定Python3.12.14、SymPy1.14.0、python-flint0.9.0，对m=1/2/4/8/16/32/64用两条精确多项式积分路线核验h_m，主代理重放稳定字段一致。第三条分部积分组合式与200项有理log2级数独立核七个h_m的Arb256区间；两条临界调度符号残差为零。一般强费用来自解析证明，不从有限值拟合；所有计算均未认证未付D比较。

[研究包](research.zip)含详细LaTeX、独审、必要脚本/结果、独立复核及SHA256SUMS。解压复现：

```sh
/Users/zhengzhihao/.local/share/math-research-tools/venv/bin/python numeric.py replay.json
/Users/zhengzhihao/.local/share/math-research-tools/venv/bin/python verify.py RESULTS.json
```

重放比较仅排除command、runtime_seconds、peak_rss_bytes_macos。不编译PDF，不归档重复草稿。维持1个Astra和3个GPT-6.1 Sol，任务完成后停止。

ZIP SHA-256：`9032d73104088e1c5c2e6c9872bde81e03a44228ba48f47d7f6f0dbf2612df55`。
