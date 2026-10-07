# 非零 mesh 误差：独立补充终态

本补充在主批 169515 项终态之后另行冻结，未改变、重跑或累计进主批计数。

三轮 n=4,64,256，每轮 ξ=1/2,1/4,1/8，共9记录。完整原来源 μ=δ0，receiver x=(7/10,0,…,0)，原连续窗口 [1,2]。闭 cube 的原响应在 R<7/5 时为0，在 R≥7/5 时为R^(-n)，所以唯一连续赢家严格为 R=7/5。取原 Δ=ξ/(2n)，向上最近 mesh 点 r=1+ceil((R−1)/Δ)Δ。所有本次参数的 Δ 分母是2的幂，R含因子5，故 r>R，误差严格非零。精确检验 Mmesh=r^(-n)<Mtrue=R^(-n)≤(1+ξ)Mmesh，并保存完整 source/receiver、赢家、mesh尺度与全部精确比值。

注册后新 script 一次运行，exit0，**91 PASS**（9×10算术检查及1script hash）；无失败，oracle0。主批9个 saved 比较仍只检验退化为相等的左端点赢家；本补充仅认证非零包络费用，不证明一般连续 maximal 定理、nearmax 或 actual FIRST。

- code SHA256 93db47b8be7445b763b8a4f658c5c3fa307b5c3311f408c4f386498b418bdf30
- registration c88eef4999f9a8cb0a7a990588cdb0519b75cec5622788c240ae0e4e26d66694
- results 2359b7e37fc640269ec35e8c953cfbf07f1e880553f2565b2d5400a680063817
- receipt 29171268a1a310c8165c32609434ebbdc9a6172da6562ef3235dfec4c5892fc8

主批和补充是两次独立注册/执行，分别报告169515与91项，不用合计冒称一次执行。

