# 第260轮：相近来源稳定性与移动载体费用

已独审：原重叠核C的全局log-Lipschitz常数12；近来源原权比exp(±δ)，δ≤3/(4sqrt(n))。共同载体条件后验TV≤tanh(δ/2)，固定辅助来源score的对数时钟切换积分≤δ。完整推导保留两个禁块、两个移动距离窗、各自实际度数与零度数交换行。

本轮没有新增统一根阶交通付款。完整后验失质量仍需与原score、history和真实接收测度联合支付；已有原硬捕获外尾不能移植到辅助后验壳。一般n^(1/2+o(1))上界尚未闭合。

[research.zip](research.zip)含LaTeX证明文本、报告、独审、数值代码与结果、复现及manifest，共7成员，15,449字节；SHA256 `4f04ebae774a3e22db3030f3723f4bac47c1b94a5e0d1c210cee2f27e0531f4d`。不编译PDF。

主代理本机实际重跑L12：485点、2408点对、60个精确有理TV例通过。40/70位25个摘要字段最大差约6.60e-41；非区间、非FULL证书。reproduce.py结果JSON逐字节一致，6文件哈希前后通过。抽象退化例有明确适用范围。

外置检查点：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round260_20261017_final_checkpoint`。解压后用Python3+mpmath运行`reproduce.py`。
