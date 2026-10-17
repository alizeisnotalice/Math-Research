# 第255轮：共同中点抵消与粗尺度来源能量

本轮证明任意非负有限L1来源、任意不交尺度窗的窗数无关平方增量界，平滑尺度为 d=a/(8√n)：

Σ||Δh*μ||² ≤ 4exp(1/96) log(b/a) √[n(log(n+1)+4)] ||h_d*μ||²。

由此得到连续octave最大平方能量 O(√n log²(n+2))||h_d*μ||²，及原来源分式抽取的根阶质量费用。一般弱型上界尚未闭合；高剩余来源的真实贡献及原交通全账仍需支付。

详细LaTeX源、正测度证明、独立审核、数值与复现材料均在 [research.zip](research.zip)，不编译PDF。解压后使用Python3+mpmath运行 `reproduce.py`：验证核心文件哈希，逐字节重放JSON，运行审计与辅助稿内嵌检验。

本轮有限检查：124个完整majorant、1184个有理Bernstein分片、6个真实L1来源/366窗；14个大维数纯公式价格，不冒称物理来源模拟。独审先筛命题再读证明。一般证明与有限证据分开记录。

包外完整复现已通过：10份核心SHA前后核验、数值JSON逐字节一致、独审35项及辅助24案例/97层/8425块均通过。归档11个成员，共57,751字节，SHA256 `ddd7ea13f1bd2eb1e76aaefb39eef9a941d157ea638c0e90cecf392acda97cbd`。包外检查点：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round255_20261017_final_checkpoint`。
