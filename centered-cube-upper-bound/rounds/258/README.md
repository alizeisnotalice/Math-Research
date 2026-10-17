# 第258轮：全软度近来源边界付款与高标签正份额

新增经独审的预算：原R223 signed账按同一时钟不变近来源指标分割后，整个原时钟区间的后验边界价≤32Λ/√3 sqrt(K_nκ)Wpair=O(√(n log(n+2)))||f||_1。覆盖全部原soft，不再乘score n或clock长度；正负、删除和有利储备同权分账，旧后缀负信用不恢复。

另一新增正份额：原软核展开中|S|≥ceil(2logn/|logζ|)的高标签尾，其次概率衰减支付原score和clock，全部soft适用，尾门可为空。原时钟score与人工核标签严格区分；没有少标签的一般score支配。

一般目标仍未闭合：近来源后缀联合恢复、低soft低标签、互补远来源及全局账仍待支付。详见[research.zip](research.zip)中的LaTeX文本、独审和复现；不编译PDF。

主代理包外完整复现PASS：3条真实C阈值有限law、15个有理中点/6个有理删除检查；289个40/70字段核验；JSON逐字节一致、10件前后哈希一致，嵌入11/12/65项检查通过。law采用log4的toy时钟，不是原octave/FULL实例；普通高精度不是区间证书。

归档11成员、50,161字节，SHA256 `1196646842dc806945152dc5564ed27b9f3384979c011da6d73479d419a47c02`。外置检查点：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round258_20261017_final_checkpoint`。解压后Python3+mpmath运行reproduce.py。
