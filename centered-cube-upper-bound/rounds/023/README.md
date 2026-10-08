# 第23轮：全窗口噪声与真实接收点信息账

2026-10-08。**一般 sqrt(n)n^{o(1)} 上界仍未闭合。**

证明任意有限正Borel来源在完整尺度窗口上的Poisson重采样径向链平方预算，包含此前全链估计未覆盖的顶部径向到达。建立原真赢家和实际重新选择赢家的两份合法接收点/history联合律，给出KL及TV控制。后者是原接收集上的值加权律，不能替代新弱水平集的均匀Lebesgue律。真实L1短区间输入说明免费正向KL迁移会因新增支撑失败。

尚缺实际新接收集迁移及保留多来源格的几何终态预算；本轮不提供新的维数上界。

## 归档与验证

[research.zip](research.zip) 仅含proof.md、proof.tex.txt、receiver_information_guard.py、results.json、inventory.json。解压后运行 `python3 receiver_information_guard.py`，只需标准库，确定性计算。

实际完成12项Fraction精确L1支撑检验、12项真实一维Lebesgue接收区间及连续尺度赢家计算，空间128/512/2048三档的差异均通过解析误差包络核验。Poisson尾和空间求积误差分别记录；普通浮点不是区间证书，数值不认证近极值、高弱型比或一般次数。

独立子Agent审核证明、计数尾与BV求积界；修正其指出的截断χ²遗漏尾概率问题后重跑通过。Python语法、LaTeX环境配对、ZIP CRC与SHA-256通过；临时目录解压重跑results.json逐字节一致。未编译PDF。

只新增本报告与压缩包，更新项目索引，不归档草稿、缓存、重复材料。历史依赖清理仍未全部完成。压缩包13211字节，SHA-256：60425b959695160197644d9f5f2dbbf64035d0da5c8c9444a9ac84c1a06ebaa0。

提交前扫描：已跟踪2646文件，最大单文件16150642字节，最大目录宽度811，Git目录357311536字节。本轮新增归档远低于100MB单文件和2GB单次推送限制；3000目录条目是建议值。[GitHub官方限制](https://docs.github.com/en/repositories/creating-and-managing-repositories/repository-limits)已于2026-10-08核对。
