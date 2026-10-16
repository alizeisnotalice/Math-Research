# 第243轮：实际接收分区与全补域二次来源价

日期2026-10-17。一般n^(1/2+o(1))一致弱型上界尚未闭合。

保留原源无关的接收纤维Ω_i，而非逐节点放大接收域，证明真实带权占用Σ(a/R_i)^(αn)O_i≤1+[1−(a/b)^(αn)]/α。配合原尺度衰减与跨节点Holder，将整个原低J补域的来源价格降为24Λe⁵ EΣ∫H_N²ν；相对旧δ=1的三次价格省去一幂。没有选择新的容易子集，没有裁剪原完整后验；旧δ加权界也保留，δ小时可能更强。

一般p<1的费用为24Λe⁵/(1−p)∫H_p²ν，原R242同一已付门降为O(log³n)||f||₁。p=1不能免费取。真正缺口仍是二次来源矩的一般根阶预算以及其它原母账分支；不能把来源价格降幂当成一般n维阶已改善。

## 实际核验与复现

[research.zip](research.zip) 含11个必要成员：详细LaTeX、主证明与一般p附件、独审、合同、数值脚本/完整JSON/说明、报告、重放驱动和manifest。不编译PDF。

- ZIP 57,910 bytes；SHA256 `8d9c2c590a03571de10a893c4dcf76957eef484c06ac8dfdf936d72580e43074`。
- 正式48项实际检查：24真实径向分区容量、12同分区联合Holder、6固定p常数、6一般p常数；覆盖初探、端点压力、扩大范围，维数至10^6。
- 40/70位主要字段最大相对变化约2.601e−34；小n由Fraction几何和直接求和交叉核对，大n采用闭式有限几何矩。普通高精不是区间证书，也不认证原FULL历史。
- 另48项附件检查实际通过。包外完整重放139,509-byte JSON逐字节一致，SHA256 `f8c6bed1bce79087760c587b010dac4355528158d528104b595e9a654277a9d8`；10项manifest前后检查通过。

解压后执行 `python3 reproduce.py --output-dir <包外不存在的新目录>`，本机runtime与mpmath路径见numeric_notes.md。checkpoint：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round243_20261017_final_checkpoint`；完整实际回执：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round243_20261017_final_replay/receipt.json`。

本轮新增README和ZIP两个tracked文件，另更新总览。自己的work已在仓库外checkpoint、独立ZIP备份和逐成员字节验证后移除；不删除原材料或他人文件。后续继续聚焦原共同来源与接收的联合几何预算。
