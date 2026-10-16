# 第196轮：同储备目标联合补偿与占用回代

一般 n^(1/2+o(1)) 一致弱型上界尚未闭合。本轮独审通过同一原 alpha 下的实际 reserve—剩余目标联合二次补偿；恢复 target 标签后精确接回原 CC22 加权占用，并将有界部分与真实尾部拆开。该补偿不把储备和目标分别取独立最坏值，也不按每个目标重领来源预算。

同时给出原 good 条件下完整分母的精确标量外放松，以及高 A 核心的可比 current 费用与真实 high 列恒等式。没有把这些局部公式登记为新增一般覆盖域。统一几何 majorant、带 w/k 和逆 target 质量的占用尾部、剩余比值平方列都仍未支付。

本机最终driver验证13个文件哈希，完整重跑两脚本，两个JSON逐字节一致。Fraction三阶段检查71179个标量约束点；联合二次界9组参数、39个标量检查、117个抽象预算回代通过。它们是有限精确代数证据，不是物理来源可实现性或一般维度阶数证明。

[research.zip](research.zip)：33,362字节、14成员，SHA256 `bff95189ec0b6ac5bd1807fa2c0fdcb28ee40dae9bb0046eb85b25de80b0f91d`。含详细LaTeX文本、最终独审、必要证明与复现材料；未编译PDF。解压后 `python3 reproduce.py --verify-only`，完整重放使用 `python3 reproduce.py --output-dir <新目录>`，仅需Python标准库。

仅新增本说明和归档，更新项目索引；展开文件已仓库外备份并逐字节核对。使用1 Astra high、2 Sol high、1 Sol medium，没有可获得的逐模型实际扣额明细。
