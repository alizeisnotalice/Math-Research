# 第211轮：原far余项的高初始参数基线与恢复付款

一般 n^(1/2+o(1)) 一致弱型上界尚未闭合。

本轮直接在原实际最终行上定义alpha0=Dq/QJ及rho0=(1-kappa/a-alpha0 J/a)+。保全部实际门、来源、时钟及已删已付子集不变，完整far正余项精确拆为该基线与正恢复；原尚可用负余额只留基线一次。无需构造合法空捕获历史。

恢复量由同源capture和far时钟积分控制，利用原完整history继承的heavy门T e_g>sqrt(ell)q，得到行价D/sqrt(ell)。在原外测度上使用独立预留的epsilon_rec槽及有限维母核补齐，恢复费用为epsilon_rec X+C_rec sqrt(n)W_A。不得复用旧捕获误差槽，或按各share重复支付。基线正支满足T Delta>D K0+(Dq/QJ)J>=D(K0+J)，仍待全局几何预算；near等未涉原块不被本结论覆盖。

辅助结果包括保完整归一化的对数恢复界、同钟链先合并后一次log，以及保真实分母的有限字段极值。一般来源在有限正窗的可积性和绝对连续性已论证。自由字段极值不宣称物理核可达。

本机Python精确计算包含9个单步模型、1个两步链、3个D校准、3个字段模型及纯代数cross-min校准；alpha0修复只复用前三个既有终态。Fraction与显式有理log尾界检验通过，有限条件模型不认证完整空间/history或一般定理。详细LaTeX不编译PDF。

复现：解压research.zip后运行python3 reproduce.py --verify-only，再运行python3 reproduce.py --output-dir <全新目录>。证明、来源链、误差额度和剩余问题见包内proof.tex.txt、report.md、audit.md及numeric_notes.md。

最终归档13个成员（12个实质文件及清单），91,547字节，SHA256 `1f49107d4bb570aac04ac474b3ee6ed68ee555e20cbb5388e70cdf282266006f`。955行详细LaTeX经全文独审与静态检查。主代理实际运行最终驱动：前后12文件哈希通过，迁移重放exit 0，完整JSON逐字节一致。展开材料及zip均已在仓库外逐字节验证备份；仓库仅提交本轮README、research.zip及项目总览。
