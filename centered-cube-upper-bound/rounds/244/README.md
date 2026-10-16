# 第244轮：共同来源—接收L²容量与交通嵌入缺口

日期2026-10-17。一般n^(1/2+o(1))弱型上界尚未闭合，本轮没有新的一般交通付款分支。

合并R243真实带权列容量与原ALL-R硬帽，得到显式masked Schur界：||TακF||²≤C0q Cα||F||²ν，其中Cα=1+[1−(a/b)^(αn)]/α。F必须是跨receiver/node共同的原来源函数；全部mask、测度与费用有明确范围。该式是已审几何容量的直接推论，不宣称解决原交通。

实际嵌入原母式及half-column后，首个未付项是逆尺度权receiver平方及完整后验归一化；原source门F_y的平方还需按β/clock/score积分。仅用总receiver预算吸收εX会有ε^(−H+1)损失。原后验槽不能替换硬来源，负信用不能重复使用。R243全补域二次来源价与旧界保持有效。

## 实际核验及归档

[research.zip](research.zip) 含11成员：详细LaTeX、主嵌入/receiver账附件、独审、合同、数值脚本/JSON/说明、报告、重放驱动与manifest。不编译PDF。

- ZIP 49,083 bytes；SHA256 `198addc87b322884a320b357c41c94f2ad98a6a699da632e66ab5e144bc89f68`。
- 8种真实连续正L1盒几何、23项算子配置，含初探、端点与扩展；14项全Fraction LDL半正定证书及有理Rayleigh核验。矩阵与完整连续L²算子由cell-mean归约严格等价。
- 普通40/70位主要字段最大相对变化约5.812e−41。数值小范数不视为一般增益；不认证原FULL/history或根阶交通付款。
- 另27项附件代数检查通过。包外完整重放103,560-byte JSON逐字节一致，SHA256 `afae3727b466464b3386ddda1fe6939cac3f6a8316bb602dd0b810c5a1b1657b`；10项manifest前后检查通过。

解压后执行 `python3 reproduce.py --output-dir <包外不存在的新目录>`，本机Python/mpmath路径见numeric_notes.md。checkpoint：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round244_20261017_final_checkpoint`；回执：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round244_20261017_final_replay/receipt.json`。

仅新增README和ZIP两个tracked文件，另更新总览。自己的work在仓库外checkpoint、独立ZIP备份及逐成员字节验证后移除。已有球/热核文献的窄读记录另保仓库外，不算新文献成果、不重复归档，也不作为立方体证明。后续重点仍是原联合平方预算的支付。
