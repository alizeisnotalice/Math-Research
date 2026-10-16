# 第239轮：共同来源证人优化与全尺度联合接收几何

日期：2026-10-17。一般 n^(1/2+o(1)) 一致弱型上界仍未闭合。

纠正第238轮过宽的范围表述：原f的receiver hardband也适用于β≤f；不能转移的是ν的rank、birth及来源密度资格。原证明上界仍有效，本轮扩大合法几何证人选择，不改原β/ν/Π。

本轮证明最低尺度接收中心的共同扩张可控制全部合法尺度；缺失质量证人的最优分式价由显式Φ阈值精确刻画；同一接收中心的联合协方差给出更强排除证书。证人仅排除接收区域，不替换原交通来源。所有节点费用仍为C*∫δ[R_D+log(1/δ)]³dν，C*=12Λe⁴(7e/8)。δR_D³≤√n的条件分支按根阶支付；一般补集及远soft质量的距离权仍未付，不宣称FULL非空比例或全局闭合。

## 冻结归档及实际验证

[research.zip](research.zip) 保存11个必要成员：完整LaTeX、两份证明附件、独审、合同、报告、数值代码/结果/说明、复现driver和manifest。不编译PDF。

- ZIP：39,782 bytes；SHA256 `93855e33e307b4e4dea4b0d8d7cc2bf0f01d7fc879c64e6ab661b309ef8833a8`。
- 命题先独筛再审证明，原hardcap资格、严格正分母/L1守卫、可测、共同协方差、zeroδ与完整回代均在所述范围通过。
- 正式40个顶层配置：3真实L1接收、3边界接收、2Φ模型（各5r/4H）、8原子域外压力、24位似扩张。全部使用Fraction精确有理计算；原子等号退化单独记录，不误归于L1定理。
- 另实际重跑3组共同中心矩的内嵌Fraction脚本。
- 仓库外完整driver复跑成功，25,802-byte主JSON逐字节相同，SHA256 `955168bad62c07aa56f27adb700a030d1dd06951a34a01cbd4a9912c90e568c1`；10项manifest复跑前后通过。

解压后使用Python3.12标准库运行 `python3 reproduce.py --output-dir <包外不存在的新目录>`。完整参数与证据边界见numeric_notes.md。本机冻结目录 `/Users/zhengzhihao/Desktop/T/output/centered_cube_round239_20261017_final_checkpoint`；复跑凭据 `/Users/zhengzhihao/Desktop/T/output/centered_cube_round239_20261017_final_replay/receipt.json`。

本轮新增tracked文件仅README与ZIP，另更新项目总览；work在外置备份与逐成员校验后移除。下一瓶颈：保留真实I0(y,z)距离衰减，控制高β(D)却不产生相应共同接收的远soft来源质量，而不是再把它压成未付的总质量矩。
