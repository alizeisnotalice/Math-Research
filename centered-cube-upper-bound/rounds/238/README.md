# 第238轮：原来源拥挤几何与真实接收占用预算

日期：2026-10-17。一般 n^(1/2+o(1)) 一致弱型上界仍未闭合。

原全尺度 hardband 强迫接收立方体避开质量过大的小来源 cube，得到真实 Lebesgue 接收占用 δ(z)≤min(1,n dH(z)/a)。与实际 soft overlap 插值、原 FIRST 和完整后验相结合，全部 j/M 节点及 soft 菜单求和给出 C*∫δ[R_D+log(1/δ)]³dν，其中 C*=12Λe⁴(7e/8)、Λ=128/119。δR_D³≤A 的原来源分支已以 C*(4A+108/e³)||f||₁ 支付；A=√n 为根阶分支。互补来源仍未付款，不宣称该分支必非空、占正比例或覆盖一般输入。

## 冻结证明与复现

[research.zip](research.zip) 保存11个成员：完整 LaTeX 文本、几何与 overlap 两份必要证明附件、独审、合同、报告、数值代码/结果/说明、复现脚本与 manifest。不编译PDF。

- ZIP：44,157 bytes；SHA256 `5c4be7e02b98a628e494c7473539376da571f960f96b86d7e56bde91516a476a`。
- 独立审核通过原 all-scale hardband 来源、量词、centered 半径、所有节点常数和费用；没有认证 FULL 原分支非空。
- 本机正式69项有限检查：24精确box、8真实有限密度接收测度、6扩大维度、20费用、5驻点、6完整级数。
- 另外实际重跑9项真实核插值与1项自适应归一压力检查。
- 外置完整复跑成功，主 JSON 44,036 bytes 逐字节相同，SHA256 `13caa12b5d2ce84855ae91b09ef6570fdc43b08741ae061f20bcc158de689b84`；manifest10项在复跑前后均通过。
- 低维例及大维 Decimal 数值均按其范围标识，不替代一般证明。

解压后使用 Python3.12 与 mpmath1.3.0，执行 `python3 reproduce.py --output-dir <包外不存在的新目录>`。本机依赖路径、命令与精度口径见 numeric_notes.md。外置冻结目录 `/Users/zhengzhihao/Desktop/T/output/centered_cube_round238_20261017_final_checkpoint`；完整复跑凭据 `/Users/zhengzhihao/Desktop/T/output/centered_cube_round238_20261017_final_replay/receipt.json`。

本轮新增 tracked 文件仅此 README 与 ZIP，另更新项目总览；自己的临时 work 在外置备份和逐成员校验后移除。新附件微信路径失效，本轮未声称读到其内容。下一瓶颈是局部 soft 高质量而 hard 不够集中的来源与接收区域的联合几何，原一般目标保持。
