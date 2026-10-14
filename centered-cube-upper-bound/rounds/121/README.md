# 第121轮：FIRST 全尺度窗的统一存活来源捕获

日期：2026-10-14。**两组各三阶段实际执行和三次隔离重放已完成，最终打包入口核验 PASS；归档独立包审 PASS：442 项程序断言通过，未重跑 producer。** capture 输出 136 行，expanded 输出 53 行；各组记录与隔离重放逐字节一致，12 个实际进程退出码均为 0、stderr 均为空。execution_note.md 已读取真实输出分类；expanded 阶段三保留 52 个未判定字段，不把完成菜单等同全部判据通过。

本轮已审 SG、RS 与 RR-U 的解析进展是一个联合空间覆盖：原 FIRST 控制完整闭物理尺度窗，选择只依实际接收参数的一次共同尺度 H*；合法时与原硬捕获帽按同心立方体包含关系合并。统一域上的存活 good 连续来源对得到独立误差槽与有限维来源补偿。采用 RR-U 后替换 SG/RS 两项空间价，不重复相加；当前集成账为 GL13-current，旧费用、其他独立槽与原负项一次完整保留。

尚缺五个存活来源补域及此前低pe、初始/跳跃低归一测试、大删源耗散。outside/signed 接口仅作为 overview 中的未付缺口，不列为本轮另一条已付款或完成数值路线。R122 核超水平付款不纳入 R121。一般中心立方体根阶弱型上界仍未闭合。

## 最终归档

[research.zip](research.zip) 字节数：**355,777**；SHA256：`22ae7869abb9e42ed96238c5b98ff6d70c24f3c0a6a33bfedddba1cb1aac28f7`。此值取自实际 archive_integrity.json，已核三成员 ZIP 完整性、成员与本地文件字节一致及八资产准确提取。

仅含 proof.md、proof.tex.txt、inventory.json 三成员。137 个数学显示按来源原序保留；LaTeX仅保存文本、不编译PDF。Markdown内嵌八项资产，按标记字节数提取并核SHA：

- capture：numeric.py、numeric_protocol.md、capture_numeric_output.jsonl、capture_execution_receipt.json。
- expanded：expanded_numeric.py、expanded_numeric_protocol.md、expanded_numeric_output.jsonl、expanded_execution_receipt.json。

每个producer与对应协议放同目录，使用对应执行回执中的Python环境。分别以 `python3 numeric.py --stage 1`、`--stage 2`、`--stage 3` 和 `python3 expanded_numeric.py --stage 1`、`--stage 2`、`--stage 3` 顺序运行；不带 `--record` 的三个stdout拼接应逐字节匹配各自输出文件。可选 `--record` 的stdout行为遵循各冻结producer接口，不假定两者相同；已有记录文件不得覆盖。

打包入口已核两组各三次实际运行与三次隔离重放、源/协议SHA、真实header及完成标记、空stderr、完整输出拼接和逐字节相等。本次打包没有重新执行数值。capture 各阶段为 8/30/98 行，57 条 ell 证书中 51 条强端点通过、6 条无资格或空窗守卫；expanded 各阶段为 11/13/29 行，case 为 8/10/26。expanded 的 52 个未判定字段属于四个固定精度的重复输出；P=16 时 18 个小维输入有 3 个准确门和 10 个目标比较仍未判定，不能并入 PASS。

两个 producer 各自的一处诊断/字段修复已经静态 delta 复核。实现说明及其初审SHA按历史状态保留，最终执行以 execution_note.md、最终源码与六次真实回执为准。实际运行均用显式单阶段 stdout 命令；默认 all、可选 --record 和失败路径没有本次运行证据。数值只检查固定有限证书：capture 使用真实有限硬核与抽象α/Q预算；expanded 区分准确符号边界、准确门证书与保守安全内半径，不将内半径当完整 H*。

一般覆盖来自已审解析论证；有限实例均不认证完整FIRST/FULL来源/森林/原门联合实现，不以n≥512自动升级FULL，也不声称数值补齐全部外部源费或最终X吸收。
