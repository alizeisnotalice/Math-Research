# 168：158 实际复用证书的独立整数审计

唯一 run001 已完成。主授权锁定 plan `667917d040a854d59133db1aa23f3795ac1c853d13c4019205f1c0dd5da5f0cf`、binding `55ddc45a117f9d054f3059fa5d022448e548a9a45b659551b45f215861ba4bc1`。外授权 `round019-authorizations/primary168.json` SHA `2a042895e0d90ab60adbbe172f9ca3d47183b755248f9f0bb47210e2b6200ad1`。启动前真实检查无 run001、无活 audit168 进程；先前消息中断不记为数学失败。原准备文件和 165、157 冻结件均未更改。

实际 wrapper PID 20291、master PID 20299、15 workers PID 20362–20376。唯一 subprocess 返回 0，耗时 8.842573667 秒，checker 8.067415917 秒；无 wall cap、无重试。退出证据形式是已审 controller 对全部 worker 逐个 join 并 assert exitcode 全为 0，以及全部 worker 结果 complete；原 driver 没有另行序列化 per-worker exitcode 表，本报告不补造该表。终检探测 16 个 PID 均已不存在。

实际独立验证 19661 个已测试单 source 整盒证书，共 314576 新角点，以及 first-PASS 后的 1833 条 metadata skip；合计固定 21494 对。4769 个 source 的全部原 155 alias/event/body/raw-record、q/Delta/gradient、同 source 的 Delta110 恒等式及公共 log hash、1240 receiver 的原 149 闭 box/node、完整 111201 原 margin 投影均按已冻来源门绑定。130 个静态输入、98 个真实候选输出和 4 份代码的 SHA 执行前后全部相同。

实际 18799 个新 source/receiver 的 margin110 为负，全部合法匹配；它们仍按实际顺序进入 best 操作。整数路径先清除共同正分母，再对 16 个整数分子取整盒最小值，与 producer 的 Fraction 路径独立。公共 log 只继承 157 已通过的正有理证书与 SHA；没有再次求值 log 或 q。

新增 120 个 PASS115，现有 110081 条完整域叶通过 115，1120 条仍 pending115。159 个 receiver 的整盒 margin110 改善，其中 39 个仍未通过 115；失败 source 因而确实保留在 best 比较中。其余 109961 个旧 PASS115 叶保持原证书。首 PASS 后 skip、日志 hash-chain、checkpoint、全部 best 文件及新 pending 文件逐项匹配。

全 111201 叶的正 margin110 域完整保留。准确全域最小值没有改变，详 `verified_margin168.json`，与 157 原值逐有理数相同；witness 为 leaf 103387，node `[51742,"0111"]`，canonical `d1efc83ccd13c673b40c533d26c3664b283ea61d4b338f8401d8d5619ac50261`，种类 original155。这里完整 110 域与仍 pending 的 115 目标严格分开，不把部分 PASS115 宣称为全域 115 证书。

结果 `run001/result.json` SHA `fd1f6defb927907b08d7751536a74c1a9ca8bf440419dfc821d25261f443f61e`。45 份原始 worker JSON 共 2288728 字节，以必要白名单封入 `worker_evidence001.tar.gz`（699833 字节，SHA `696b0a7815df36331e4a4358194007e0177ba90c45f20801250c7a61a3b213ca`），逐成员解包后与原字节、长度和 SHA 全匹配；原文件仍保留本地。归档只做字节重建，没有重放数学。源锁与精简文件白名单见 `execution_manifest.json`。

严格适用性及负 margin、Delta 恒等式、整盒 max、H64/H128 到同一 H128 面下界的证明见原 `report.md`。本次 0 新 H/R/gradient/q/log、0 旧 corner 基线、0 P8/exp/g/Cn；无 MATLAB、Astra、下界、提交或推送。主保留全局接纳与常数合成责任。
