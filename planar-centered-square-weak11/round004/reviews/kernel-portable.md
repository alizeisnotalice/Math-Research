# R004 portable-v1 便携层终核

2026-11-05。结论：**PASS，无修改要求，可以按当前白名单接入。** 本轮只审便携化差异、身份、历史receipt与移位行为，没有重审已签数学结论或重跑C4。只写本审查目录，原包及数学产物未修改。

审查对象为 `scratch/round004-kernel-public-package/round004/kernel`。外部白名单 `package-freeze.json` 的实际SHA256为 `de445e1c271a83be438d1946ef50c59b903b673bbff9c375b4121984cacff659`，与指定值相同。实际文件集合、逐项SHA和大小均与白名单一致：16文件、789091字节。内部manifest为 `775d3dbbbf6b86d21af44036c60f8499714761b266bfa58498e7fc51b9bd94f1`；适配verifier为 `df87ec8c6cb4c8a6da22e7b62316851e6be82951464c617c889081cbf31cc8c3`；receipt为 `3fef586ca85a44b995aed1a2be1e0c9b69ae3d022cb8ddb1b9fa702a49a486e7`；execution为 `438f6b0b457344ffed70826fc2c0b0caeeb3e204e6bf2c7d50516403cef26777`，全部实核一致。

已全文读取适配脚本并对原已审 `computation/round004-kernel-direction/verify_result.py` 作逐段diff：候选精确身份/a00变化、其他数学字段、C4 expected集合与计数/无重复、每worker历史provenance、443757叶、实际192位C1/C2、独立完整轨道质量及精确端点、射线负点控制全部保留。增加了13项static和5项共享文件哈希、9项原字节copy绑定、stdout逐cell核对、样条多项式有理积分、实际cell/leaf计数和C1端点。不存在以打印float替代证明或弱化全域条件的变化。

输入定位全部经脚本文件位置找到 `project/round004/kernel` 与 `project/kernel-rescaling`；脚本没有打开历史绝对路径，没有要求本机FLINT二进制摘要或路径等于旧主机。历史worker摘要仍对历史主进程provenance核验，当前环境另录。`-O`在读取/计算前拒绝。默认只打印；显式输出保护所有静态输入、5共享文件、manifest及原receipt/execution。

`source-bindings.json`中9项原字节副本均与源实际字节相同，原verifier与适配版SHA正确；原数学说明和独审文件的SHA也与本地已签版本一致。13项静态输入与5项共享依赖集合明确，缺少原全部freeze、旧C4、探索点、监控wrapper均已说明，未冒称重验遗漏项。README保留“128位新完整C4记录/192位默认便宜检查”的区别，明确完整报告是逐cell汇总而非逐叶证书，并给完整C4可选重放命令。

历史execution与receipt的manifest/verifier/result哈希、自报加载路径及当前/历史provenance差异一致；receipt的候选、报告、C4数量、精确C2下端点、质量上下端点、claim、半质量差和负点上端点逐项等于已签数学运行。非交付移位runner字节哈希也与外部白名单一致，代码确实构造仅需交付static及5共享文件的project、从`/tmp`执行并核结果。它们是便携运行记录，不是又一次C4运行。

为直接核便携行为，本次又在本审查目录构造独立project副本：仅16个交付文件和5个绑定共享文件，从`cwd=/tmp`以`-I -B`运行默认verifier。4.48秒返回0，全部上述数学字段与原已签运行和历史portable receipt完全一致；实际加载checker/core路径确位于新副本，当前provenance与历史不同。另实际运行`-O`，返回1并给出预期拒绝信息。未重跑C4。运行前后原16文件哈希保持不变。

证据为本目录 `portable-audit-result.json`、`portable-independent.stdout.json`及空stderr日志，核验脚本为 `audit_portable.py`。独立smoke stdout SHA256为 `785ec2402a002401c0620c5db6b7ebaa186f070da9f5f200fede649486f95134`。

接入只复制白名单中的 `round004/kernel`，共享5项从已发布 `kernel-rescaling`复用；外部builder、relocation runner、临时project和本审查smoke证据不属于交付kernel目录。原数学主稿与审报由项目另行保留原字节并映射发布路径。便携包不改变已签的一般L1、全部正尺度上界 \(C_2<3.590865312933\) 的范围。
