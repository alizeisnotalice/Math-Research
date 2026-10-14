# 第122轮：合法共同测试的核超水平集与统一捕获

日期：2026-10-14。基础 UP-E/UC 两稿与四份数学审计、数值协议及其审计、实现说明与代码审计均已冻结。三阶段实际记录和三次隔离重放已完成：六进程退出码均为 0、stderr 均为空，阶段原字节与隔离重放一致；完整输出为 142,112,481 字节。最终执行说明已纳入，打包入口自身验证 PASS；独立包审已完成：执行证据与发布包均 PASS，无阻断或待核项。三个阶段分别为 10/27/43 条 JSONL，总计 80 条。阶段三 n=32 两个 η 配置三档均记录 NP2 预定充分证书失败，抽象付款资格未认证，不计为通过；该失败不否定一般引理或证明整个盒在 E 外。

最终登记 UC 一笔 ε_G/C_G 空间价格，整体替换 RR-U 的 C 与 UP-E 的 E\\C 两笔费用。原完整负项保留一次，其他旧费用与独立槽继续保留。两类低核密度存活补域、低 pe、初始/跳跃低测试与大删源耗散尚未付款；一般根阶目标仍 UNPROVED。有限几何与抽象交通不认证完整 FIRST/FULL 或实际新增流。

## 归档格式与复现

package_round.py 只读取已有证据，不执行 numeric.py。研究 ZIP 仅含 proof.md、proof.tex.txt、inventory.json 三成员；LaTeX 仅为文本，不编译 PDF。proof.md 按原始字节数与 SHA256 内嵌六项资产：numeric.py、numeric_protocol.md、kernel_superlevel.md、union_capture.md、kernel_numeric_output.jsonl、kernel_execution_receipt.json。两数学稿是 producer 启动时的哈希依赖，必须随源码和协议放同目录；执行回执给出实际 Python 环境。

准确 CLI 为：

```sh
python3 numeric.py --stage1 > stage1.jsonl
python3 numeric.py --stage2 --prior stage1.jsonl > stage2.jsonl
python3 numeric.py --stage3 --prior stage1.jsonl > stage3.jsonl
```

阶段二、三均引用阶段一，三个 stdout 按序拼接应匹配 kernel_numeric_output.jsonl。隔离运行阶段一输出保存为 prior_stage1.jsonl，后两次读取同一 prior；它是文件名，不是 --prior_stage1 参数。逐条 record_sha256、真实源哈希、完成标记、六次回执与原字节拼接均由打包入口核验。零退出、完成菜单、原字节复现仅是执行事实，不能合并未判定或资格未认证结果为 PASS。

只发布基础 kernel_superlevel/union_capture 范围；后续 weighted/suffix/softness 材料不纳入正文或资产。归档必须小于 100,000,000 字节；若完整证据压缩后超限，只报告超限，不拆包或截断资产。实际打包结果见 archive_integrity.json，已核三成员 ZIP 完整性、ZIP/本地成员一致和六项资产准确提取。

## 最终归档

[research.zip](research.zip)：**44,070,449 字节**，小于独占 100MB 上限；SHA256：`43f963f2c2aea00fb17a0308c0599f4a6bc0732caea74925b9a422723817a1cc`。正文 70 个数学显示按原字节及原序保留，proof.md 完整保存 142,112,481 字节原 JSONL 资产，未截断、拆包或改动冻结证据。打包没有重新执行数值，没有 Git 或 PDF 操作。
