# 第123轮：后缀来源对消去与共同测试耗散软付款

日期：2026-10-14。最终 current_ledger_audit.md 已判定 G123-current PASS，主代理已核对概览数学范围与证据分类。作者稿的起草待审字句保留原字节，待审状态已由最终独审解除。package_round.py 仅核现存证据，不重新运行任何 producer。源文件读取自 R122 scratch，衍生归档只写本 R123 目录，不改已发布 R122 研究包。

## 正文与独审对应

- 主线 MC：mixed_test_capacity_private.md；独立 mixed_capacity_statement_audit_private.md / mixed_capacity_proof_audit_private.md。
- 主线 CG：complement_geometry_next_private.md；suffix_sign_statement_audit_private.md 包含先行独立筛查与追加作者全文审查。
- 主线集成：weighted_suffix_integration_private.md；其实际独审范围已在 mixed_capacity_proof_audit_private.md 中登记，不另造不存在的审稿。
- 最终总账：本目录 current_ledger.md 和必需 current_ledger_audit.md；只按审稿实际范围登记。
- 工具边界：remaining_mixture_geometry_private.md，softness_mixture_statement_audit_private.md / softness_mixture_proof_audit_private.md，softness 数值说明与独审。
- 未独审失败定位：remaining_suffix_geometry_private.md 仅摘录开头状态及最终停止点，不称其结构论断为本轮已证主结论。

## 数值正文与复现资产

weighted 协议、协议独审、implementation note、代码独审、实际执行说明和 weighted_numeric_result_audit_private.md 全部列为必需。suffix 和 softness 各包含源码、数值说明及独审。十二项资产按原始字节数与 SHA256 内嵌：

1. weighted_capture_numeric_private.py、weighted_capture_numeric_protocol_private.md、weighted_numeric_output.jsonl、weighted_execution_receipt.json。
2. suffix_numeric_private.py、suffix_numeric_note_private.md、suffix_numeric_output.jsonl、suffix_execution_receipt.json。
3. softness_numeric_private.py、softness_numeric_note_private.md、softness_numeric_output.jsonl、softness_execution_receipt.json。

后两组的独立 output/receipt 资产从对应 note 的实际代码块原字节提取；不重新序列化 receipt，不冒称新运行或独立保留的六份隔离 stdout。原 note 也作为完整资产保存，支持证据溯源。weighted 校验六份独立 stdout/stderr/直接回执、聚合回执、真实 header/完成标记、input/observed摘要、原始字节拼接和隔离重放；不移用基础 kernel 的逐条 record_sha 字段。

weighted 源码与协议必须同目录，使用回执 Python 环境并设置 PYTHONINTMAXSTRDIGITS=0，三个 CLI 为 --stage1、--stage2、--stage3，无 prior 依赖。suffix / softness 均仅依源码，以 --stage 1、--stage 2、--stage 3 分别运行；各组 stdout 按阶段顺序拼接应逐字节等于内嵌 output。大有理数不截断。

## 费用与范围

最终选择先 Z signed 消去、再 Zc 正项付款；N_Z 已消费，N_Zc 和 N_del 各保留一次。数值协议的 signed 加权负项消费是互斥展示，不能拼接其剩余负项与当前总账。统一空间费替换旧费，不追加旧 UC/UP-E/RR-U/SG/RS 价格。局部工具、固定有限数值和更细余域都不声称一般目标闭合；一般目标仍 UNPROVED。

ZIP 仅 proof.md、proof.tex.txt、inventory.json 三成员。详细数学显示保留原字节及原序，复现资产完整内嵌；归档须小于100,000,000字节，超限只报告，不截断/拆包。少量归档衍生适合 Git，原大 stdout 不另入 Git；本打包任务不执行 Git。LaTeX仅保存文本，不编译 PDF。归档实际大小、SHA与内建完整性结果见 archive_integrity.json；具体数据在下节登记。

## 最终归档

[research.zip](research.zip)：**4,749,700 字节**；SHA256：`472306d1187cd9f6127c607c705882eb0e99f7d0eaaf8f2fbc1f5f3584d7fde6`。内建校验 PASS，严格三成员、十二项复现资产原字节提取和 SHA 核验通过，187 个数学显示原字节/原序保留。ZIP 小于100,000,000字节，没有截断或拆包。

最终作者主账 SHA256：`488834fb6c1ecef84b1f9e2987ea60f7c2a5b6b9a1466ad0f09d8c1a0f2f138d`，与独审绑定一致。打包没有重新运行 producer、Git 或 PDF；独立归档审查已通过：只读校验器实际运行退出码0，核验三成员、187个显示和12项资产逐字节一致；没有重跑数值producer或扩大数学认证范围。
