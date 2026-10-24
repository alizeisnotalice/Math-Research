# 第39轮最小复现

实际终态为首深度4未决，360接受叶、218pending；不是完整覆盖。执行前计划及实际计数见PLAN.md/result.json/extend.run.json，独立核验见saved_verification.json。

建立空临时根，按前轮复现布局放置：round01/scaled_certificate.json（第20轮归档），round17/fixed_checker/hlcert.py（第36轮补丁重建），round18/leaf_cover.json和385个final模型（第37轮saved_models.tar.xz），round19/new_model.json（第38轮），upstream-v4/checker/independent/hl_core.py与hl_check.py（锁定公开原件）。本目录源码放round20，Python3.12/python-flint0.9.0。

重新计算只能在空round20输出目录执行`python -B run_once.py`，300秒硬截止，程序拒绝覆盖既有证据。旧partial不是必需文件；只接受来源哈希一致的旧最终模型，不重算补齐。

仅审保存结果：本轮saved_models.tar.xz解压至round20，含leaf_cover.json及61个final模型；calls.jsonl.xz解压为calls.jsonl，复制result/run/source及前述历史来源。`python -B verify_saved.py`是标准库Fraction/hash/tree/分区核验，不调用Arb或积分；必须未存在saved_verification.json。`diagnose_saved_once.py`仅对保存node661的16行作准确误差分解，拒绝覆盖既有failure_breakdown.json。

archive_manifest记录62个成员及原始字节哈希；压缩结果均已逐字节往返核验。final模型足以审查，partial留本地不重复归档。历史模型与修正源码均引用前轮，不复制；原始私人资料不包含在成果包中。
