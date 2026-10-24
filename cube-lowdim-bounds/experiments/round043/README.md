# 第43轮复现与证据

固定B1、父(52,1)，纯有理四象限限制与192bit Arb模型。实际唯一运行27.561042459秒，232新模型盒、591最终接受树叶、1278证明区域；不是新C₂界。

将本轮程序置于临时根round24。输入沿用round01基线、round17修正版hlcert、原v4 independent源码，及round18至round23最终模型/补充证书。round23/quadrants_once.py纯Fraction函数是固定SHA依赖，仅导入函数不调用main。各准确输入SHA见result.json；旧模型归档从对应轮次展开，历史partial不要求存在。

空输出副本运行`python -B run_once.py`，外层300秒硬限，已有输出拒绝覆盖。不得在原证据目录重跑。saved_models.tar.xz含leaf_cover.json、calls.jsonl和232个新最终模型，共234成员；archive_manifest.json提供每件SHA及大小，全部已逐字节比对。没有复制旧大模型或私人原件。

只核保存证书时，解压归档并提供上述历史来源，在没有saved_verification.json的副本运行`python -B verify_saved.py`。独立核验采用幂基仿射代入，检查树覆盖和内部proof pieces覆盖；不调用积分、不使用生成器de Casteljau。
