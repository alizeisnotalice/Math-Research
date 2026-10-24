# 第二优化轮最小复现

本轮仅浮点发现，未产生候选，未运行新的严格核验。实际账本见 result.json；参数及退出条件见 PLAN.md。

在新的本地临时目录中复制本目录的 search.py 和 formula_modules（不复制 result.json），将上一轮 round020/round01/scaled_certificate.json 复制为 input_baseline.json，将同目录 scaling_arithmetic.json 复制为 input_scaling_arithmetic.json，然后运行 `python -B search.py`。依赖：Python3.12、numpy2.5.4、scipy1.18.1、highspy1.15.1。程序拒绝覆盖既有结果；此次真实运行无需重复。

三个 formula_modules 源文件原样取自固定 Tomar 上游 bca281a706bbf6375b04d601ad19ba4ad0bf28b6，许可见 ../round020/upstream-v4/LICENSE。它们仅用于浮点探路；不是严格核验器。基线与许可证引用上一轮，不重复归档。

HiGHS time_limit=8 是求解器软限时，实际run最高9.006714秒，不能声称严格8秒终止；外层每组25秒进程硬截止和总300秒预算均未触发。无额外重试。实际计时不包含顶层模块导入0.130822秒；即使加上也远低于总预算。
