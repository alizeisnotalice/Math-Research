# 11Q 极窄加锐附审登记

2026-10-25。仅固定三个已审物理 case 的 A=a³N、D 系数比值及锚点 (1,1,1,0)，不重跑完整 g/表示批，不改源或账本，不提交扩招。沿用原 035 的谱表示、符号、非负 cone 覆盖和边界连续性审核。

输入为主 rounds/011q-registration.md、src/round011_n4_coefficient_ratio.py、results/round011_n4_coefficient_ratio.json、proofs/round011_n4_minimum_sharpening.md；共同 global_data SHA43be4de041c9da823d192ff6d8ef2440c3668f88b13fc52dbff320dd0011d96f。运行时逐项记录实际 hash。

唯一有限批：三个准确多项式代换，检查 A、D 所有系数及支撑并集，对 Dα>0 计算全部准确 Aα/Dα 并核最小值；对 Dα=0 单独核 Aα≥0。用差多项式 A−(17/640)D 全系数作交叉验证。一个锚点另从 1152 锥的准确时间有理函数合并及 Laurent 简单极点系数计算 C=π³G，避免仅接受“matches_known_point”布尔字段。最后只推该固定模型在 cube 的准确 minimum，不推最终 maximal 常数。
