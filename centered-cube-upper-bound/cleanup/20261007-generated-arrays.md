# 2026-10-07 数值产物清理

移除 280 个可由已保留 guard 脚本生成的全量 NPZ 数组，共 211,737,524 字节。报告、生成代码、JSON 摘要和每个系列一个代表样本保留；未删除数学报告或代码。移除文件已备份在仓库外的本地归档目录，完整清单及校验值见同名 JSON。

涉及四个系列：projection_overlap_spatial_core、projection_overlap_spatial_core_strengthened、ordered_source_weighted_dual、bernoulli_weak_endpoint。其 guard 脚本直接调用 np.savez_compressed 生成这些产物。需要完整数组时，在独立输出目录按原报告的固定参数运行生成脚本；不要直接覆盖已有结果。

本次清理采用新提交删除文件，不改写 Git 历史。旧数组仍可能保留在历史提交中，因此这次主要解决当前目录条目过多及后续产物膨胀的问题，并不宣称已减小全部历史体积。
