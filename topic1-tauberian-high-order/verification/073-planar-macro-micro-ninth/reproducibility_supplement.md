# R073 独立 Fraction 小核的可复现性补充

2026-11-05。仅补配对代码和一次新运行记录，不新增数学审查，不修改冻结终审。

此前 micro_increment_fraction.json 对应已执行的工具 snippet；原来只落盘输出，未保存配对脚本。现将同一有限有理计算核心导出为 micro_increment_fraction_recheck.py，并实际新运行。新脚本以 Fraction 重算所有数学字段，再与历史 JSON 的完整 payload 精确比较；没有把新 PID 或时间补写成旧执行记录。

实际执行目录：/Users/zhengzhihao/Desktop/T/topic1-tauberian-work。

实际命令：

    python3 scratch/r073-review/micro_increment_fraction_recheck.py > scratch/r073-review/micro_increment_fraction_recheck.log 2>&1

Python 版本：3.15.0 (main, Oct  9 2026, 11:16:41) [Clang 21.0.0 (clang-2100.3.34.2)]。新 PID：86148；脚本实际 UTC 时间：2026-11-05T07:48:10.816418+00:00；工具返回 exit code 0。全部有理断言通过，historical_payload_matches_exactly=True。

| 本次文件 | SHA256 | Bytes |
|---|---|---:|
| micro_increment_fraction_recheck.py | e31f77487d21a94c66e79e1a55d01fc7229fa5a23f1411e14ee82f1eea4126eb | 2306 |
| micro_increment_fraction_recheck.json | a47968a0c2d7e5a9ef8f47b9c14e0da0798ad9bded66a7073365de53bdce20b9 | 1022 |
| micro_increment_fraction_recheck.log | a47968a0c2d7e5a9ef8f47b9c14e0da0798ad9bded66a7073365de53bdce20b9 | 1022 |

历史 JSON 保持 SHA256 97554a4974e40146bde0cf33a068c06e7297a7687164ea222a7c52aa5acdd90d；冻结 audit.md 保持 SHA256 c0edf166efd8db1b4865e1e435948bd23be1437ddaafa8a46045c5337c992520。两者均再次按实际字节核验，没有覆盖。

复现时将脚本及历史 micro_increment_fraction.json 复制到一个新目录，再运行脚本；新 JSON 使用 exclusive creation，已有输出不会被覆盖。只需 Python 标准库，无第三方依赖。脚本内部数学计算不依赖历史答案；历史 JSON 只用于最后的独立结果一致性比较。

该运行仍仅核有限有理归一化：κ0=16/3、κ1=−1088/81、Λ宏线性项=−2560/81、κ2九阶乘数=0、微孔增量 8/27 与 −496/243、新增首保护单项 −40/243。几何及严格实际构造仍由未改的终审报告承担。
