# 161 纯归档附录

本包取冻结 manifest_v2.json 的95项与 execution002_minimal.tar.gz 的178项之必要联合。
210个原相对路径映射到165个 SHA-256 内容 blob；所有保留文件均保留原字节。
run001 原 prefix 与 run002 prefix_reused、重复 coefficients/factor 只存一次。
这是归档整理；没有 H、mean、旧 power、DFS 或其他数学重算。

`python3 portable.py verify` 校验元数据、压缩包、全部 blob、路径与大小。
`python3 portable.py restore --root /absolute/new-directory` 先全量校验，再恢复到尚不存在的目录。
恢复入口拒绝覆盖、绝对映射、路径穿越、软链接与非普通 tar 成员；不会执行所归档代码。
不同路径写独立文件，避免硬链接别名导致后续修改相互影响。

`paths.json` 为原相对路径映射，`manifest.json` 为附录清单。
原绝对源路径仍是历史证据身份，不在归档时改写。
不重复全 raw H/cache/source packets；这些外部源的身份仍由原 binding 和执行清单锁定。
省略的历史失败进度文件及 SHA 单列在 omitted_failure_files.json。
本包保留完整数学记录和真实执行记录；因省略进度／控制文件，
不宣称能原样重演失败控制流，也不自行重放任何数学计算。
