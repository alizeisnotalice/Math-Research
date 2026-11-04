# 第 25 轮归档审计裁决

## 已采纳并落实

1. **计数语义。** `round25.md` 的计数列改为“**两 \(\tau\) 合计 DE 评估数**”；
   `audit_search.py` 分别打印每个 \(\tau\) 的评估数和实际原子数。三族历史总数
   仍为 7,083,648，无 diagnostic free-v1 重复计入。
2. **真实覆盖域与终止原因。** 报告及 `search_parameters.json` 记录每族原子数、
   位置盒、\(d\) 盒、质量盒和自由族的 \(10^{-12}\) 近零拒绝；明确 32 条历史
   记录均达到 `maxiter`，不是收敛声明。盒边界贴合也已披露。
3. **历史可追溯性。** 修改缓存前的三份原始脚本、原始摘要和原始日志已复制到
   `history/pre_signature_20261105/`。`historical_manifest.json` 给出相互关联文件
   的 SHA-256、32 条记录与 7,083,648 次评估；`historical_results_minimal.json`
   是不含大原子/光晕列表的最小结果表。它们显式标为 `historical_unsigned`。
4. **新版缓存。** 当前三份搜索程序使用 `run_provenance.py`，把源码和被导入
   公共源码哈希、Python/NumPy/SciPy、任务、模型盒写入每个 checkpoint。签名不合
   则将旧文件保留到 `stale/` 并重算，而不再无条件 `if exists` 复用。新版最小
   冒烟运行和一次参数签名不匹配重算均实际完成；其输出与历史日志分开。
5. **环境与运行证据。** 历史使用 Python 3.15.0、NumPy 2.5.4、SciPy 1.18.1、
   4 worker（而非 15）；复现命令和旧搜索 README 路径载入 `search_parameters.json`。
6. **MATLAB 留痕。** `weak_atomic_sanity.m` 已实际重跑，生成
   `matlab_sanity_output.json`；命令窗逐行输出也保存于 `matlab_run_output.txt`。
   Wolfram 15 的本次重跑输出保存于 `wolfram_output_v2.txt`。
7. **文字修复。** 修复了漏反斜杠的 `quad`，并将“近零位置产生高权光晕”降为
   分子层面的解析可能性，不把它说成当前数值样本观察。

## 独立裁决：`.pyc` 意见

“`__pycache__` 会静默使用旧实现”作为普通 CPython 的一般断言**不成立**：标准
timestamp pyc 会检查源码的修改时间和大小，正常源码变更会重新编译。因此没有把
它误报为已发生的数学计算错误。

但 pyc 也不是可审计的计算来源，且异常时间戳、复制和不同导入环境会损害可追溯性。
归档复现命令已加入 `PYTHONDONTWRITEBYTECODE=1`，新版签名还包含
`multiscale_search.py` 的哈希；因此不依赖 pyc 来判定结果版本。现有本机 pyc 是
可再生缓存，不建议提交。

## 归档边界

建议提交：三份当前搜索源码、`run_provenance.py`、`audit_search.py`、
`export_minimal_historical_results.py`、`search_parameters.json`、
`historical_manifest.json`、`historical_results_minimal.json`、`round25.md`、
两份独立符号/数值检查程序及 MATLAB/Wolfram 输出。保留本机而不提交：完整历史
checkpoint、完整原子/光晕 JSON、运行日志、`signed_cache/`、`__pycache__/` 和
free-v1 诊断输出。后者从不进入 7,083,648 的统计。

## 未被掩盖的数学状态

本轮只有有限盒、有限原子、浮点 DE 的负向探索；它没有给出 \(W_{1/4}\) 或
\(W_{1/2}\) 的全局上界或精确范数。新版最小运行验证的是归档机制而非重新获得
7 百万次历史搜索的数学证据。


## 主代理精简归档与独立复算

Flash只读job5ee067f0f489的计数/来源意见经主代理实际重求和与读源码核验，拒绝其普通pyc必然错误的泛化。最终audit_search.py改读可移植的32个候选输入，不依赖未提交的完整历史summary；重新计算每个候选的光晕积分/成本，复核评估计数公式及总和。历史三份源码快照有明确哈希必要性，保留以解释无签名记录；完整日志和轨迹仍仅本机。当前默认15worker，历史命令仍为实际4，归档小样本只6任务且每族2并发，不夸为满15进程研究。签名不匹配后2任务实际重算、旧记录进stale，再次匹配复用时缓存内容哈希保持。
