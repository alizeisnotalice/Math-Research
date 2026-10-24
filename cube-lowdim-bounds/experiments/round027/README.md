# 第27轮：一次全点池半余量LP

目前仅候选，第一全域核验失败；不得将准确成本3.59064620372910349562视为弱型上界。已证B1<3.59091不变，未证明3.4。最终认证终态见rounds/027.md与两个checker运行记录。

## 实际计算

21缓存hash验证，全部51574固定点、288固定方向、半基线余量；一次MATLAB dual-simplex-highs LP、两批103148个完整K/AK值和103148次各目标判定。数学Python NumPy/SciPy用于展开与回代；MATLAB和之前Python求解器都使用HiGHS，不称独立求解核心。候选唯一2^-50有理化。enclose_cost.py用整数40次根/Fraction证明准确成本下降，未证明正性。

## 复现依赖

沿第26轮README恢复scratch round01/04/05/06/07。须有原18项cache、point_manifest、round06/search.py与依赖、round07/result.json和local_artifact_manifest.json、最终K/A_directions_288.local.npy、ipc.local/result_06.json；hash不匹配即停，不自行重建或重新求解替代带计时的历史JSON。round08放本轮search.py、batch_solver.m、run_once.py，执行 `python -B run_once.py`。本机实际MATLAB路径/Applications/MATLAB_R2026a.app/bin/matlab；Python环境同第26轮。本次已运行，不需为查结果重复执行。

准确成本命令：`python enclose_cost.py candidate.json --baseline ../round01/scaling_arithmetic.json --out cost_enclosure.json`。两个全域checker命令及600秒各自预算记于round08_first.run.json/round08_second.run.json；上游版本与算术依赖沿用第20轮。两实现共享Arb和解析公式，准确成本另外独立核算。

## 保存范围

保留最终代码、事前PLAN、源码审查、结果审查、候选、准确成本、实际紧凑日志、进程退出记录及本地大文件hash清单。完整MAT/原矩阵/全值数组约95MB和第一核验336MiB叶报告仅本地。local_artifact_manifest是运行后文件清单，不是事前承诺或数学证书。原私人附件不在归档内。
