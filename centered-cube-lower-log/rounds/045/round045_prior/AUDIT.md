# Round045 prior：最终独立浮点回放复核

2026-10-29。首次复核新增 `replay_summary.py`、`REPLAY_SUMMARY.json` 与本文件；随后归档精简仅修改replay脚本与本文件，新增 `NEAR_SHA256.json`。没有修改既有PLAN、FROZEN、probe、SUMMARY或RESULT，首次REPLAY_SUMMARY也保持不动。实际回放用本机NumPy 2.5.4完成，耗时12.24秒，不重跑训练。

复核通过。按PLAN中的test_seed+n与原batch顺序，实际重新生成n32/128/512各32768个whole-vector proposal查询。它们是原先独立于训练的冻结holdout的确定性回放，并非新增独立重复样本。未读取query/features缓存。

独立实现使用“完整周期质量+剩余弧交叠”计算真实周期box窗口；它与原probe的周期primitive差分采用不同的计算路径。稳定logsumexp仅用NumPy实现。对每个向量，q始终为uniform项加48个整向量anchor乘积之和；没有逐坐标mixture乘积，也没有只除所抽中component的密度。所有proposal component计数与原记录完全一致。

对全部3维×5模型×11阈值重新累加sums、sums2、命中数、raw/guard/near、RN质量，并重算样本方差、经验Bernstein半径与截断端点。原SUMMARY没有直接存sums/sums2，因此对照量分别由N×mean及(N−1)variance+N×mean²恢复。收费方向数重新核实为336；使用原α=.05/336。

全部165组guard命中数精确一致；全部near query索引逐项一致，总数分别12316、33995、64753。全部模型的完整阈值profile最大值与最大值所在阈值也一致，不只检查训练选中的h。

|n|完整profile的mixture最大值 (h)|全部fixed完整profile最大值 (模型,h)|差值|
|---:|---|---|---:|
|32|.958931007 (0)|.971039562 (fixed1,.5)|−.012108555|
|128|.827906267 (.5)|.912358405 (fixed2,.5)|−.084452138|
|512|.650224855 (2)|.684050368 (fixed3,4)|−.033825513|

这些是冻结有限profile的描述性比较，不能改称连续阈值最大值；也没有依据holdout重新选π。复核只读取FROZEN中的π，并核验它恰对应已保存训练网格argmax；三维训练选择均为π=(0,0,.5,.5)。

独立计算路径相对原SUMMARY的最大绝对差如下：

|量|最大绝对差|
|---|---:|
|IS均值|7.55×10⁻¹⁵|
|sums|2.48×10⁻¹⁰|
|sums2|1.80×10⁻⁸|
|样本方差|5.47×10⁻¹³|
|EB半径|2.67×10⁻¹⁵|
|EB端点|9.33×10⁻¹⁵|
|RN均值|3.22×10⁻¹⁵|
|RN的EB端点|3.34×10⁻¹⁵|

sums/sums2是32768次累加量，其绝对差与均值/方差的极小变化一致。这些差值只是两个浮点计算路径的实际差异范围，不是针对真实实数积分的向外误差界。

h=0的raw计数没有精确重现，三维单格最大差271/786/2032；差异全部限于h0的已知plateau，不能将raw当作稳定结论。m2来源周期1/2，故mL为整数时a=1严格成立。回放另强制所有此类结构恒等分支，重新计算完整q与全阈值统计。全部guard计数仍不变，修正后的raw计数也与guard计数完全相同；独立路径修正前后IS均值变化≤3.34×10⁻¹⁶。h>0没有near，raw计数与guard精确一致。这支持原报告采用guard及恒等回放，而不采用h0 raw假阳性。

本复核仍不是严格概率证书：没有向外区间包住kernel、logq和阈值事件，没有连续query格下包，也没有证明PRNG满足理想iid条件。既有EB解释保留这些数值与随机假设；本次一致性不能升级为C_n新值或accepted √n来源的资格证明。

归档精简后最小可移植复核包共7个输入/脚本文件，只依赖Python标准库与NumPy：

- `replay_summary.py`
- `PLAN.json`
- `SUMMARY.json`
- `NEAR_SHA256.json`
- `FROZEN_n32.json`、`FROZEN_n128.json`、`FROZEN_n512.json`

原3个NEAR_INDICES文件合计canonical UTF-8字节数1,406,846，共111,064个near索引。已在NEAR_SHA256冻结3个原文件SHA-256、各文件canonical字节数和索引数。canonical格式为外层dict依次插入`indices`与`scope`，scope精确等于`All near query indices, regenerated from frozen seed; no query/features cache.`，再按`json.dumps(..., indent=2) + '\n'`序列化；dict键插入顺序、空格与末尾换行都参与hash。

新版replay仍实际重建全部nearids。若原near文件存在，则核验原文件SHA并逐项比较索引；若缺失，则按上述完全相同格式序列化重建索引，核对冻结SHA及canonical字节数，不需要生成near数组文件。SHA是对已冻结原证据的计算一致性校验，不是额外数学或随机性证明。

实际无数组可移植测试已执行：把上述7文件复制到临时目录，确认没有任何NEAR_INDICES文件，再运行`python3 replay_summary.py`，三维全部以`canonical_sha_only`通过。实际耗时12.82秒，3个canonical SHA和总字节数全部匹配，均值最大误差仍为7.55e-15。运行后临时目录亦没有生成near数组；程序标准输出仅小摘要。原PLAN/probe/SUMMARY/3份FROZEN的前后SHA全部不变。原near数组没有在本次责任范围内删除，归档时可直接省略。

将这7个文件放到同一目录执行`python3 replay_summary.py`可生成完整REPLAY_SUMMARY。若同目录保留原probe.py，脚本额外核验其PLAN hash；它不import该文件，亦不需要SciPy、RESULT、全query/features缓存或训练重跑。可连同本AUDIT和已生成的REPLAY_SUMMARY交付，便于直接阅读结论。
