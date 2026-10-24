# 第33轮：最小拒绝证据与复现

本轮取得完整12方向连续模型，运行一次MATLAB，但候选未通过预定half-margin门槛；没有全域checker或新常数。保存的是失败路线证据，不能把MATLAB的最优状态当有限LP严格最优性证明。

## 最小独立算术核查

将verify_rejection.py和rejection_witness.json复制到新的空目录，执行 `python -B verify_rejection.py`。仅标准库Fraction，不依赖Arb/NumPy/MATLAB。核对12个已保存方向模型的对应Bernstein端点、midpoint、误差penalty和该单行准确参数，得到half-margin残差<0而原安全行>0。原本地输出拒绝覆盖；若核验归档副本，将新输出置于临时目录，不覆盖原运行记录。

## 完整模型与LP复现

源代码按原scratch邻目录组织：round01为第20轮已证基线；round05为第24轮refinement_identity；round06为第25轮search.py及缓存构建结果；round07为第26轮288列完整缓存；round10为第29轮result/浮点参数；round13为第32轮解压后的完整接受cover和核验记录。公共hl_core/hl_check取第20轮归档，固定commit和SHA由代码检查。各原始目录可按旧实验README重建；旧缓存是数值发现输入，不是证书，不应重新命名为严格端点。

在无输出的新round14内运行数学Python `python -B run_stage.py build`；成功后 `python -B run_stage.py solve`。build至多12204 AS模型/300秒；solve至多一次MATLAB dual-simplex-highs20秒soft/180秒hard。依赖Python3.12、python-flint0.9.0、NumPy/SciPy及MATLAB R2026a Optimization Toolbox。只调用新启动的batch，不改变GUI。

实际完整rows.json约52MB、原始12方向JSONL约147MB和MATLAB矩阵均在本地保留，不上传；它们可由最终builder重建。directions.json保留每文件SHA，build_result保留完整rows SHA；本归档保留拒绝门槛所必需的单行和12原模型，不声称仅凭它独立复核全部16272行。完整计数由原始运行记录支持；单行诊断与全域证明明确分开。

证据提取曾一次用相对'.'的parent导致找不到round13，在任何诊断输出生成前修正为resolve；未执行新模型/LP，也没有覆盖已运行证据。其后一次单行Fraction诊断通过。本轮不调整舍入、不混合基线、不重试solver。
