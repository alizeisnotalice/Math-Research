# 第24轮：全局粗细样条有限发现的接口阻断

此轮没有生成候选或新弱型上界。程序语义终态是`failed_or_rejected`，不能把外包装exit0当成搜索通过。准确基准仍为第20轮B1<3.59091，尚未证明C₂<3.4。

实际一次运行2.825331417秒：K26604/A24970，共51574个已证基线点残差全部通过；192方向共9902208个列条目完成；LP调用为0。在第一次addRows返回非kOk时停止，原日志没有保存具体HighsStatus，不能断言错误还是警告。没有last_model.mat，因此MATLAB脚本只静态检查、未求解。MATLAB已实际确认interior-point不使用MaxTime，故计划采用支持MaxTime的dual-simplex-highs，与Python共享HiGHS内核。

四倍三次B样条细化恒等式已独立用SymPy实际精确核验16开区间和17端点；这证明核族展开与质量公式的适用性，不证明任何新核满足全域正性。完整数值基线通过仍是有限点结果，不替代第20轮严格全域证书。

## 复现

依赖Python3.12、numpy2.5.4、scipy1.18.1、highspy1.15.1、SymPy（版本见refinement_identity.json）。在新临时目录复制search.py/run_once.py/refinement_identity.json；第21轮search.py复制为dependency_round02_search.py，同轮formula_modules放在该目录；第20轮round01/scaled_certificate.json复制为input_baseline.json；第23轮points.json复制为input_diagnostic_points.json。执行`python -B run_once.py`。已有真实运行无需重跑。

若复现SymPy恒等式，使用refinement_identity.py；原16段17端点结果已经完整保存，本轮不重复执行。matlab_crosscheck.m是本轮未运行求解的预备脚本，依赖成功导出的last_model.mat；本轮不存在该输入，不可把脚本存在视为实际MATLAB求解。

全量可再生点表（2.3MiB）保留本地，不重复归档；矩阵缓存未持久化，下一轮为定位失败须有理由地重建一次并在API调用之前保存。源码和运行记录保持原样，不在原轮修改或重新运行。公式模块/许可证/已证基线引用第20、21轮，避免重复存储。
