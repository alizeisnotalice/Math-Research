# 第25轮：装载审计通过，Python数值错误与MATLAB有限解

本轮不更新一般L¹弱型界，仍为第20轮B1<3.59091；尚未证明C₂<3.4。

唯一Python发现运行4.726954333秒，基线51574点及192方向完成；实际addRows返回kWarning，3329项≤1e−9被忽略。载入后逐项矩阵与所有界/成本核对通过；仅这些删除项不同。Fraction独立核对最大scaled行盒误差约4.441828812e−9。删小项会改变自由变量LP，不能把它称为原模型的准确等价变换。

唯一Python LP调用约0.96秒返回kError/Not Set、无有效primal，原日志明确dual simplex ratio test失败；非超时/非不可行。随后唯一MATLAB有限模型求解约3.17秒，exitflag1，成本3.590583140728088，30968行的最低slack为−8.01069668e−10。这是有限矩阵数值解，未评价其全固定池核残差、未理性化、未运行全域核验。MATLAB与Python共享HiGHS内核，不称完全独立solver。

## 最小复现

依赖Python3.12、numpy2.5.4、scipy1.18.1、highspy1.15.1；MATLAB R2026a Update4及Optimization Toolbox。仅在新临时目录复制search.py/run_once.py。第21轮search.py复制为dependency_round02_search.py，同轮formula_modules放在此目录；第20轮round01/scaled_certificate.json复制为input_baseline.json；第23轮points.json复制为input_diagnostic_points.json；第24轮refinement_identity.json复制同名。

还需input_previous_point_manifest.json：原始2.3MiB表本地保存。可仅调用第24轮search.py的pools(input_diagnostic_points路径)，取返回keys，使用其dump写`{kind:[[str(u),str(v)] ...]}`生成同字节表；这仅构造准确点坐标，不执行核公式或LP。表的SHA256须吻合result.json依赖字段。之后`python -B run_once.py`。现有唯一运行不需重复。

MATLAB脚本直接复用第24轮matlab_crosscheck.m；复制到新运行目录，在last_model.mat产生后执行一次。实际源/输入/输出hash见matlab.run.json；它记录的是30秒求解器软限，不声称外部硬限。该解不能自动作为本轮核候选输出。

interface_probe.py是前置固定3次装载观察，10秒限、0LP/0核点。normalized结果只将Infinity界转为字符串，其raw源hash保留；不倒推上一轮未保存的实际返回。

## 归档边界

保留必要源码、计划、API/运行日志、3329删除项的准确证据及审计、MATLAB数值解和hash。约170MB B/L/方向缓存、稀疏输入/加载矩阵、各MAT快照和可再生坐标仅本地留存；哈希与最小复现说明保留，不上传私人原件。下一轮可复用这些缓存和有限解，避免重复核公式构建。
