# 第31轮最小复现

没有新C2界。预定构建在24个基线模型后早停，0个T模型、0次LP、无候选。不能把1536项上限写成实际次数。

代码依赖Tomar v4 checker/independent/hl_core.py和hl_check.py（commit bca281a706bbf6375b04d601ad19ba4ad0bf28b6；代码内固定SHA）。放于../upstream-v4/checker/independent或本目录public_src。已证核取第20轮scaled_certificate.json，放../round01；第29轮result.json放../round10，均由程序核hash。python-flint0.9.0，Python具体版本见build_result.json。运行 `python -B run_stage.py build`，300秒hard；已存在输出即拒绝重复。

本次baseline失败，程序未进入T方向分支，所以无需执行/恢复该分支的旧refinement和方向依赖。不要启动search；本轮也未实现其求解器。备用成本程序未运行且未归档，避免把准备工作误作计算成果。

cells.json保存12盒、192个准确基线行和24模型hash。build_result.json保存24模型实际计数及失败门槛。model_manifest登记全部模型；只归档最坏box1的ψ/S1两模型以核查明确的负充分行，其余22模型本地保留、可从已归档确定性builder恢复。完整模型不影响已有K1证书，本轮不是重证它，也不推翻它。

两个审读文件职责分开：source_review/independent_result_audit为独立源码及终态核账；diagnostic_breakdown是从已保存Q端点/余项做准确算术分解，无新函数评价。只保留最终分工稿，未归档并发暂名result_review.md。linear_cell_lemma给LP安全行证明。

目标仍为一般L1和全部尺度C2<3.4。本轮仅建立连续约束接口、发现粗包围门槛，不构成最优性/不可行性定理。后续需新的预定范围和预算，不在此结果后自动细分或扩算。
