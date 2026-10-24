# 第29轮复现与证据范围

没有新弱型界。当前已证B1<3.59091；候选成本3.590691293368599223...虽严格下降，两全域核验均在生成元条件失败。

`candidate.json`含核参数结构s/T及完整核系数；`exact_cost.json`为独立整数根/Fraction成本和结构结果。原核基线取第20轮scaled_certificate.json，对应scratch round01；成本基线scaling_arithmetic.json同轮，程序固定hash。运行：`python enclose_cost.py candidate.json --baseline ../round01/scaling_arithmetic.json --kernel-baseline ../round01/scaled_certificate.json --out new_cost.json`。现有成本已计算，不为查看结果重复运行。

LP恢复scratch round06/07缓存，按第25/26轮的来源、point_manifest及hash加载；使用本轮search.py/batch_solver.m/run_once.py，Python3.12+NumPy/SciPy与本机MATLAB R2026a Optimization Toolbox。唯一`python -B run_once.py`最多180秒，已有运行记录拒绝重跑。三个模型都51574行128列；旧B/L缓存和大模型仅本地，local_artifact_manifest.json登记其本轮引用与产物。result.json和三份MATLAB结果保留实际解/残差/状态。MATLAB与Python HiGHS不是两个独立优化器。

两全域checker来自Tomar v4 commit bca281a706bbf6375b04d601ad19ba4ad0bf28b6。对candidate只保留alpha/R/h/N/c/a六字段输出candidate_public_schema.json（beta准确等于alpha）；schema_binding.json固定新旧hash，不重复存相同1181系数。按CERTIFICATION_PLAN与*.run.json命令各运行一次，各600秒限。run_bounded.py复用已有公开核验包装。二者均exit1，不应再次调深或重复。两实现共享Arb/解析公式，成本整数根实现另立。

完整失败叶JSON约8.6MiB和574KiB仅本地；full_checker_summary.json保留来源hash、计数和少量失败盒。第二checker大量逐单元日志不上传，运行记录保存完整日志hash。负的充分下界不是严格负点。解析附件中Bregman正余项是下一机制，没有新积分或数值收益证据。
