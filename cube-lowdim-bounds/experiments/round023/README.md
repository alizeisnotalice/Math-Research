# 第23轮：独立点反例与工具交叉核验

本轮证明第22轮候选不满足核支配及生成元正性；不是C2的新下界或上界，目标3.4尚未实现。原严格上界见第20轮。

## 复现依赖与布局

Python3.12、python-flint0.9.0；符号转换另需SymPy（版本见refinement_sympy.json）。MATLAB R2026a Update4及Symbolic Math用于部分独立计算，记录层发生错误，实际仅1K点与1成本；未重跑原脚本。不得把静态Code Analyzer无问题等同运行成功。

在新本地临时根目录建立round01、round03、round04，复制本目录程序到round04（不要复制既有输出）；从第20轮round01复制scaled_certificate.json到round01，从第22轮复制candidate.json及first_checker_summary.json到round03。将第20轮归档hlcert.py及independent/hl_core.py按同样层级置于round04/public_src。这只复制必要公开依赖，不需要原私人附件。

points.json是固定有理输入，可直接复制；prepare_points.py用于重建坐标，不评价核。`python -B run_once.py`最多48个Arb标量调用、180秒外层截止；`python -B fraction_K.py`另做14个精确分数K点包围；`python -B refinement_sympy.py`检查8段多项式及9端点。程序拒绝覆盖既有输出，本次实际运行无需重复。

`python export_matlab_input.py`重建MATLAB输入，再执行matlab_crosscheck.m会复现本次结构体记录错误。matlab_partial_result.json是从当时工作区提取的已完成首点及成本，不包含未完成的9点。大输入JSON和恢复MAT仅本地保留；不将有缺陷的记录程序当作完成的认证器。后续修复应保存新版本及单独预算。

## 证据独立性

两Arb路径使用不同公开实现，但共享python-flint/Arb与解析公式。fraction_K.py只用Python整数/Fraction重新实现样条并以整数40次根包围幂项，独立于Arb。MATLAB的有理样条独立实现，80位vpa不是向外区间；仅作为交叉检查。SymPy用于精确样条转换恒等式，不能证明任何候选全域正性。

所有严格反例只针对固定第22轮候选，不排除其他核达到3.4。原始调用逐行追踪、MATLAB恢复文件在本地，仓库保存足以核验的最终记录。
