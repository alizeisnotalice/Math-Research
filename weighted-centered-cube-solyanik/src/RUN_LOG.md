# 实际运行回执

运行日期：2026-11-04（本机 Asia/Shanghai）。以下为已完成命令的压缩回执；完整的可重复输出由同目录脚本重新生成。

| 检查 | 环境 | 结果 |
|---|---|---|
| `exact_certificate.py` | Python 3.15.0 | `m=2,3,4,8` 均有 `s(lambda_star)=1, Phi=0, PASS`；九个 `n=1,2,3` 和 `alpha=1/2,3/4,9/10` 的星型有理构造均为 `PASS`。 |
| `scan_independent.py` | Python 3.15.0 | 24 个临界两侧样本 (`m=2,3,4,8`，三种 `epsilon`) 全为 `PASS`；45 个近端点星型浮点样本均严格低于目标值并严格超阈。 |
| `symbolic_check.wl` | Wolfram 15.0.0 for Mac OS X ARM | 返回四个精确阈值恒等式和三阶级数；星型严格不等式返回 `True`。 |
| `matlab_check.m` | MATLAB R2026a Update 4 | 24 个临界两侧检查与 45 个星型样本都输出 `pass=1`。 |
| `star_geometry_check.py` | Python 3.15.0 | `n=1,2,3, delta=1/10` 分别输出 2、4、8 个卫星 `PASS`。 |
| `water_filling.py` | Python 3.15.0 | 四个权向量案例中，粗有理网格值均不超过精确 KKT 值。 |
| `py_compile` 与 MATLAB Code Analyzer | Python / MATLAB | 所有 Python 文件编译成功；MATLAB Code Analyzer 返回空问题列表。 |

数值回执只确认脚本所编码的有限命题；完整的数学含义及限制见 `report.md`。
