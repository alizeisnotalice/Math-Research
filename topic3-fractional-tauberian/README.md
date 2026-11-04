# 课题3：分数次极大算子的尺度不变光晕

作者：唐浩天。基线：2026-08-25 原稿，53 页；修订日期：2026-11-04。

当前版本包含原稿全部实质章节的可编辑重写、72 个有编号命题/定义/注/问题的原页映射，以及完整证明。旧修订稿只作为线索重新核查。固定输入极限与全体集合上确界始终分开。

一维精细化已完成：全阶一致 Taylor 余项、有限阈值转换曲线及临界窗口四阶展开。

固定凸体已推广到任意有界满维、含锚点的凸体，得到尺度和平移联合切向极限；另有固定方体奇性权首项。

固定方体的有限接触尺度公式及二维四分支高水平精确面积、三阶系数已完成。

## 文件

- `paper/main.tex`、`paper/main.pdf`：统一中文正文与当前 PDF；章节按引言、预备知识、主要证明、应用组织。
- `paper/references.bib`、`references/verification.md`：19 条核验题录及实际阅读范围。
- `audit/`：逐项原编号、原页、假设、证据、状态和新标签；`source/` 保存原 PDF 与 SHA-256。
- `verification/`、`results/`：必要验证代码与实际运行摘要。
- `rounds/`：每轮结果、未解决问题和复现信息。代码与数值结果的作用以各报告为准。

## 复现

在本目录运行：

```sh
python3 -m pip install -r requirements.txt
python3 verification/check_manuscript.py
python3 verification/baseline.py
cd paper
latexmk -xelatex -interaction=nonstopmode -halt-on-error main.tex
```

LaTeX 使用 XeLaTeX、BibTeX，以及 macOS Songti SC / Heiti SC / PingFang SC 字体；异机可在主文导言区替换为已安装的中文字体。运行环境实测为 Python 3.12.13、SymPy 1.14.0、NumPy 2.5.4、SciPy 1.18.1、mpmath 1.3.0，MATLAB R2026a，Wolfram 15.0。

所有一般结论由正文证明支持；有限数值实验、优化搜索和独立实现吻合只作为所列范围的计算证据。高维全体集合的中心最优首项、近极值分类、一般多输入有限阈值等号及 PDE 反向 good-lambda 仍待研究。

## Mathematica / MATLAB 复现

```sh
wolframscript -file verification/one_dimensional.wl
wolframscript -file verification/polytope_tangent.wl
wolframscript -file verification/fixed_cube.wl
python3 verification/one_dimensional.py
python3 verification/fixed_cube.py --workers 12
```

在MATLAB中，从项目目录分别运行：

```matlab
run('verification/one_dimensional_matlab.m')
run('verification/polytope_tangent.m')
run('verification/fixed_cube_matlab.m')
```

前两类符号脚本核对系数、精确顶点和积分；MATLAB独立求根、LP和边界积分。结果摘要保存在 `results/`。最新PDF为41页，各轮具体成果见 `rounds/001.md` 至 `rounds/004.md`。
