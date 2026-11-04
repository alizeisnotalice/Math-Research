# 课题3：分数次极大算子的尺度不变光晕

作者：唐浩天。基线：2026-08-25 原稿，53 页；修订日期：2026-11-05。

当前版本包含原稿全部实质章节的可编辑重写、72 个有编号命题/定义/注/问题的原页映射，以及完整证明。旧修订稿只作为线索重新核查。固定输入极限与全体集合上确界始终分开。

一维精细化已完成：全阶一致 Taylor 余项、有限阈值转换曲线及临界窗口四阶展开。

固定凸体已推广到任意有界满维、含锚点的凸体，得到尺度和平移联合切向极限；另有固定方体奇性权首项。

固定方体的有限接触尺度公式及二维四分支高水平精确面积、三阶系数已完成。

## 当前研究顺序

按用户2026-11-05提供的 [Pro研究指令](source/Pro_Project_3.md) 继续；原文件及SHA-256已归档。
先核验和优化无权单线性的定义、归一化、光晕及全体Tauberian常数，再推进高维、一般输入与高阶匹配；随后处理带权，最后核验多线性sharp信息的传递。
当前新增研究聚焦前两项：第9轮已闭合中心方体全体可测输入的精确首项与紧参数统一余项；全体下一系数及一般几何基仍待研究。

## 文件

- `paper/main.tex`、`paper/main.pdf`：统一中文正文与当前 PDF；章节按引言、预备知识、主要证明、应用组织。
- `paper/references.bib`、`references/verification.md`：24 条核验题录及实际阅读范围。
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

所有一般结论由正文证明支持；有限数值实验、优化搜索和独立实现吻合只作为所列范围的计算证据。中心全体集合的下一系数、近极值分类、一般几何基的全体首项、一般多输入有限阈值等号及 PDE 反向 good-lambda 仍待研究。

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

前两类符号脚本核对系数、精确顶点和积分；MATLAB独立求根、LP和边界积分。结果摘要保存在 `results/`。各轮具体成果见 `rounds/001.md` 至 `rounds/009.md`。

第5轮确定中心对称、β=n/2的固定模型首项与椭球极值。复现命令为 `wolframscript -file verification/half_order.wl`；MATLAB运行 `verification/half_order_matlab.m`。优化范围是资料固定E=K且基与K共同变化；一般凸体基对全体E的上确界仍为后续问题；中心方体基的首项已在第9轮确定。

第6轮完成二维中间临界点的两层一致过渡：β=1+λt 的四阶复合公式及 β=1−t/2+μt² 的分段四阶系数，余项分别一致为 O_Λ(t⁵)、O_M(t⁵)。运行 `wolframscript -file verification/critical_square.wl`，MATLAB运行 `verification/critical_square_matlab.m`，随后 `python3 verification/critical_square.py --workers 8` 可复现符号、独立边界积分与跨工具比较。

第7轮在二维固定方体上建立任意阶统一边界截断，给出显式适用阈值和误差常数，允许β(t)取紧区间内任意路径；完成β=1/2、3/2附近的显式三阶过渡。任意维方体外部中心的精确顶层为max(1/2,2^(β−n))。复现：`wolframscript -file verification/uniform_square.wl`；MATLAB运行 `verification/uniform_square_matlab.m`；`python3 verification/uniform_square.py`。138组N=1至6积分近似、80组侧窗口及60组MATLAB交叉检验已实际运行。

第7轮PDF为46页。第7轮还证明：固定β∈(3/2,2)时，N阶截断的首个误差系数严格为正，t^(N+2)阶一般不可提高。

第8轮把固定方体推进到任意维数：2n个唯一径向接触根、带显式Bell递推余项的任意阶统一表示、全部临界参数的显式次首项T。首项后的t^(n+1)阶在每个固定β处都不可提高；任意N阶径向截断的t^(n+N)误差阶在外扩区也严格sharp。运行 `wolframscript -file verification/higher_dimensional.wl`，MATLAB运行 `verification/higher_dimensional_matlab.m`，再运行 `python3 verification/higher_dimensional.py`。最新文献核验见 `references/search-log-20261105.md`。

第8轮PDF为49页，其当轮范围见 `rounds/008.md`。

第9轮的 gl:main 证明中心方体全体常数 C(1−t)=L t^n+O_J(t^(n+1))，L=2^n/[floor(2c)!(2c)^(n−floor(2c))]、c=min(β,n−β)。共同核质量修改、任意可测面层容量上界和Brunn–Minkowski凸纤维平均共同提供匹配上界，允许输入随t变化。正次首项区还给出全体余项的sharp阶；中心限制型高矩首项等于n!L。复现 `wolframscript -file verification/global_leading.wl`，MATLAB执行 `verification/global_leading_matlab.m`，再运行 `python3 verification/global_leading.py`。数值核对不替代一般证明。

第9轮当前PDF为54页；正文、逐项审计、独立复审及最小复现材料同步归档，详见 `rounds/009.md`。

第10轮在n=2、β=1证明矩形输入类的精确高水平上确界R(α)=(α⁻¹−α)²/2，严格超过固定正方形测试值；精确差为(1−α)⁴/(2α²)。给出任意阶精确余项，但尚未证明全体C=R或全体下一系数。复现 `wolframscript -file verification/rectangle.wl`、MATLAB执行 `verification/rectangle_matlab.m`、`python3 verification/rectangle.py`；当前PDF55页，见 `rounds/010.md`。
