# 项目2：加权中心方体 Solyanik 估计

本项目研究内禀权平均、外层加权光晕及共同中心多集合传递。

- `paper/paper.tex` 与 `paper/extension.tex`：当前核验稿，含固定权一阶、全权族精确公式及不等指数水位传递。
- `src/`：必要的 Python、MATLAB、Wolfram 脚本与历史实际运行摘要；有限计算不替代算子证明。
- `rounds/`：逐轮报告。每轮结束后独立提交并推送，不改写既往历史。
- `sources/Pro_Project_2.md`：用户提供的继续研究方向，作为项目输入保存。

本地第一至第五轮在GitHub逐轮推送要求之前完成，没有逐轮Git记录。第六轮建立可复现Git基线，此后按轮提交。

复现：在src目录运行 `python3 exact_certificate.py`、`python3 star_geometry_check.py`、`python3 water_filling.py`；Mathematica运行 `wolframscript -file symbolic_check.wl`；MATLAB运行 `matlab_check.m`。在paper目录用XeLaTeX编译paper.tex两遍，需要宋体与Times New Roman。PDF及56页OCR、渲染和完整本地历史包保留于本机，不重复上传可再生产物。
