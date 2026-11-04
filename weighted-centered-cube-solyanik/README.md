# 项目2：加权中心方体 Solyanik 估计

本项目研究内禀权平均、外层加权光晕及共同中心多集合传递。

- `paper/paper.tex` 与 `paper/extension.tex`：当前核验稿，含固定权一阶、全权族精确公式及不等指数水位传递。
- `src/`：必要的 Python、MATLAB、Wolfram 脚本与历史实际运行摘要；有限计算不替代算子证明。
- `rounds/`：逐轮报告。每轮结束后独立提交并推送，不改写既往历史。
- `sources/Pro_Project_2.md`：用户提供的继续研究方向，作为项目输入保存。

本地第一至第五轮在GitHub逐轮推送要求之前完成，没有逐轮Git记录。第六轮建立可复现Git基线，此后按轮提交。

复现：在src目录运行 `python3 exact_certificate.py`、`python3 star_geometry_check.py`、`python3 water_filling.py`；Mathematica运行 `wolframscript -file symbolic_check.wl`；MATLAB运行 `matlab_check.m`。在paper目录用XeLaTeX编译paper.tex两遍，需要宋体与Times New Roman。PDF及56页OCR、渲染和完整本地历史包保留于本机，不重复上传可再生产物。

第7–10轮已新增显式端点分段、有限质量临界发散、一维无限质量构造和任意维乘积权临界闭合。当前状态以STATE.md最后一节和paper正文为准。

第11轮进一步闭合全部无原子Radon测度的内禀临界问题；当前正文采用统一证明。前几轮的有限质量/乘积权片段留作对应轮次的可复算构造附件，不重复插入正文。

第12轮给出小阈值与线性弱型范数的精确关系及不等指数临界留数；一般权的线性弱型范数本身并未因此获得显式数值。

第15轮闭合与方体几何匹配的径向幂权外层首项，给出一维部分参数范围的第二项，并传回不等指数多集合模型。证明见`paper/power_endpoint.tex`，本机计算证书见`src/power_endpoint/`，审阅和适用范围见本轮报告。

第16轮通过薄壳层析把二阶端点统一到任意维数的明确参数范围，正文见`paper/power_second_order.tex`，并替代重复的一维二阶证明。

第17轮将匹配径向幂权端点推广到任意有限维范数球；通过packing与Aldaz统一弱型定理验证了所需几何条件，无需额外弱型假设。正文见`paper/norm_geometry.tex`；该结果要求权与平均球使用同一范数。
