# R074：三维edge局部容量与真实inside合同

[最终v3](../../proofs/074-edge-inside.md)和[非作者数学重建](../../reviews/074-edge-inside.md)证明局部全点容量常数一，并在实际同N9合同下给inside领先2。它不确定C₃九阶系数，也不自动移到未重建的新递归层。

`python3 finite_checks.py`是作者实际Fraction核，2016有理点/常数/量纲仅作辅助；`python3 edge_regions.py > 新输出.json`是独审实际24主组、276对仿射包络的连续参数有限分类核。其原输出保留，常数和几何的真正证明见完整源。

独审发现原q合同未控制已删除原源D0，v3明确q=旧raw+实际扩张、Y=whole-normal零次名义上包。另v2排除了原本误留的九阶边界候选（实为v¹¹）。两份原始JSON补丁保全部旧新文字，`python3 contract_checks.py`从公开v3在内存逆构造v2/v1并核所有SHA，随后执行独审原有的有理单项式/删除导数断言。为独立复现，公开脚本仅改文件发现/逆构造前缀，数学断言逐字未改；原脚本/输出SHA见`adaptation.json`，原本机版本保留在仓库外。公开适配已实际再跑通过，输出为此处contract_checks.json。

所有脚本仅Python标准库，没有区间算术或浮点认证声明。补丁中机器旧路径是来源记录，不作运行依赖。`manifest.json`绑定冻结源、独审和实际证据。
