# 欧氏球中心极大算子研究

状态：维数无关 weak-(1,1) 上界及随维数发散的下界均未在本项目证明。

- [第46轮研究记录](rounds/046/report.md)：真实前序亏损反例、共同阈值的首次记录区间及平均交叉流；整体无维预算仍未证明。[验证摘要](rounds/046/verification.json)。复现：`python3 rounds/046/verify.py`。
- [第45轮研究记录](rounds/045/report.md)：平均秩亏损的等价性、层内最优重排至多节省 X/6，以及真实高维双原子基准；主无维上界仍未闭合。[验证摘要](rounds/045/verification.json)。复现：`python3 rounds/045/verify.py`。
- [第44轮研究记录](rounds/044/report.md)：连续占用分配、平移后的来源矩与尾收缩；维数损失仍在几何偏移。[验证摘要](rounds/044/verification.json)。复现：`python3 rounds/044/verify.py`。
- [第43轮研究记录](rounds/043/report.md)：全局净流、观察端无维平方预算、密度向量漂移与谱隙障碍。[验证摘要](rounds/043/verification.json)。复现：`python3 rounds/043/verify.py`（标准库有理算术）。
- [第42轮研究记录](rounds/042/report.md)：有符号交换流、远层正流与内半球费用；两个局部质量加强版反例。[验证摘要](rounds/042/verification.json)。复现：`python3 rounds/042/verify.py`（需要 NumPy；高维部分仅为浮点探索）。
- [第41轮研究记录](rounds/041/report.md)：来源二阶接口、远层行预算及纯层差衰减反例；[验证摘要](rounds/041/verification.json)。复现：`python3 rounds/041/verify.py`。
- [第40轮研究记录](rounds/040/report.md)：一般避让约束、点态覆盖反例、原 band 正超额有理证书。
- [复现核验](rounds/040/verification.json)；[来源哈希与可移植性修改](rounds/040/provenance.json)。
- `python3 verify_round40.py`：在临时目录从零执行三批28例并核验；`--regenerate` 是兼容别名。

每一轮单独提交并推送；保留连续编号。先前各轮完整材料未在本次导入。AGENT_ROUND_SUBMISSION_PROMPT.md、第三方文献库、缓存和重复打包文件不纳入本项目提交。

文件预算：本轮不超过20个文件，后续每轮默认不超过10个新增文件；复现依赖复用，报告合并，可再生数据不提交。提交前检查具体文件数及体积。GitHub目录宽度推荐不超过3000项、单次diff最多显示300文件；单文件100 MiB硬上限并非整个仓库文件总数上限。参见 [GitHub repository limits](https://docs.github.com/en/repositories/creating-and-managing-repositories/repository-limits)。本项目采用更严格的单文件10 MiB内部预算；超过预算先调整存储方案。
