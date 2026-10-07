# 欧氏球中心极大算子研究

状态：维数无关 weak-(1,1) 上界及随维数发散的下界均未在本项目证明。

- [第41轮研究记录](rounds/041/report.md)：来源二阶接口、远层行预算及纯层差衰减反例；[验证摘要](rounds/041/verification.json)。复现：`python3 rounds/041/verify.py`。
- [第40轮研究记录](rounds/040/report.md)：一般避让约束、点态覆盖反例、原 band 正超额有理证书。
- [复现核验](rounds/040/verification.json)；[来源哈希与可移植性修改](rounds/040/provenance.json)。
- `python3 verify_round40.py`：在临时目录从零执行三批28例并核验；`--regenerate` 是兼容别名。

每一轮单独提交并推送；保留连续编号。先前各轮完整材料未在本次导入。AGENT_ROUND_SUBMISSION_PROMPT.md、第三方文献库、缓存和重复打包文件不纳入本项目提交。

文件预算：本轮不超过20个文件，后续每轮默认不超过10个新增文件；复现依赖复用，报告合并，可再生数据不提交。提交前检查具体文件数及体积。GitHub目录宽度推荐不超过3000项、单次diff最多显示300文件；单文件100 MiB硬上限并非整个仓库文件总数上限。参见 [GitHub repository limits](https://docs.github.com/en/repositories/creating-and-managing-repositories/repository-limits)。本项目采用更严格的单文件10 MiB内部预算；超过预算先调整存储方案。
