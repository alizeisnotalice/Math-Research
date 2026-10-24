# 第28轮最小复现

**没有新弱型上界。** 本轮严格证明第27轮候选有负核值，并实际完成10个SymPy恒等式、18个不同表达的生成元积分（23106回调）。数值积分不作为全域证书。

## 复现

使用数学Python：Python3.12、SymPy/SciPy/mpmath，准确版本见result.json。恢复scratch目录round08、round06、round09：round08/candidate.json取前轮归档，round06/point_manifest.json按第25轮恢复且hash必须为513c8d953695d9ddd1c7a097ee9af2ab624f79b875330d02e0d42ad2dcc70ebd。round09放diagnose.py/run_once.py，运行 `python -B run_once.py`；180秒硬限，存在运行记录即拒绝重复。本次已完成，无需为查看结果重新执行。

独立小证书negative_point_certificate.json只需16个准确系数与权重向量(8,93,60,1)/162即可检查负值。它来自同一次Fraction求和结果；没有宣称另一个程序重新计算了该点。原候选完整hash与先前归档绑定，支持截断部分准确为零；连续性把负点扩展到开邻域。

新族公式与分布桥在第27轮next_kernel_family.md/theory_review.md，本轮仅实测有限公式与归一化。源码静态审读、实际结果审计独立完成，未重算积分。两数值路线分别SciPy原始对称差和mpmath二阶导数含负delta；beta=alpha用Gamma式作控制，beta>alpha检查与凸性参考相容，均不称严格包含。

## 精简记录

result.json保存全部值/分段/回调/误差估计。calls_ledger.json从实际36事件trace精简抽取状态和计数，完整trace的hash固定且原文件仅本地，避免再次归档重复的完成行数值。point_selection.json是旧数组选点来源，不是本轮新函数评价。代码、证书、实际退出/日志和审查文件归档；无新增大矩阵或私人原件。
