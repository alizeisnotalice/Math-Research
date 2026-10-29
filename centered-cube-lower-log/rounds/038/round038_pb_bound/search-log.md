# 定向原始定理检索记录

2026-10-29，按 `academic-literature-search` 技能进行窄问题检索。目标是独立非同分布 Berry–Esseen 的一个明确可用常数，不是研究综述，不以检索数量代替适用条件。以下均为实际web检索/阅读。

| 查询或动作 | 结果与采用范围 |
|---|---|
| `Tyurin non identically distributed Berry Esseen constant 0.56 theorem original paper pdf` | 定位原始arXiv论文、MathNet出版记录及其他转述；不采用百科或学位论文作常数依据。 |
| `arxiv Tyurin Lyapunov theorem 0.5606 0.56 independent` | 定位arXiv0912.0726。其摘要明确区分一般0.5606与iid0.4785。 |
| 打开 https://arxiv.org/abs/0912.0726 | 核验题名、作者、2009-12-03 v1提交日期、19页、arXiv DOI。PDF内部打印日期2018不作为提交年份。 |
| 打开 https://arxiv.org/pdf/0912.0726 | 实际阅读§1定义与式(1)、第6页定理7、其第15–17页证明/数值常数说明。已核对“独立但无需同分布”及三阶绝对矩定义。不是只读摘要。未重做其全文所有中间定理或计算机辅助常数认证。 |
| MathNet英文PDF直链（rm9337_eng.pdf） | 实际打开重定向至MathNet主页，未取得该PDF；故不把该正式发表全文算作已读，也不据不完整出版字段写BibTeX。 |

唯一用于定理引用的已读原始源是上述Tyurin v1。它足以支持本轮选用C=0.5606（或更松C=1）；未主张这是当前最佳常数。检索中出现其他较小数值的结果，未核读其完整条件，故不采用。直接将0.5606向下舍入0.56不受本轮已读来源支持。`references.bib`仅含已核验的原始论文元数据。

证据区分：引用文献支撑一般BE定理；本轮57组NumPy DP仅检查有限输入；本轮Fraction小例仅认证其有限代数输入；本报告中的分箱、宏观p范围、正权重均值预算为从已列假设逐式推导。
