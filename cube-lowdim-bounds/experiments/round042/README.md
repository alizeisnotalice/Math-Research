# 第42轮最小复现

仅标准库Fraction，固定4个中心符号证明区域、64条行、60秒硬限。实际全通过；这是一份node663局部证书，不是新一般弱型界。

临时根放round01/scaled_certificate.json，以及第41轮归档的round22/leaf_cover.json和models.local/node_00663.json。将本目录程序放round23，空输出目录执行`python -B run_once.py`。输入SHA锁定，既有结果拒绝重跑。旧模型及私人附件不重复上传。

仅审保存结果时，复制quadrants/result/run/log及上述输入，在没有saved_verification.json的副本运行`python -B verify_saved.py`。生成器用de Casteljau；核验器用独立幂基展开、仿射代入和Bernstein回转，核64行及4块完整分区，不调用生成器算法、Arb或积分。两者共享Python整数/Fraction及旧区间模型，不把这种独立性扩展到旧积分评价。

旧tree depth保持4。pieces是新的代数证书区域，未来cover应保留原节点身份并显式引用本证书，不回写原失败结果。
