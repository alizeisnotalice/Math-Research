# 第41轮最小复现

实际只处理2新盒，16次积分，第二盒node663在原depth4未接纳即停。输出不完整，不能用于一般常数的新声明。

保留前轮布局：round01基线、第36轮修正源码于round17、round18/round19/round20/round21历史cover、最终模型和补充证书；上游independent hl_core.py/hl_check.py来源固定。所有final和supplemental强制hash；历史partial文件非必需。

本目录源码放新临时根round22，Python3.12/python-flint0.9.0。空输出目录运行`python -B run_once.py`，300秒外部硬限，已存在结果则拒绝重复。每新盒一次AS及中心6Hint+Hab2；8请求不必等于8积分，零延拓分支分别计数。

仅审保存证据：saved_models.tar.xz解压到round22（当前cover和2个final模型），复制result/run/log/calls及来源。`python -B verify_saved.py`只用标准库Fraction核算术、来源、树及分区，尚无saved_verification.json时方可运行。归档包逐字节往返核验，archive_manifest记录SHA。旧模型不重复上传，私人材料不包含。
