# 第32轮复现与证据范围

本轮只构造三个父格的基线连续约束，不是新全域证书或新C₂界。完整范围、预算见预运行PLAN.md；该文件保留预运行时态，实际结果以result.json和两个run.json为准。

## 先核查保存证据（无Arb）

在新临时目录创建round13，将本目录文件复制进去，解压accepted_cover.tar.xz（包含leaf_cover.json及1017个接受叶原始模型）。把第20轮round01/scaled_certificate.json复制为临时目录的round01/scaled_certificate.json。在round13移走已有cover_verification.json、verify.run.json、verify.log后运行 `python -B run_verify_once.py`。程序拒绝覆盖已有证据，使用标准库Fraction，90秒硬预算。原始本机记录不要删除或重跑。

核验器检查全部模型原始SHA、准确行、参数及几何映射，并用有理x切片独立检查无缝覆盖，不依赖导出的树结构。它不独立证明Arb端点正确，也不调用超越函数。物理覆盖面积仅3/1024。

## 从模型公式重建（需要python-flint0.9.0）

在另一新目录布置同样round01；将historical_scaled_second.json.xz解压为round13的上一级scaled_second.json。将第20轮upstream-v4/checker/independent的hl_core.py及hl_check.py放round13/public_src（固定版本与SHA在代码中）。只复制export_cover.py/run_once.py，执行 `python -B run_once.py`。300秒、1500盒、深度4固定。该运行与已存模型核查是两种复现路线，日常不重复已有计算。

accepted_cover.tar.xz保存独立准确核查所必需的完整1017接受模型及行；15个被细分父盒的模型留在本地，没有重复上传。数值字符串仅诊断，核验使用exact man_exp导出的分数。归档未包含私人原件。完整baseline已在第20轮核验，本轮没有重复该核验。

首次修改verify路径失误只产生FileNotFoundError，未修改文件也未运行verify；随后在首次verify前修正Taylor模式大小写。导出与核验均只有一次实际运行，两个源码前后SHA一致。
