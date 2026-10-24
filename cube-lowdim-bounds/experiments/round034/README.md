# 第34轮复现与证明范围

本候选未通过全域核验。准确成本较小和局部行通过均不能用于更新C₂。

## 直接核验这个失败候选

在新临时round15解压kernel_inputs.tar.xz（plain候选与含结构的核元数据），解压second_report.json.xz/second.log.xz可查看完整实际报告/日志。候选SHA f46a31bb2b2356651a2e255ad63f058ddd4ce857e14f1365d07d8cabd24de9dd。

用第20轮已归档公共源码于相邻upstream-v4，按full_check_plan两命令可复现；run_stage first/second还要求result.json/exact_cost.json等绑定文件，已有run/log/output先放另外参考目录，wrapper拒绝覆盖。所有stdout与report均为该候选实际运行。两个checker共享python-flint/Arb及解析公式，不是完全独立算术实现。原版first的质量项仅输出；严格成本比较另由整数根证书完成。

## 精确修复与成本复现

完整repair依赖第33轮本地rows.json（约52MB，可按第33轮builder重建）、参数/方向/终态和旧51574点池缓存；代码固定所有SHA，读取相邻round14/round01/round06/round07/round10目录。只有一次λ=1023/1024准确收缩，90秒，0新AS/0LP。不得将缺少旧缓存的环境描述为已经完整复现。

准确成本可以单独用第29轮enclose_cost.py对structured_kernel.json核验：kernel-baseline取第20轮scaled_certificate.json、baseline取其scaling_arithmetic.json；输出另用新文件。该程序是标准库分数和整数40次根，无Arb。exact_cost.json绑定structured原始hash及plain候选hash；不能把该成本证明冒充核正性。

## 精简归档

first_report完整5.6KB保留，second完整报告/日志压缩保留，kernel两种schema压缩在一个包中。没有重复保存第32/33轮大模型/点缓存或私人附件。full_checker_summary提供可读摘要；C1总392取second.log，JSON仅保存50条见证。摘要补录日志总数是只读整理，没有重新运行核验。
