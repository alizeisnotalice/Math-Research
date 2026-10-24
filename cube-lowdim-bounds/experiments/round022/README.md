# 第三优化轮：26维有限方向，数值候选与严格核验分列

候选candidate.json不是已证比较核；查看认证终态报告再使用。准确成本较小本身不推出弱型上界。原已证结果见第20轮。

## 最小重现

依赖Python3.12、numpy2.5.4、scipy1.18.1、highspy1.15.1；全域核验另需python-flint0.9.0。

在新临时目录复制search.py及run_once.py（不复制既有result或discovery日志）；将第21轮search.py复制为dependency_round02_search.py，并复制该轮formula_modules子目录；将第20轮round01/scaled_certificate.json和scaling_arithmetic.json分别复制为input_baseline.json、input_scaling_arithmetic.json。运行`python -B run_once.py`。这只是重现说明，当前唯一实际运行已完成，不需重跑。

准确成本：`python enclose_cost.py candidate.json --baseline input_scaling_arithmetic.json --out new_cost.json`。此程序仅整数/Fraction成本核验，不检查正性。

完整核验复用第20轮已归档的两套源码及依赖，将候选路径换为本目录candidate.json，其他参数和全域条件保持；本轮实际命令及终态见CERTIFICATION_SUMMARY.json。两个核验器共享Arb/python-flint及数学公式，独立实现不等于完全独立可信基。

归档不重复拷贝公式模块、公开核验源码、基线证书和许可证，均引用第20、21轮；不保存可再生矩阵缓存。唯一候选保留是为了重现失败路径，不作新上界证书宣传。
