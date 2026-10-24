# 第38轮：单盒多项式下界比较

只对第37轮node615作一次新Psi模型，复用原AS端点。实际0.1427655秒、1模型/12Hint/12真实积分均完成；16准确下界全>7/100。它不是父52,1完整覆盖，也不是新一般弱型界。

复现布局：新临时根建round01/round17/round18/round19。第20轮scaled_certificate.json放round01；用第36轮make_patch.py及已锁定原件在round17生成fixed_checker/hlcert.py；第37轮saved_models.tar.xz解压至round18（至少需原leaf_cover.json和models.local/node_00615.json）；本目录compare_once.py/run_once.py放round19。Python3.12、既有python-flint0.9.0；各输入SHA由源码锁定。

空输出目录运行 `python -B run_once.py`，60秒硬截止、192bit、只1次模型和至多12积分。已有result/new_model/calls或运行元数据会拒绝重复。终态后 `python -B verify_saved.py` 仅以标准库Fraction重组16行并核50个调用事件，不调用积分或Arb。

只审既有算术可复制result/new_model/calls/run及对应源码和旧模型后运行verify_saved.py，输出saved_verification.json必须尚不存在。原模型/cover/修正版源全部引用前轮，不重复上传。代理统计严格转交原函数和参数；记录1模型和12积分的区别。

端点来自Arb精确lower/upper的man_exp，不能以显示小数替代。新Psi模型返回已升阶4×4 Bernstein且未乘c；16条组合只乘c1一次，AS直接复用准确端点。新旧同基底，15条改善、1条变差，并非新模型逐项支配旧模型；但最小行从负变为严格正。

proofs/exterior_integral_taylor.md与analytic_audit.md给出积分余项和C1零延拓论证。首次数值实验不跨坐标R；不可把解析适用范围当作已实际运行的多盒覆盖。所有未来扩展需另立预算。
