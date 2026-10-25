# 直接谱核的 32 点路线诊断

对每个 generic 点，计算

\[
G(x)=\frac1{128\pi^5}\sum_R\alpha_R
\operatorname{DD}_5(\operatorname{sgn};b_{R,1},\ldots,b_{R,6}),
\]

其中每个 point 有 21 个 base simplex 与完整 46080 个 signed-permutation
像，共 967680 个项。节点 zero/repeated 的判定使用整数有理式，不依赖浮点。

大素数分母的 32 点全部为 `GENERIC_OK`。实际总计时 1.30604 秒，单点最多
0.04566 秒；最大绝对项/有符号和 cancellation ratio 为约 \(1.40481\cdot10^6\)。
首点的完全独立 `fmpq` 参考值是

\[
[6.3984213567624446970612568299985983373603382701502085420992921376686675008584333
\cdot10^{-5}\ \pm1.50\cdot10^{-85}],
\]

而 C++ 值为 \(6.39842135673982659\cdot10^{-5}\)，相对差
\(3.535\cdot10^{-12}\)。这只校准本次 generic direct-sum 实现，不能转化为
连续误差界。

令 \(q=6/5\)，此固定 32 点的算术样本给出

\[
\operatorname{mean}(G^q)=7.73539157547646\cdot10^{-6},\quad
\min G=2.29461160356569\cdot10^{-5},
\]

和形式样本成本

\[
6^{6/5}\,\operatorname{mean}(G^{6/5})/(\min G)^{6/5}
\approx24.5139209381.
\]

这不是 \(I\) 的上界、不是 \(g\) 的下界，也不是 \(C_6\) 上界。点集排除了
face 边界，且没有 QMC/积分/极小值误差预算；最近 \(S_2\) 距离约 0.14387，
样本最小点距离约 0.34251，均只是避免已知奇异层的诊断。

最小复现命令：

```sh
python generate_integer_input.py --denominator 1000000007 --seed 6003004 --points-output face_points_generic.txt
c++ -O3 -std=c++20 direct_kernel.cpp -o /tmp/cube6_direct_kernel
/tmp/cube6_direct_kernel 32 face_points_generic.txt > direct_kernel_generic_results.json
python exact_generic_reference.py
python analyze_direct_diagnostic.py
```

后续若推进，只能以连续有理/区间 certificate 支付从有限点到全 face、全 cube 的缺口。

平台精度审计：本机C++ long double的LDBL_MANT_DIG=53且sizeof=8，实为double精度；本报告中的类型名称不表示80-bit扩展精度。浮点相消/零与形式成本均不是证书，只有明确列出的fmpq/Arb点值另有准确参考。
