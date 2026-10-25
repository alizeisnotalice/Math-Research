# 第003轮局部低值诊断：结果

121 个登记点全部通过整数有理节点门控：每点均有 967680 个接受项、零节点和
重复节点计数均为零。实际命令及退出状态为

```sh
c++ -O3 -std=c++20 -Wall -Wextra ../direct_kernel.cpp -o /tmp/cube6_direct_kernel_local
/tmp/cube6_direct_kernel_local 121 local_points.txt > local_direct_results.json
# exit 0; wall 5.087 s
python exact_local_reference.py 24
# exit 0; wall 23.383 s
```

最小 long-double 诊断值在 index 24：

\[
x=(1,998836847/D,892738545/D,807464540/D,467291939/D,96968636/D),
\quad D=1000000007,
\]

其中 \(G_{\rm ld}=2.25013385968857252675\cdot10^{-5}\)。完整 21 cones 与
46080 signed permutations 的独立 `fmpq` 累加及 Arb384 值为

\[
G(x)\in[2.2501338596882661014596145325523546206069246978571478779463719021875719946855087
\cdot10^{-5}\ \pm3.97\cdot10^{-85}].
\]

二者相对差 \(1.3618091611\cdot10^{-13}\)。该点在第一坐标的正方向搜索
端点（第 12 层），故这 121 点恰恰不支持“已达局部最小”这一说法。

这只是一个准确的**单点** Green 函数值和数值路线诊断：它不提供连续极小值、
\(g\) 的下界、\(I\) 的上界、模型成本或 \(C_6\) 的任何上界。搜索未增加
积分样本，且本批到此停止。若继续，必须另行冻结新的有限点集与停止条件。

平台精度审计：本机C++ long double的LDBL_MANT_DIG=53且sizeof=8，实为double精度；本报告中的类型名称不表示80-bit扩展精度。浮点相消/零与形式成本均不是证书，只有明确列出的fmpq/Arb点值另有准确参考。
