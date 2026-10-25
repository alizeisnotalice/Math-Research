# 第005轮：\(H_{1/1000}\) 有限 Sobol 积分路线诊断

固定四个 scrambled Sobol 5-face 组（每组 128 个有理量化点），直接计算完整
21-fan × 46080 signed-permutation 群的截断 Hermite 势 \(H_{1/1000}\)。C++
节点相等性为精确有理判定，但本机 `long double` 是 53-bit、8-byte binary64，
故 512 个 C++ 值只作数值探路。

先验 8 点检查全部为正（取消比约 \(4.4\cdot10^4\) 至
\(5.8\cdot10^4\)）；随后完整 512 点也全部为正。总逐点 CPU 计时为 38.372 秒，
单点最多 0.08395 秒，且未超过 20 分钟批次限制。

令 \(q=6/5\)。固定样本的**形式**平均给出

\[
\overline{H^q}=8.18422230801594\cdot10^{-6},\qquad
(\overline{H^q})^{1/q}=5.76520181480976\cdot10^{-5}.
\]

作为路线目标，\(1/16000=6.25\cdot10^{-5}\)，其相应 \(I_H\) 门槛是
\((1/16000)^{6/5}=9.01687441192008\cdot10^{-6}\)。形式样本在这两个数值
门槛之内，但这绝非积分误差认证或上界。

四组形式范数分别为 \(5.81285,5.80028,5.81496,5.63230\) 乘 \(10^{-5}\)。
最大取消比 \(1.4722\cdot10^5\)，最近 \(S_2\) 诊断距离为 0.065442；最大贡献点
距 \(S_2\) 为 0.133029。前 1%（6 个）贡献了形式总量的 5.1203%，因此没有将
低维 32 点样本当作积分代理。

最大贡献点 index 90 的 C++ 值为 \(2.65642527976120\cdot10^{-4}\)。对同一个
点以 Arb384 重算完整 967680 项，得到

\[
H_{1/1000}(x)=[2.6564252797057752541730863832535935278151977295950896850382916632189202261274537004
\cdot10^{-4}\pm1.96\cdot10^{-100}],
\]

与 binary64 探路值相对差 \(2.08656\cdot10^{-11}<10^{-6}\)。这是一个点的实现
校准，不验证 512 点的统一浮点误差，更不验证 Sobol 求积误差。

本批不含 \(G-H_\epsilon\) 尾、全局 \(g\) 下界、连续积分封闭或一般 \(L^1\)
传递。因此没有得到 \(I\)、上界成本或 \(C_6\) 的严格新界。

最小复现：

```sh
python generate_sobol.py
c++ -O3 -std=c++20 -Wall -Wextra probe_kernel.cpp -o /tmp/cube6_integral_probe
/tmp/cube6_integral_probe ../../round003/direct_kernel/fan_integer_input.txt sobol_pilot8.txt > pilot8_results.json
/tmp/cube6_integral_probe ../../round003/direct_kernel/fan_integer_input.txt sobol_4x128.txt > full512_results.json
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python analyze_probe.py
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python arb_max_reference.py 90
```
