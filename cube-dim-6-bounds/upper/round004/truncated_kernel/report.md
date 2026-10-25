# 第004轮：截断核 Hermite 诊断报告

按 theory round-002 的公式，在完整 21 个 base cones 和每个 cone 的 46080
个 signed-permutation 像上计算

\[
 H_\epsilon(x)=(2\pi)^{-6}\sum_R\alpha_R
 \operatorname{DD}_5[\arctan(\cdot/\epsilon);b_{R,1},\ldots,b_{R,6}].
\]

节点排序和相等性通过整数有理交叉乘法完成。每个重节点区块使用 confluent
Newton 表以及
\(F_\epsilon^{(k)}(b)/k!\)；没有跳过 signed 项或重节点项。

九个 C++ 探路值如下（本机 `long double` 为 `LDBL_MANT_DIG=53`、
`sizeof=8` 的 binary64，不是证书）：

|代表点|\(\epsilon\)|\(H_\epsilon\) binary64|相邻重节点事件|取消比|
|---|---:|---:|---:|---:|
|minimum-near|1/100|2.18121930038250e-5|0|7.094e2|
|minimum-near|1/200|2.24445349483255e-5|0|4.254e3|
|minimum-near|1/500|2.25004858573742e-5|0|2.160e4|
|two-level \(S_2\)|1/100|1.41544993281379e-4|684720|1.429e2|
|two-level \(S_2\)|1/200|2.64585493721610e-4|684720|6.994e2|
|two-level \(S_2\)|1/500|3.91263147067125e-4|684720|6.431e3|
|three-level Green anchor|1/100|3.23376477432167e-5|492480|5.519e2|
|three-level Green anchor|1/200|3.42594427654935e-5|492480|3.528e3|
|three-level Green anchor|1/500|3.47121897564843e-5|492480|2.104e4|

对有大量重节点的 \(S_2\) 点，\(\epsilon=1/100\) 的独立 Arb256 计算得到

\[
H_{1/100}(x)=[1.4154499328135003546703235477110860106112337060490332436078358810099238852908856
\cdot10^{-4}\pm1.90\cdot10^{-73}].
\]

它与 C++ 值的相对差为 \(2.05256\cdot10^{-13}\)。这仅给出这一点的高精度
实现校准，不为其邻域提供误差界。

minimum-near 点已有独立精确 \(G=2.250133859688266101\ldots\cdot10^{-5}\)
锚。将 C++ \(H_\epsilon\) 与该锚相比较，所观察到的相对短差依次是 3.06269%、
0.252446%、0.00378973%。这是选择未来认证中 \(\epsilon\) 尺度的有用诊断，
不是截断尾界，也不是 \(g\)、\(I\) 或 \(C_6\) 的结论。

实际命令：

```sh
python make_points.py
c++ -O3 -std=c++20 -Wall -Wextra truncated_kernel.cpp -o /tmp/cube6_truncated_kernel
/tmp/cube6_truncated_kernel ../../round003/direct_kernel/fan_integer_input.txt points.txt epsilons.txt > exploratory_results.json
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python arb_reference.py 1 0
```

前两步退出 0、wall 1.418 秒；Arb 校准退出 0、wall 18.283 秒。本有限批到此
停止，未启动积分或全空间网格。
