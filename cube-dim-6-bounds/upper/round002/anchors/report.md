# 六维 body frame 的有限截时锚点结果

`truncated_anchor_batch.py` 使用理论接口

\[
H_\epsilon(x)=(2\pi)^{-6}\sum_R\alpha_R
\operatorname{DD}_5[\arctan(\cdot/\epsilon);b_{R,1},\ldots,b_{R,6}]
\]

对 12 个非 \(S_2\) 三层锚点给出 384-bit Arb 区间。完整区间写在
`truncated_anchor_results.json`；这里列出各点 \(\epsilon=1/1000\) 的中值，
仅方便定位后续连续认证的候选区域：

|\(k\)|\(a=1/4\)|\(a=1/2\)|\(a=3/4\)|
|---:|---:|---:|---:|
|1|2.385705947595738e-4|1.655427765360434e-4|1.346719481631153e-4|
|2|2.001558526053102e-4|1.173174817483790e-4|7.086723498160538e-5|
|3|1.776632937912525e-4|9.064396627437515e-5|4.459524268505753e-5|
|4|1.627857026339174e-4|7.742793922325810e-5|3.480797238324517e-5|

最小显示值是 \((1,(3/4)^4,0)\) 的约 \(3.4808\cdot10^{-5}\)。它只是一个
有限 \(H_\epsilon\) 样本；虽由正性给出该点 \(G\) 的下界，但不能作为全 cube 的统一下界、minimum、或
成本中的 \(g_0\)。所有区间半径及另一个 \(\epsilon\) 结果见 JSON，不以表中
小数作任何比较。

原点校准在两个 \(\epsilon\) 下均通过。首锚点另作群压缩复核：完整群有
46080 元，轨道有 120 点、稳定子 384；21 个 base simplex 的未合并 orbit 和
含 2520 项，合并为 2512 个精确有理节点 tuple，两个 \(\epsilon\) 的差区间均
包含零。该复核验证了本批的群作用和节点聚合，未验证连续 Green 论证。

运行：pilot 约 0.90 秒，全 12 锚点约 15.58 秒，聚合复核约 1.36 秒；三个
进程均 exit 0，无保留句柄。

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python truncated_anchor_batch.py --phase all --precision 384
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python verify_anchor_aggregation.py
```

下一步若推进连续认证，应从这些非奇异层附近的有理 box 着手，并使用已经证明
的 tube/Chebyshev 误差账；不得把本批插值为全 cube 结论。
