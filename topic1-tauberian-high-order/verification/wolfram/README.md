# Mathematica 精确核验

实际运行环境：Wolfram 15.0.0，macOS ARM。三个脚本各从独立内核运行，不依赖 Notebook 状态。

- `gram.wl`：带完整实参数假设的 Resolve/Reduce。
- `energy.wl`：0<ε<1 上的精确分段积分；缓存正部符号，检测区间内部根。
- `endpoint.wl`：n=2,…,8 的形式系数提取与代回比较。
- `*-result.txt`：实际成功运行的关键输出。

运行方式：`wolframscript -file verification/wolfram/gram.wl`，其余文件同理。
符号结果支持各自的代数步骤；一般几何证明、极限存在性及余项控制仍须按论文的证明链核对。

## 第4轮：真实二阶证明的标量核验

`sharp-check.wl` 为独立审查者的精确检查，主代理在本机 Wolfram 15.0.0 ARM 重跑；16项均 True，Failures=0。输出见 `sharp-check-result.txt`。脚本无运行时限，显式列出每个代数断言的假设。它支持正文解析证明中的常数和符号核对，不替代输入全域、可测选择和极限传递的证明。大变量联合CAD的内存未决运行未计作定理证据，见 `rounds/004.md`。

## 第5轮：有限阈值障碍与任意铺砌等孔

`threshold-obstruction.wl` 核对有限尾分布、预算及阈值成本，`tiling-holes.wl` 核对角接触导数、面邻孔的全维余量及部分填孔优化。两者在本机 Wolfram15 实际运行，输出分别保存为 `*-result.txt`。几何全半径、全部阈值及有限输入的证明见 `proofs/005-*.md`；符号等式不自行认证几何可实现性。
