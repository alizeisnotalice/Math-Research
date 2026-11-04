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

## 第6轮：空间缓冲

`spatial-buffer.wl` 的14项精确断言在 Wolfram15 ARM 实际全部通过，主代理重跑输出一致。它检查有限路径、端点插值、局部预算及指数均衡。实际集合几何及全输入证明为 `paper/sections/global.tex`，另有两份独立审查与最终转写复核。

## 第7轮

conditional-third.wl、finite-boxes.wl、tail-moments.wl、tail-gap.wl、third-scale.wl均从归档路径重跑，输出与原成功运行逐字一致。satellite-budget.wl在15子内核上实际15项check=True。全捕获形式系数须另有几何证明。conditional-third输出的False指整数违反集为空。tail-gap列出较强的辅助参数域，完整替代被中断CAD的同域数学断言，没有声称原CAD成功。

## 第8轮：异常回收与多线性第三项

`global-third.wl`由必要独立审查者实际执行；所有违反集为空，固定常数及分母恒等式符合解析证明。`multilinear-third.wl`在主代理本机执行，最终Failures=0，核对一致三阶Taylor余量及二维系数66−2/m。第一次简化留下未决布尔式后改用等价的一元Reduce，未决运行不计通过。各自成功输出保存为`*-result.txt`。`../constant_hierarchy.py`以15进程和精确有理数检查n=2,…,16，输出`../constant_hierarchy.json`，不是全部维数或全部几何的有限枚举证明。

当前`paper/sections/global.tex`已改用第8轮证明；第6轮空间缓冲原稿在Git提交bec584f及其历史审查中恢复。

`next-scale-ledger.wl`实际5项True，给出两分支抽象标量信息的精确minimax三次损失；不声称这些标量自动由真实输入实现。
