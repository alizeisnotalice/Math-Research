# 第004轮 spatial：归档口径与独立几何/标量复核

本目录只归档本代理的空间连续化接口；它不是完整的六维上界证书，也不宣布
任何新的 (C_6) 界。

## 已严格闭合的接口

令 (V=\{(1,s_2,\ldots,s_6):s_j\in\{\pm1\}\})，(A c=\sum_{v\in V}c_vv)，
并令 (H_\epsilon=\int_\epsilon^\infty p_t\,dt)。报告 005--006 严格给出：

1. \(\operatorname{conv}(\pm V)=[-1,1]^6\)，从而参数平移的最小
   \(\ell^1\) 成本为 \(\|h\|_\infty\)，并有
   \[
   e^{-\|h\|_\infty/\epsilon}H_\epsilon(x)
   \le H_\epsilon(x+h)\le
   e^{\|h\|_\infty/\epsilon}H_\epsilon(x).
   \]
2. \(AA^T=32I\) 的 coarea 表示给出
   \(H_\epsilon(x)e^{\|x\|_2^2/(32\epsilon^2)}\) log-convex，故可由 face
   box 顶点严格上界积分上界 \(H_\epsilon^{6/5}\)。
3. 沿任一坐标 \(\partial_{ii}H_\epsilon\le(129/128)\epsilon^{-2}H_\epsilon\)。
   双端点 hyperbolic-sine 比较在五个自由坐标迭代，给局部 cube/face 的严格
   下界；这是用于本轮局部 box 的几何接口。

本轮的复参数相对管也只在逐一核验
\(\left|(A^Th/32)_v\right|<\epsilon\) 时使用；没有把单轴宽度错误推广成
多变量 max-norm 管。

## 实际运行的独立标量复核

`local_box_independent_audit.py` 由本机
`/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python` 实际运行，
退出码 0、用时约 0.1 秒。它读取主项目（只读）的 Arb384 输入
`arb384_results.json`，输入 SHA-256 为
`175a8de838642c003ebbef99ae5d62e0ff9c92e6803aa231c904d96061573708`。

审计的 box 为
\[
 E=\{x_1=1,\ x_2,\ldots,x_5\in[.749,.751],\ x_6\in[0,.001]\},\qquad
 \epsilon=.01.
\]
四个 .749/.751 坐标的五类代表有 \(\binom4k\) 个排列；第五坐标的 \(\pm\)
反射使完整 box 为32个顶点，故顶点平均正是
\(\sum_{k=0}^4\binom4kH_k^{6/5}/16\)。先在完整第五轴 box 使用顶点不等式、
再按反射取正半部，得到因子 \(1/2\)，而不是重复或漏掉一个32因子。

审计输出：
\[
 \cosh(\sqrt{129/128}/10)^5=1.025472114\ldots,
 \qquad G>3\cdot10^{-5}\quad\hbox{on }E,
\]
以及
\[
 \int_EH_{.01}^{6/5}<6.5551507514\cdot10^{-20}<7\cdot10^{-20}.
\]
后一个量**仅**是局部截时势的积分上界，绝不可替换为 \(I=\int G^{6/5}\)
或 \(G\) 的积分上界。

作为纯代数的条件检查，若未来独立证明
\(g\ge1/60000\) 且 \(\|H_{.001}(1,\cdot)\|_{6/5}\le1/16000\)，再接第002轮的
whole-face tail，成本为 \(61.6895427961\ldots<64\)。这两个全局前提都尚未证明，
所以该条件句不是上界结果。

本代理输出 SHA-256：

- `local_box_independent_audit.json`：`a6ba28768302865fa655279a52d849c014fcc2f9224d33d0b20d1986c4222205`；
- 读取的主项目局部证书 JSON：`54b032d21d4e4c907840d556bb3e9370081ae6655b7843265f0a41fabf667588`。

## 交叉实现的只读记录

主项目实际运行的单角点第二个显式实导数 DD 实现报告
`PASS_second_formula_corner_crosscheck`（10.732085 秒，输入 SHA-256 同上）；
其输出 JSON SHA-256 为
`4ca253de6f9e394f1c5f93db4e087d9073a8626a1796ffed4472a6a64a4f32f5`。
整数 fan 至 Fraction fan 的射线、rate、det 缩放桥报告
`PASS_exact_ray_rate_determinant_scaling_bridge`，其输出 JSON SHA-256 为
`52d04e224b82dbf492f7330e1cfa1e043e983bb3fa367c540271367065d36a6b`。
它们是主项目的运行证据；本代理本轮未重跑其 10 秒 DD 计算，也没有把它们
升级为全局连续认证。

## 本轮未闭合的必要全局对象

仍需一个覆盖完整 cube 的 \(g\) 下界，以及正 face 的
\(\|H_{.001}\|_{6/5}\) 全域严格上界，并将其与短时 tail 合并。局部 box、
少数顶点和条件成本均不能替代这两项。均匀网格按照 005 的 cell 计数仍不经济；
下一步应使用完整的自适应覆盖/误差账，而不是删除未采样区域。
