# 第十一轮：一般无原子 Radon 测度的临界构造

## 已证主结论

对任意非零、局部有限、无原子的 Radon 测度 \(\mu\) 于 \(\mathbb R^n\)，
以及任意正指数向量 \(\boldsymbol a\)，共同中心、公共交集归一化的内禀
常数满足

\[
C_{c,\mu}^{\boldsymbol a}(\lambda)=\infty
\quad\text{for }0<\lambda\le P_{\boldsymbol a}(1),
\qquad
P_{\boldsymbol a}(1)=\max_{\sum d_i=1}\prod_i(1-d_i)^{a_i}.
\]

因此一般非负非零局部可积权 \(w\,dx\) 已包含在内；等指数端点为
\(((m-1)/m)^m\)。`general_critical.tex` 给出可直接 `input` 的完整证明。

## 构造和依赖次序

先把有界的共同集合 \(F\) 全部置于一个正比例色 \(j\)。第 \(k\) 级的
所有见证方体先包含旧区域 \(D_{k-1}\)。对每个中心选择一个充分大、
边界 \(\mu\)-零的半径；其局部中心扰动满足 \(\mu(Q\triangle Q')\to0\)。
有限覆盖后，只在新区域的有限 Boolean 原子上按最优比例 \(\boldsymbol d\)
精确分色。

由此，偏差 \(\mu(A_i\cap D_k)-d_i\mu(D_k)\) 在所有层级恒为
\((1-d_j)\mu(F)\)（第 \(j\) 色）或 \(-d_i\mu(F)\)（其余色）。每个参考
方体故有精确余量 \(d_i\mu(F)\)。若
\(\mu(Q\triangle Q')<\varepsilon<\mu(F)\min_{d_i>0}d_i/4\)，实际方体
仍有至少 \(d_i\mu(F)-2\varepsilon>0\) 的余量。后续只在见证方体外
着色，故见证永久有效。

无原子分割引理在正文中自证：用可数开基的三进制编码给出 Borel 单射，
推前测度的连续分布函数可精确切出每个 Boolean 原子的指定比例。

## 阈值和反例

结合前面已证的超临界传递与统一弱型上界，得到
\(0<\lambda<1\) 时 \(C_{c,\mu}^{\boldsymbol a}(\lambda)<\infty\) 当且仅当
\(\lambda>P_{\boldsymbol a}(1)\)。无原子性不可删除：\(\mu=\delta_0\) 时
每个允许公共交集和光晕均有质量 1，所以 \(C=1\)（\(0<\lambda<1\)）。

## 本机计算

- `python3 general_critical_certificate.py`：以 `Fraction` 精确核验三次
  新原子按比例分色后的偏差不变量、参考方体余量 \(d_if\)，以及
  \(d_if-2\varepsilon>0\)。
- `wolframscript -file general_critical_check.wl`：独立符号化简相同余量。

这些程序仅是有限代数证据；一般测度、半径选择、有限覆盖及可测分割均由
`general_critical.tex` 的证明负责。


## 主代理验收

主代理实际重跑Python与Wolfram，均通过；独立审查另行确认一般无原子Radon版本。主文保留严格的局部半径余量和有界中心邻域，避免把wdx的统一平移连续性不加条件用于奇异测度。无原子分割用开基三进制Borel单射和连续CDF自证。有限质量、乘积权旧结果由本轮涵盖，当前正文不再重复输入；旧构造脚本及片段为相应轮次复现所需而保留。

整合稿9页，XeLaTeX双遍通过，无未定义引用或Overfull；新增临界章节和末页已渲染核看。旧的特殊权构造从正文移出，其独立片段作为历史轮次复现附件保留。
