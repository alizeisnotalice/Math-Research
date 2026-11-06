# 反馈控制熵下界：独立审计

日期：2026-11-06。审计结论：所述熵变测度路线在明确“策略对过去可预测、轮廓为单调且已完成对数变量平移归一”的条件下成立；未发现自适应性导致的方向或常数漏洞。

设 \(r_i\) 是过去增量的可预测函数，条件增量为 \(Z_{r_i}\)。取
\[
 d=(\sigma-\varepsilon)\sqrt{n\log n},\quad
 \omega=(\log n)^{1/4},\quad
 \delta={d+\omega\sqrt n\over n}.
\]
当 \(m(r)>\mu-\delta\) 时，将该步以 \(e^{-\lambda Z_r}\) 倾斜至条件均值 \(\mu-\delta\)；其余步不变。近极值性强制
\(r\) 距 \(B^{-1}\) 为 \(O(\sqrt\delta)\)，故 \(v(r)=\mu+O(\sqrt\delta)\)。一致累积量展开给每步
\[
 D(Q_i\Vert P_i)\le {\delta^2\over2\mu}(1+O(\sqrt\delta)),
\]
并由链式规则
\[
 D(Q\Vert P)\le {n\delta^2\over2\mu}(1+o(1)).
\]
在 \(Q\) 下条件均值至多 \(\mu-\delta\)，条件方差一致有界；鞅 Chebyshev 因而给
\(Q(S>n\mu-d)=O(\omega^{-2})\)。二元熵不等式（保留 \(+\log2\) 的常数项）推出
\[
 P(S\le n\mu-d)\ge
 \exp\left\{-{n\delta^2\over2\mu}(1+o(1))-O(1)\right\}.
\]
由于 \((\sigma-\varepsilon)^2/\sigma^2<1\)，该量严格大于
\(\sqrt{\log n/n}\)，足以压倒近阈值轮廓的 \(O(\sqrt{\log n/n})\) 期望超额。

自适应不会破坏远尾估计：条件密度满足
\(f_r(z)\le2e^{-z}\)，所以任何可预测策略的联合密度至多
\(2^ne^{-\sum z_i}\)，从而总和密度有同一 Gamma 型上界。条件 mgf 的一致二次界也可逐步条件化得到 Chernoff 上界。故原固定比例证明的高轮廓远尾步骤可移植。

仍需在正式定理中保留两项定义性工作：先由坐标重排比较把一般 \(C_n\) 归约到 \(T^n\) 的单调一维模型，再用平移把一个 endpoint level 归一到控制过程的初态。固定比例 \(Q_n\) 的方差结论不能直接替代这条熵论证；本论证正是为反馈漂移而设计。
