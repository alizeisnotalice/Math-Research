# 反馈熵下界的 \(O(\sqrt n)\) 精化审计

日期：2026-11-06。结论：精化成立，但似然比鞅事件必须取正确方向。

令 \(q=n\mu-d-2\log n\)，
\[
 \delta={d+2\log n+\omega\sqrt n\over n},
\]
其中 \(\omega\) 是充分大的固定常数。对每一个可预测控制，条件地向下倾斜均值超过 \(\mu-\delta\) 的增量。令 \(Q\) 为改变后的路径律，\(P\) 为原律。条件均值、方差和三阶矩的一致控制给
\[
 K_{\max}:=D(Q\Vert P)
 \le {n\delta^2\over2\mu}(1+C\sqrt\delta).
\]
在 \(Q\) 下总均值比 \(q\) 小 \(\omega\sqrt n\)，故以固定 \(\omega\) 选择
\(Q(A)\ge7/8\)，其中 \(A=\{S\le q\}\)。

若 \(\lambda_i\ge0\) 是向下倾斜的强度，
\[
 \log{dQ\over dP}=\sum_i k_i-N,
 \qquad N=\sum_i\lambda_i(Z_i-m_{Q,i}).
\]
这里 \(N\) 是 \(Q\)-鞅，\(E_QN^2\le Cn\delta^2\)。**符号关键点**：为上界 \(\log(dQ/dP)\)，必须在事件
\(N\ge-C_1\sqrt n\delta\) 上工作（或用 \(|N|\) 的双侧事件）；控制 \(N\le+C_1\sqrt n\delta\) 的方向不足。修正后，交集有 \(Q\)-概率至少 \(3/4\)，并给出
\[
 P(A)\ge {3\over4}\exp\{-K_{\max}-C_1\sqrt n\delta\}.
\]

对 \(0\le d\le2\sigma\sqrt{n\log n}\)，展开得到
\[
 K_{\max}+C_1\sqrt n\delta
 \le {d^2\over2\mu n}+C\left({d\over\sqrt n}+1\right).
\]
这是因为 \(n\delta^2\sqrt\delta=o(1)\)，而 \(2\log n\) 的交叉项也被右端吸收。故当
\(d=\sigma\sqrt{n\log n}-K\sqrt n\) 时，低尾下界的指数为
\[
 -\tfrac12\log n+\left({K\over\sigma}-C\sigma\right)\sqrt{\log n}-O(1).
\]
选择固定充分大的 \(K\) 后，它严格大于近阈值轮廓上界要求的 \(O(\sqrt{\log n/n})\)，从而给出 \(O(\sqrt n)\) 余项级别的反馈模型上界。该结论仍依赖先前的 \(T^n\) 全局重排归约。
