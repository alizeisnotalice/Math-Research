# 停止小均值策略的自适应 Gaussian 低尾：审计

日期：2026-11-06。此记录审计停止/耦合部分；相对 Gaussian 中偏差结论仍须以所选 Fan--Grama--Liu 型定理的完整假设为准。

令 \(D_k=\sum_{i\le k}(\mu-m(r_i))\)，并令 \(\tau\) 为首次
\(D_k>A\sqrt{n\log n}\) 的时刻。停止策略 \(\pi_K\) 在 \(\tau\) 后恒取
\(r_*=B^{-1}\)。于是
\[
 D_n^{(K)}\le A\sqrt{n\log n}+\mu-1,
 \qquad |\langle M_K\rangle-\mu n|
 \le C\sqrt{nD_n^{(K)}}=O(n^{3/4}(\log n)^{1/4})
\]
逐路径成立。所有条件 Bernstein 常数一致：\(Z_r=Y+cJ\) 中稀有大移位的矩由
\(c^ke^{-c}\) 的有限上确界控制。故若所用鞅中偏差定理接受这种逐路径二次变差误差和条件 Bernstein 条件，则在 \(x\le2\sqrt{\log n}\) 有一致相对 Gaussian 低尾。

耦合方向必须如下处理。让原策略 \(\pi\) 和 \(\pi_K\) 在 \(\tau\) 前共用同一增量；无停止时二者完全相同，且
\(M_K\le-d\) 蕴含原策略失败 \(S_n\le n\mu-d\)。停止后不应直接把 Gaussian 事件扣除。相反，原策略若已停止却成功，则
\[
 M_n=S_n-n\mu+D_n>D_\tau-d
 \ge(A-2\sigma)\sqrt{n\log n}.
\]
原策略的逐步条件正 mgf 一致界给该事件概率至多 \(n^{-3}\)，只要固定 \(A\) 足够大。因此
\[
 P_\pi(S_n\le n\mu-d)
 \ge P_{\pi_K}(M_K\le-d)-n^{-3}.
\]
这一步保留了停止事件本身的失败质量，避免错误地损失一个可能很大的停止概率。若中偏差定理的条件核验完成，该式可把停止策略的相对 Mills 下界统一转给所有可预测策略。
