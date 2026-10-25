# 11J 高阶有限批本地执行登记

2026-10-25。只新建 cc11_cheb_*，不覆盖首批锁证据。固定 τ=1/3、ε=1/50、16×16 方格中 136 个 i≥j patches；每坐标 degree 20 的 Chebyshev roots，ellipse 参数 3。

先执行四个端到端代表 patch (15,7)、(15,15)、(8,3)、(0,0)，分别覆盖 edge 数值候选、body 方向奇线/diagonal、一般内部、zero/repeated 投影。成功后执行固定全部 136 patches，不新增网格或参数。

初始算术选 384-bit Arb；Chebyshev root 以分母 2^128 的精确 dyadic 数替代（成对反对称，中央节点准确零），原节点位移由全纯 Cauchy 导数界支付。数值节点 enclosure 通过准确 rational-point Hermite divided difference 求得；不以 Arb 区间重叠判定准确重节点。

用 Arb DCT 与准确有理 Chebyshev→Bernstein 矩阵得到二元 Bernstein enclosure。每 patch 先以全部 Bernstein 系数下端认证；未通过时用 de Casteljau 严格二分，交替分 y/z，最大 binary depth=8、最多 511 个节点/256 个最终叶子。每 leaf 必须满足 polynomial 下界减统一插值余项严格大于 13/2500；超限则记录失败，不机械扩大分支深度。节点位移与所有舍入误差传播在 Arb 球内，不再另扣两次。

只允许一次由明确失败原因触发的表示或精度修订；必须先记录原因。初始 384-bit 已为 Bernstein 转换和近重节点的相消预留精度，不以首批 1e-12 节点误差作为充分算术目标。

可再生节点缓存置仓库外 `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/cc11_cheb_cache/`；最终保留必要 checker、严格误差推导及每 patch/leaf 的证书。实际 runtime、版本和 hash 随交付锁定。全域成功后送主及原 independent_audit，自己自检不代替独审。
