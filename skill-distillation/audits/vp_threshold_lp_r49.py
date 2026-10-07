# -*- coding: utf-8 -*-
"""R49: VP (Vysochanskii-Petunin) 不等式门槛的数值裁决。
离散化 LP：max P(|X-mu|>=r) over unimodal 密度（从众数 nu 起非增），
约束 E[X]=mu, Var=σ²=1。扫描 k=r/σ 与平移 d=|mu-nu|，
对照 (a) Gauss/Camp-Meidell 型 4/(9k²)（VP 目标形式）与 (b) 精确门槛。
"""
import numpy as np
from scipy.optimize import linprog

def solve_max_tail(k, d, N=3000, L=60.0):
    """unimodal: mode at 0, support [0, inf) 两侧由 even-of-mode 处理:
    一般 unimodal 密度 = 在 nu 两侧各自非增。为 LP 简化取 even (对称) 于 nu:
    这给出双侧 tail 的最坏情形族（对称化只增双边 tail）。
    网格: x in [-L, L], mode at 0; f 非增 in |x| (even unimodal)。
    d = |mu - nu|: 约束 E[X] = d (mu = d)。
    """
    # even unimodal: 变量 h_j = f 在第 j 个同心环 (|x| in [t_j, t_{j+1})), h_j 非增
    xs = np.linspace(-L, L, N)
    t = np.abs(xs)
    # 环索引: j = floor(t / dt)
    dt = t[1] - t[0]
    J = int(L / dt) + 1
    j_idx = np.minimum((t / dt).astype(int), J - 1)
    # 变量 h_0..h_{J-1} (密度高度), 非增 h_0>=h_1>=...
    # f(x) = h_{j_idx(x)}; 概率 = f * dx (近似)
    # 约束: sum f dx = 1; sum x f dx = d; sum (x-d)^2 f dx = 1 (Var=1)
    # 目标: max sum_{|x-d|>=k} f dx
    # 非增 → 用差分变量: h_j = s_j - s_{j+1}, s>=0 (标准单调化)
    M = J
    nv = M  # s_0..s_{M-1}, h_j = s_j - s_{j+1}, s_M=0
    # 矩阵 A_eq (3 x nv), A_ub (单调自动满足), 目标 c
    w = dt * np.ones(N)
    x = xs
    # h 映射矩阵 H (N x M): H[i, j] = 1 if j_idx[i]==j
    H = np.zeros((N, M))
    H[np.arange(N), j_idx] = 1.0
    tail_mask = (np.abs(x - d) >= k)
    # f = H h, h = D s where D = I - shift (h_j = s_j - s_{j+1})
    D = np.zeros((M, M))
    for j in range(M):
        D[j, j] = 1.0
        if j + 1 < M:
            D[j, j + 1] = -1.0
    F = H @ D  # N x M
    A_eq = np.vstack([
        (F.T @ w),            # total mass
        (F.T @ (w * x)),      # mean = d
        (F.T @ (w * (x - d) ** 2)),  # var = 1
    ])
    b_eq = np.array([1.0, d, 1.0])
    c = -(F.T @ (w * tail_mask))  # maximize tail
    res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=[(0, None)] * M, method='highs')
    if not res.success:
        return None
    return -res.fun

# 网格扫描: 对每个 k, 最优平移 d*(k) 下 max tail vs 4/(9k^2) 与 Cantelli 单侧
print(f"{'k':>6} {'maxP(双侧)':>12} {'4/(9k^2)':>10} {'Cantelli双':>10} {'k*门槛(4/9k^2=实际)':>8}")
for k in [1.0, 1.1, 1.1547, 1.2247, 1.3, 1.5, 1.633, 2.0, 2.5, 3.0]:
    best = None
    for d in np.arange(0, max(0.01, k - 0.5), 0.05):
        v = solve_max_tail(k, d, N=2400, L=40.0)
        if v is not None and (best is None or v > best[1]):
            best = (d, v)
    if best:
        d, v = best
        vp = 4.0 / (9 * k * k)
        print(f"{k:>6.3f} {v:>12.5f} {vp:>10.5f} {'—':>10}  d*={d:.2f}")
