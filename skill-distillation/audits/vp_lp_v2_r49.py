# -*- coding: utf-8 -*-
"""R49 v2: VP 门槛裁决——一般 unimodal（mode 两侧独立非增），mode=0。
变量: 左高度 l_j (x in [-(j+1)dt, -j dt]), 右高度 r_j。各自非增。
约束: mass=1, E[X]=d (即 mu=d), Var=1。
目标: max P(|X-mu|>=k) [双侧] / P(X-mu>=k) [单侧]。
对照: VP 双侧 4/(9k^2)；Camp-Meidell 单侧 4/(9(k^2+1))。
"""
import numpy as np
from scipy.optimize import linprog

def solve(k, d, N=1600, L=25.0, oneside=False):
    xs = np.linspace(-L, L, N); dt = xs[1] - xs[0]
    J = int(round(L / dt)) + 1
    # 环索引：右支 j = floor(x/dt) for x>=0；左支 j = floor(-x/dt) for x<0
    right = xs >= 0
    jr = np.where(right, np.minimum((xs / dt).astype(int), J - 1), -1)
    jl = np.where(~right, np.minimum((-xs / dt).astype(int), J - 1), -1)
    M = 2 * J
    # 变量: [r_0..r_{J-1}, l_0..l_{J-1}]
    Wr = np.zeros((N, M)); Wl = np.zeros((N, M))
    mr = right & (jr >= 0); ml = (~right) & (jl >= 0)
    Wr[np.arange(N)[mr], jr[mr]] = dt
    Wl[np.arange(N)[ml], J + jl[ml]] = dt
    W = Wr + Wl   # 概率质量算子: W @ v, v=[r;l]
    tail = (np.abs(xs - d) >= k) if not oneside else ((xs - d) >= k)
    # 单调: r_j >= r_{j+1}, l_j >= l_{j+1}
    A_ub = np.zeros((2 * (J - 1), M))
    for j in range(J - 1):
        A_ub[j, j] = -1.0; A_ub[j, j + 1] = 1.0            # 右支非增
        A_ub[J - 1 + j, J + j] = -1.0; A_ub[J - 1 + j, J + j + 1] = 1.0  # 左支非增
    b_ub = np.zeros(2 * (J - 1))
    # 矩约束: 每环质量 / x 贡献 / (x-d)^2 贡献, 用 bincount 精确求和
    m_r = Wr[:, :J].sum(axis=0); m_l = Wl[:, J:].sum(axis=0)
    # 环内 x 的和: 用 bincount
    sx_r = np.bincount(jr[mr], weights=xs[mr] * dt, minlength=J)
    sx_l = np.bincount(jl[ml], weights=xs[ml] * dt, minlength=J)
    s2_r = np.bincount(jr[mr], weights=(xs[mr] - d) ** 2 * dt, minlength=J)
    s2_l = np.bincount(jl[ml], weights=(xs[ml] - d) ** 2 * dt, minlength=J)
    A_eq = np.vstack([
        np.concatenate([m_r, m_l]),
        np.concatenate([sx_r, sx_l]),
        np.concatenate([s2_r, s2_l]),
    ])
    b_eq = np.array([1.0, d, 1.0])
    # 目标: tail 集合的环内权重
    t_r = np.bincount(jr[mr & tail], weights=dt * np.ones((mr & tail).sum()), minlength=J)
    t_l = np.bincount(jl[ml & tail], weights=dt * np.ones((ml & tail).sum()), minlength=J)
    c = -np.concatenate([t_r, t_l])
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                  bounds=[(0, None)] * M, method='highs')
    if not res.success:
        return None
    return -res.fun

if __name__ == '__main__':
    # 可行性快测: d=0.5, k=1.5
    v = solve(1.5, 0.5)
    print('feasibility test d=0.5,k=1.5, oneside tail max =', v)
    v = solve(1.5, 0.5, oneside=False)
    print('feasibility test d=0.5,k=1.5, 双侧 tail max =', v)
