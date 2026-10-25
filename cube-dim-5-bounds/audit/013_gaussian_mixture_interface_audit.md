# 013 — Gaussian 混合全域 Green 下界的结构审计

审阅对象：主 clone `round004/gaussian_mixture_lower.md`，本次 SHA-256
`6443d519bda7e1db1b48d09f562a7ea6af47e6e3685e53ec352a2a5421da5e48`。
本报告只判定结构引理；没有估计所出现的正期望，因而不把它称作数值 \(g\)
证书或 \(C_5\) 改进。

## 结论

结论成立，须在正式版本补写一个 determinant-invariance 的加权 Jensen
步骤。归一化常数、32 个符号平均、Jensen 方向、Tonelli 与 \(\infty\)
情形均正确。

令 \(V=\{(1,\pm1,\pm1,\pm1,\pm1)\}\)，\(|V|=16\)，独立
\(N_v,M_v\sim N(0,1)\)。\(N_v/|M_v|\) 与 \(N_v/M_v\) 同分布，均为
标准 Cauchy；故 \(\sum_vtvN_v/|M_v|\) 的 characteristic exponent 是
\(t\sum_v|v\cdot\xi|\)，即恰为既定模型的 \(p_t\)，没有遗漏 Cauchy
scale 或第 17 个方向。

条件于 \(M\)，其 covariance 为
\[
t^2A,\qquad A=\sum_{v\in V}{vv^T\over M_v^2}.
\]
所有 \(M_v\ne0\) a.s.，且方向族张成 \(\mathbb R^5\)，所以 \(A\) 正定
a.s.。写 \(a=x^TA^{-1}x>0\)，条件密度是
\[
(2\pi)^{-5/2}t^{-5}(\det A)^{-1/2}e^{-a/(2t^2)}.
\]
直接符号检查（本机 SymPy）给
\[
\int_0^\infty t^{-5}e^{-a/(2t^2)}dt={2\over a^2};
\]
故 Tonelli 给
\[
G(x)=2(2\pi)^{-5/2}
 E\left[{(det A)^{-1/2}\over(x^TA^{-1}x)^2}\right],\quad x\ne0.
\]
非负性允许等式值为 \(+\infty\)，没有非法交换。\(x=0\) 时同一条件表达式
的时间积分发散，故 \(G(0)=\infty\)。

## 对称性与 Jensen

每个对角符号矩阵 \(S\) 使 \(Sv\) 等于方向族某元素或其负元；外积不变且
16 个权重 iid，因此 \(A\stackrel d=SAS\)。不过正式证明必须保留
\((\det A)^{-1/2}\) 权重：它同样在 \(A\mapsto SAS\) 下不变，所以
\[
 E[D(A)f(q_A(x))]
 =E\left[D(A)2^{-5}\sum_S f(q_A(Sx))\right],
\]
其中 \(D=(\det A)^{-1/2}\), \(q_A(y)=y^TA^{-1}y\), \(f(z)=z^{-2}\)。
现在才可对每个固定 \(A\) 用凸的 \(f\) 作 Jensen。符号平均消去交叉项：
\[
2^{-5}\sum_Sq_A(Sx)=\sum_i(A^{-1})_{ii}x_i^2.
\]
所以对 \(x\in[-1,1]^5\setminus\{0\}\)，其不大于 \(\operatorname{tr}A^{-1}\)，
且因 \(z^{-2}\) 单调递减，得到文中的方向
\[
G(x)\ge2(2\pi)^{-5/2}E\left[{(det A)^{-1/2}\over(\operatorname{tr}A^{-1})^2}\right].
\]

期望的 integrand a.s. 严格正且有限。它的期望先允许为 \(+\infty\)；从已知
\(G\in L^{5/4,\infty}\) 选一个 cube 内 a.e. 有限的非零点，再用上述点态
下界，才推出该期望有限。这个顺序正确，也不需假定 \(G\) 在奇异子空间
连续或有限。

## 限制

该结构引理只将全域 \(\inf_xG(x)\) 问题变成一个正期望的严格下界问题。它
尚未提供该期望的有理/区间下界，也未给同一 \(G\) 的面 \(I\) 上界；因此
当前不产生 \(g\)、成本或 \(C_5\) 的认证数值。MC 仅可作量级诊断。
