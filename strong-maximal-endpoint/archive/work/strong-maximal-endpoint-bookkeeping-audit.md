# 强型稳定性到端点估计的四项账本：审计

日期：2026-11-06。令 \(Z=M_Sf\)，\(A=\int f^p\)，\(D=\int\Phi_n(f)\)，并在
\(L=|\{Z>\tau\}|>0\) 时定义
\[
 E={\tau^pL\over\|Z\|_p^p},\qquad I={A\over D},
 \qquad K=n\log c_p-\log{\|Z\|_p\over\|f\|_p}.
\]
若 \(A\le D\)，则 \(0<E,I\le1\)，且以 \(b_p=p\log c_p\) 有精确恒等式
\[
 nb_p-p\log\tau-\log{L\over D}
 =pK-\log E-\log I.
\]
将 \(M_Sf\le M_n\cdots M_1f\) 的逐坐标正规化缺口记作
\(\sum_i\log(c_p/a_i)\)，并令
\(G=\log(\|M_n\cdots M_1f\|_p/\|M_Sf\|_p)\)，则
\[
 pK-\log E-\log I
 =p\sum_i\log{c_p\over a_i}+pG-\log E-\log I.
\]
四项均非负。故仅追踪 \(L^p\) 稳定性或真实矩形几何差，都不足以推出端点 halo 的定量缺口；Markov/Chebyshev 效率 \(E\) 和 Orlicz 输入效率 \(I\) 是独立必要项。

缩放反例证实这一点。固定 \(n,p,H,\tau\)，取有界紧支撑的一维 \(L^p\) 近极值函数，并张量化。则真实强极大算子等于逐坐标复合，范数比趋近 \(c_p^n\)，各坐标 rising-sun 缺口与几何差趋于零。再乘足够小的常数使
\(\|f\|_\infty\le\min(H,\tau,1)\)。范数比和上述稳定性缺口不变，但
\(M_Sf\le\tau\)，严格 halo \(\{M_Sf>\tau\}\) 为空。此时不能把 \(E\) 代入对数恒等式（\(E=0\)）；它正是说明不存在仅由范数/几何缺口强迫正 halo 的统一下界。
