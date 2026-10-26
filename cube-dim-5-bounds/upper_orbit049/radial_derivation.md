# 径向 Gamma 精确消去：运行前独立推导

固定029的16个独立正 stable \(1/2\) 变量，Laplace 归一化为
\(\mathbb E e^{-uS_i}=e^{-\sqrt u}\)。若 \(Y_i\sim\Gamma(1/2,1)\)，则
\(S_i=(4Y_i)^{-1}\) 正是此归一化：对 Gamma 密度积分得到
\(\pi^{-1/2}\int_0^\infty y^{-1/2}e^{-y-u/(4y)}dy=e^{-\sqrt u}\)。

令 \(R=\sum_iY_i\)、\(Z_i=Y_i/R\)。Gamma–Dirichlet 换元的 Jacobian 为 \(R^{15}\)，联合密度分解为
\[
 \frac{R^{8-1}e^{-R}}{\Gamma(8)}\,
 \frac{\Gamma(8)}{\Gamma(1/2)^{16}}\prod_iZ_i^{-1/2}.
\]
因此 \(R\sim\Gamma(8,1)\)，\(Z\sim\operatorname{Dir}(1/2,\ldots,1/2)\)，并且两者独立。
定义
\[
 H(Z)=\sum_iZ_i^{-1},\qquad
 P_i=\frac{Z_i^{-1}}{H},\qquad M=\sum_iP_iv_iv_i^T.
\]
原 stable 总尺度是
\[
 T=\sum_iS_i=\frac{H(Z)}{4R},\qquad T^{-1/2}=2R^{1/2}H(Z)^{-1/2}.
\]
这里 \(T\) 依赖 \(R,Z\)，不能把 \(T\) 与 \(P\) 设为独立。

另一方面 \(T\) 的 Laplace transform 是 \(e^{-16\sqrt u}\)，故
\[
 Z_{\rm tilt}=\mathbb ET^{-1/2}
 =\frac1{\sqrt\pi}\int_0^\infty u^{-1/2}e^{-16\sqrt u}du
 =\frac1{8\sqrt\pi}.
\]
对于非负可测 \(F(M)\)（或绝对可积的有符号 \(F\)），Tonelli 与已证 \(R,Z\) 独立性给
\[
 \begin{aligned}
 \mathbb E_{R_{1/2}}F
 &=\frac{\mathbb E[T^{-1/2}F(M)]}{Z_{\rm tilt}}\\
 &=16\sqrt\pi\frac{\Gamma(17/2)}{\Gamma(8)}
   \mathbb E_{\rm Dir}\left[H^{-1/2}F(M)\right].
 \end{aligned}
\]
因此精确权重为
\[
 w(Z)=C H^{-1/2},\quad C=16\sqrt\pi\frac{\Gamma(17/2)}{\Gamma(8)}=\frac{6435\pi}{256},
 \qquad \mathbb E_{\rm Dir}w=1.
\]
由 \((\sum_i Z_i)(\sum_iZ_i^{-1})\ge16^2\) 得 \(H\ge256\)，所以 \(0<w\le C/16\)。

本轮直接抽样 Dirichlet 后使用此权重。样本均值不再除以样本权重均值，否则会换成自归一化估计。
径向变量没有留下抽样误差；Dirichlet 形状积分、有限求积和浮点矩阵计算仍全部是非认证诊断。

群能量沿038：\(r_\Gamma=\operatorname{avg}_U\rho_U\)，
\[
 \rho_U=\mathbb E_{t\sim\operatorname{Beta}(5/4,5/4)}
 [\det M\det(tM^{-1}+(1-t)UM^{-1}U^T)]^{-1/2}.
\]
全局 \(U\) 与 \(-U\) 的共轭相同，所以1920个代表与3840个 signed permutations 的平均完全一致。
沿用同一个 \(c_1=1/(32\pi^3)\)、\(q=5/4\) 和 face 系数 \(\pi^2/60\)；不改动029的body半群与全尺度资格。
