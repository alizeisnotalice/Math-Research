# 多面体切向泛函：反射对偶、LP 支持函数表示与本机核验

设 \(K\subset\mathbb R^n\) 是含原点于内部的满维凸多面体，外法向量为 \(\nu_i\)，相应面面积为 \(A_i>0\)，原点支持数为 \(h_i=\langle y,\nu_i\rangle>0\)（\(y\) 在第 \(i\) 个面上）。记 \(V=|K|\)。对 \(0<\beta<n\)，固定模型 \(E=K\) 的一阶边界缺陷为

\[
D_K(v,a)=\sum_iA_i\bigl(-\langle v,\nu_i\rangle-ah_i\bigr)_+,
\qquad
\Phi_{\beta,K}(v,a)=(n-\beta)a+\frac{D_K(v,a)}V,
\]

并令 \(N_{\beta,K}(v)=\inf_{a\in\mathbb R}\Phi_{\beta,K}(v,a)\)。以下证明只涉及这个固定模型；不涉及全体可测集的 Tauberian 上确界。

## 反射对偶

凸多面体的平衡恒等式为

\[
\sum_iA_i\nu_i=0,
\qquad
\sum_iA_ih_i=nV.
\]

由 \(x_+=(-x)_++x\)，有

\[
\begin{aligned}
D_K(-v,-a)
 &=\sum_iA_i(\langle v,\nu_i\rangle+ah_i)_+\\
 &=D_K(v,a)+\sum_iA_i\langle v,\nu_i\rangle+a\sum_iA_ih_i\\
 &=D_K(v,a)+naV.
\end{aligned}
\]

因此

\[
\Phi_{n-\beta,K}(-v,-a)=\Phi_{\beta,K}(v,a),
\qquad
N_{n-\beta,K}(-v)=N_{\beta,K}(v).
\]

这给出参数 \(\beta\leftrightarrow n-\beta\) 下的精确反射对偶；它不要求 \(K\) 中心对称。

## 有限维 LP 对偶

考虑原始线性规划

\[
\begin{aligned}
 \text{最小化}\quad &(n-\beta)a+\frac1V\sum_iA_iu_i,\\
 \text{满足}\quad &u_i+h_ia\ge-\langle v,\nu_i\rangle,\qquad u_i\ge0,\quad a\in\mathbb R.
\end{aligned}
\]

给定 \(a\) 后，最优 \(u_i\) 正是 \(( -\langle v,\nu_i\rangle-ah_i)_+\)，所以原始最优值为 \(N_{\beta,K}(v)\)。其 LP 对偶为

\[
\begin{aligned}
 \text{最大化}\quad&\left\langle v,-\frac1V\sum_iA_iq_i\nu_i\right\rangle,\\
 \text{满足}\quad&0\le q_i\le1,\qquad \sum_iA_ih_iq_i=(n-\beta)V.
\end{aligned}
\]

定义

\[
P_{\beta,K}=\left\{-\frac1V\sum_iA_iq_i\nu_i:0\le q_i\le1,\ \sum_iA_ih_iq_i=(n-\beta)V\right\}.
\]

原始问题可行且有有限最优值：\(a\to+\infty\) 与 \(a\to-\infty\) 时目标分别由 \((n-\beta)a\) 和 \(-\beta a\) 强制发散。于是有限维 LP 强对偶给出

\[
N_{\beta,K}(v)=h_{P_{\beta,K}}(v).
\]

取 \(q_i=(n-\beta)/n\) 可见约束超平面穿过立方体的相对内部，并且其像包含原点。由 \(N_{\beta,K}(v)>0\)（\(v\ne0\)）可知 \(P_{\beta,K}\) 满维且原点在内部。故切向体满足

\[
S_{\beta,K}=\{v:N_{\beta,K}(v)<1\}=\operatorname{int}(P_{\beta,K}^{\circ}).
\]

对偶还能在多面体层面看见反射：\(q\mapsto1-q\) 把 \(P_{\beta,K}\) 双射到 \(-P_{n-\beta,K}\)。这与上面的函数恒等式完全一致。

## 实际本机核验

运行文件：

- `verification/polytope_tangent.wl`：Wolfram Language 15.0.0，完全有理/代数精确算术；枚举切片立方体顶点、计算二维凸包和极体面积，并对五个精确方向检查 LP 原始值与支持函数。
- `verification/polytope_tangent.m`：MATLAB R2026a Update 4，`linprog` 双单纯形法，在 64 个均匀圆周方向上独立检查原始 LP、对偶支持函数和反射对偶，并用 `convhull`/`polyarea` 计算二维边界与面积。

| \(K\) | \(\beta\) | Wolfram 精确 \(|S_{\beta,K}|\) | \(|S_{\beta,K}|/|K|\) | MATLAB 最大原始--对偶误差 |
|---|---:|---:|---:|---:|
| \([ -\tfrac12,\tfrac12]^2\) | \(1/2\) | \(4\) | \(4\) | \(3.331\times10^{-16}\) |
| \([ -\tfrac12,\tfrac12]^2\) | \(1\) | \(2\) | \(2\) | \(2.220\times10^{-16}\) |
| \([ -\tfrac12,\tfrac12]^2\) | \(3/2\) | \(4\) | \(4\) | \(1.110\times10^{-16}\) |
| \(\operatorname{conv}\{(-1,-1),(2,-1),(-1,2)\}\) | \(1/2\) | \(18\) | \(4\) | \(6.661\times10^{-16}\) |
| 同上 | \(1\) | \(9\) | \(2\) | \(4.714\times10^{-13}\) |
| 同上 | \(3/2\) | \(18\) | \(4\) | \(1.110\times10^{-16}\) |

Wolfram 的五个精确方向上，原始--对偶及反射误差均严格为零；所有 MATLAB 方向上的反射对偶误差至多 \(5.551\times10^{-16}\)。三角形的 \(\beta=1\) 结果确为 \(|S|=9\)、\(|S|/|K|=2\)。这些数值交叉检查验证了实现与有限例子，不替代上面的 LP 证明。

完整机器输出保存在 `results/polytope_wolfram_results.txt` 和 `results/polytope_matlab_results.txt`；脚本与输出的 SHA-256 值由本轮终端记录，便于复算时核对。
