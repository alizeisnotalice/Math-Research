# r040：右屏幕与 corrected inner 的 crossing 审计

日期：2026-11-06。范围仅为无限终端、固定 \(\mu\) 参考的固定状态比较；不声称
full transport 或全局 Bellman supersolution。

令 \(\alpha=1+\varepsilon\)，令 old inner 为
\(V_t=v_t\exp[(\alpha J_t+s_tf_t)/v_t]\)，并以
\(\widetilde V_t=V_t\exp[tF(\eta)]\) 修正，其中
\[
F=I_\mu-I_*-D=K_4\eta^4+O(\eta^5),\qquad
D=A_3\eta^3+O(\eta^4).
\]
取
\[
\eta=C(\log t/t)^{1/3},\quad y=\Lambda_\mu'(\eta),\quad r=ty.
\]
则可严格得到
\[
\log\frac{\widetilde V_t(ty)}{O_t^+(ty)}
=\left(\varepsilon A_3C^3-\frac16\right)\log t
-\frac13\log\log t-\log C-rac12\log(2\pi\mu)+o(1).\tag{1}
\]

理由如下。固定小正倾斜的 saddlepoint density 与 uniform Mills 主项给
\[
 f_t(ty)=\frac{e^{-tI_\mu(y)}}{\sqrt{2\pi t\mu}}(1+o(1)),
\qquad {v_t(ty)\over f_t(ty)}={1\over\eta}(1+o(1)).\tag{2}
\]
这里所需的相对误差足够弱：
\[
{v\over f}={1\over\eta(1-\eta)}\left[1+O\left((t\eta^2)^{-1}+(t\eta)^{-1}+t^{-1}\right)\right].\tag{3}
\]
于是 \(J_t/v_t=tD+o(1)\)：调和源贡献为 \(O(\eta\log t)=o(1)\)，
Mills 相对误差产生 \(tD/(t\eta^2)=O(\eta)=o(1)\)，有界收敛 source 给
\(s_tf/v=O(\eta)=o(1)\)。又 \(tF=O(t\eta^4)=o(1)\)。代入
\(\log v=-tI_\mu-\frac12\log(2\pi t\mu)-\log\eta+o(1)\)，并利用
\(O^+=e^{-tI_*}\)，即得 (1)。

故临界常数为
\[
C_{\rm crit}=(6\varepsilon A_3)^{-1/3}.
\]
当 \(C<C_{crit}\) 时 inner/屏幕趋零；当 \(C>C_{crit}\) 时趋无穷；在等号处，
\(-\frac13\log\log t\) 仍使 inner/屏幕趋零。这精确说明 screen 在何处开始
选择，并不涉及自然边界的 \(\eta/t\) transport 精度。

必须使用 uniform Mills 而非 fixed-\(\eta\) 的 \(O(\log t)\) 截断；主积分位移是
\(u\asymp1/\eta\). WolframScript 15 实测正态 Mills 代理
\(\int_0^\infty e^{-\eta u-u^2/(2\mu t)}du\)；在上述临界 \(C\) 下，
\((Q/(f/\eta)-1)t\eta^2\) 于 \(t=10^6,10^9,10^{12}\) 分别为
\(-.7821269,-.7821836,-.7821879\)，收敛到 \(-1/\mu\)，与 (3) 一致。

## 一致 Mills 主项的充分证明

以 \(g_t\) 记 \(\eta\)-倾斜后、移至均值的 \(t\) 步密度。对
\(0\le\eta\le\delta<1\)，Exp 因子和紧参数 Bernoulli 因子给
\(|\psi_\eta(u)|\le e^{-cu^2}\) 于原点邻域，\(|\psi_\eta|\le\rho<1\)
于其外，且 \(\int u^2|\psi_\eta(u)|^mdu\) 一致有界。Fourier 反演给
\[
 g_t(0)\ge c t^{-1/2},\quad \|g_t''\|_\infty\le Ct^{-3/2},
 \quad g_t'(0)=O(t^{-3/2}).
\]
最后一项须用高斯主项的奇性消失和三阶 cumulant 项（或中心 Edgeworth 展开）；
绝对值 Fourier 估计本身只会给较弱界。Taylor 定理遂对所有 \(u\) 给
\[
 |g_t(u)/g_t(0)-1|\le C(|u|+u^2)/t.
\]
对 \(v/f\) 的两端权 \(e^{\eta z}\)、\(e^{-(1-\eta)z}\) 积分，即得到 (3)，
其中相对式只在 \(0<\eta\le\delta\) 主张；\(\eta=0\) 时 \(h_0\) 发散，不能按
该式理解。这提供 crossing 所需的主项，而不使用错误的 \(O(\log t)\) 尾截断。
