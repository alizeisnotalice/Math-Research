# r036：对数凹 Doob 核结构审计（2026-11-06）

## 已确认的结构引理

令
\[
\phi_N(r)=\begin{cases}e^{-r},&r>0,\\1-r/N,&r\le0.\end{cases}
\]
它在全实轴正且对数凹：左支 \(\log(1-r/N)\) 严格凹，右支线性，且在零点对数导数从 \(-1/N\) 向右跳到 \(-1\)，是允许的向下跳跃。虽然 \(\phi_N\) 在左端不属于 \(L^1\)，这不妨碍 Prékopa 边缘化：对 \(\phi_N e^{-\epsilon|r|}\) 先应用通常的可积版本，再令 \(\epsilon\downarrow0\)；卷积在每个点有限，因为增量有有限一阶矩。故若共同参数密度 \(f_{t,a}\) 对数凹（r027 对 \(t\ge4\) 的结论），则
\[
v_{t,N}^a=\phi_N*f_{t,a}
\]
正且对数凹。它也单调递减。因此，对未来函数
\(v=v_{t-1,N}^a\)（\(t\ge5\)），
\(\delta(r)=-(\log v)'(r)=-v'(r)/v(r)\) 非减，且凹性的切线不等式给
\[
\frac{v(r-z)}{v(r)}\le e^{\delta(r)z},\qquad z\ge0.\tag{1}
\]

在自然管中，r029 的共同参数、带符号倾斜密度估计给
\(R_t\le Cf_t\)。另外，局部 Mills 比率给
\[
f_t(r)/Q_t(r)\le C\sqrt{\log t/t}.
\]
这在上管来自 saddlepoint/Mills，在下管由 \(Q_t\) 有正下界和 \(f_t=O(t^{-1/2})\)；常数须统一于 \(a\) 的固定紧邻域和 \(N\ge t\)。由于 \(v\ge Q_t\)，
\[
\delta(r)=\frac{R_t+Q_t/N}{v_t}
\le C\sqrt{\log t/t}.\tag{2}
\]
这里 (2) 不是单靠抽象对数凹性得到的，必须实际调用上述统一局部密度/Mills 引理。

对 Doob 变换核
\[
\widehat q_c(z;r)=q_c(z)\frac{v(r-z)}{K_cv(r)},
\]
单调性给 \(K_cv(r)\ge v(r)\)，而 \(q_c(z)\le2e^{-z}\)。由 (1)，当 \(\delta\le1/2\)，
\[
\widehat q_c(z;r)\le2e^{-(1-\delta)z}\le2e^{-z/2}.\tag{3}
\]
该界对所有动作 \(c\) 一致，且与 r032 的不同参数之间不存在 FOSD 完全不同。

## 加权 Lipschitz 到 Jensen 误差的条件

若一个一致有界函数 \(w\) 满足
\[
|w(r-z)-w(r)|\le Lz e^{\kappa z},\tag{4}
\]
则 (3) 能推出 \(\log\widehat{\mathbb E}e^w-\widehat{\mathbb E}w\le C L^2\)，但必须补充
\[
2\kappa<1-\delta
\quad\text{（例如在 \(\delta\le1/2\) 时取 \(\kappa<1/4\)）},\tag{5}
\]
以及 \(\|w\|_\infty\) 的统一界。事实上 (3)--(5) 先给
\(\operatorname{Var}_{\widehat q}w\le CL^2\)；有界性再把指数 Jensen 差由方差控制。若 \(\kappa\) 未受限，则 \(z^2e^{2\kappa z}\) 在 (3) 下可能不可积；仅称“bounded accessible”而不写出这两个条件不足以推出 \(CL^2\)。

因此，若后续确实证明所有需要的全局状态移位函数满足 (4)--(5)，局部导数尺度
\(L\asymp(f_t/v_t)(\log t)^{3/2}/\sqrt t\) 所产生的源项
\(Cf_t(\log t)^4/t^{3/2}\) 可由参考半群求和。当前“全局 \(w\) 移位界”尚未给出，故这一段只能登记为条件接口，不能作为全控制 Bellman 闭合。

## 补充审查：二次桥校正的全局状态移位

令 \(p=(1+e^a)^{-1}\)，\(f_t\) 为
\(\Gamma(t,1)+a\operatorname{Bin}(t,p)\) 的密度，且 \(C_t\) 为 r032 的条件桥系数。逐项使用
\[
{t\choose k}k(t-k)=t(t-1){t-2\choose k-1}
\]
及 Gamma 密度的幂次，可得下列恒等式（第二项的平移不可漏掉）：
\[
\boxed{\quad
 f_t(r)C_t(r)=t\,g_{\Gamma(t+2)+a\operatorname{Bin}(t)}(r)
 +a^2tp(1-p)\,g_{\Gamma(t)+a\operatorname{Bin}(t-2)}(r-a).
 \quad}\tag{6}
\]
因此题设的辅助密度恒等式正确。

设当前 \(r\) 在自然管，并以其共同倾斜 \(\eta\) 为基准。对 \(|\eta|\le1/32\)，原密度和 (6) 中两项在同一倾斜下的归一化常数之比均为 \(O(1)\)，其均值仅相差 \(O(1)\)。若 r029 的 Fourier 论证确实提供这些 Gamma--Binomial 混合密度及有限阶导数的**全局上确界**（在倾斜均值处有相应下界），则
\[
|U'(r-z)|\le
C f_t(r)\frac{(\log t)^{3/2}}{\sqrt t}e^{\kappa z},
\qquad
U(r-z)\le C f_t(r)(\log t)e^{\kappa z},\tag{7}
\]
对 \(z\ge0\) 成立；多项式因子可吸收到任意固定
\(\kappa>1/32\)，例如 \(\kappa=1/16\)。结合 \(v(r-z)\ge v(r)\)、
\(\delta(r-z)\le\delta(r)\)，给出
\[
|w'(r-z)|\le C\frac{f_t(r)}{v_t(r)}
\frac{(\log t)^{3/2}}{\sqrt t}e^{\kappa z}.\tag{8}
\]

可达状态上的 \(w=U/v\) 也有统一小上界：右半段用
\(f/v\le C\sqrt{\log t/t}\) 和 \((r-t\ell)^2/t=O(\log t)\)，给出
\(O((\log t)^{3/2}/\sqrt t)\)；左半段用密度的全局下尾界、\(Q\ge c\) 及 \(C_t/t=O(1)\)，给出较小的 \(O(\log t/\sqrt t)\)。故
\[
\|w\|_\infty\le C(\log t)^{3/2}/\sqrt t=o(1).
\]
由 (8)、\(\kappa=1/16<1/4\) 和前节 (3)--(5)，严格 Jensen 模块遂可用，且其单步 defect 至多
\[
C f_t(r)\frac{(\log t)^4}{t^{3/2}},\tag{9}
\]
这是参考半群可求和源项。

**更正：不需要逐点鞍点或两段拆分。** 设
\(a_0=tm(a)-r'\in[0,tm(a)]\)。紧参数集上一致有
\[
\log M_a(-s)\le-m(a)s+Cs^2\qquad(0\le s\le s_0).
\]
取 \(s=a_0/(2Ct)\)，并令 \(C\ge\sup_am(a)/(2s_0)\)，则该 \(s\) 总留在固定紧区间。记负倾斜密度为 \(g_t^{-s}\)。它是
\(\Gamma(t,1+s)\) 与有界 Bernoulli 位移的有限混合，故
\(\|g_t^{-s}\|_\infty\le C/\sqrt t\) 一致成立。恒等式
\[
f_t(r')=M_a(-s)^t e^{sr'}g_t^{-s}(r')
\]
于是直接给出
\[
f_t(r')\le Ct^{-1/2}\exp\{-a_0^2/(4Ct)\},\tag{10}
\]
覆盖全部 \(0\le r'\le tm(a)\)，包括 \(r'/t\to0\)。此前审计要求两段证明过强，现撤回。

这只完成参考 Doob/Jensen defect 的控制；全控制 Bellman transport、右尾延拓及拼接仍是独立问题。
