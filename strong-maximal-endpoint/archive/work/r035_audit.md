# r035 独立审查（2026-11-06）

## local all-action gain：桥四阶矩的先前反对已撤回；仍待统一参数引理

1. `d_t=R_t+Q_t/N=-partial_r v_t^a` 的使用是必要且方向正确。尤其
\[
\frac{(bQ_t/N)^2}{d_t}\le b^2Q_t/N\le b^2/N
\]
只用 \(d_t\ge Q_t/N\)，所以没有把下管的 \(d_t\) 过早替换为 \(f_t\)。这一项的估计可接受。

2. 参数 Taylor 的精确恒等式若成立，确实产生文中两个主项；把
\(b f_t+x f_t/(2t)\) 写成 \(f_t(r-t\ell)/(2t)\) 是代数正确的，因为
\(\ell=m(a)-2b\)。但是 (A)--(B) 在当前文字中只是推导提纲，不是可引用的严格引理。要把它用于 all-action 根定位，仍需给出：

   - 对 \(|a-\mu|=O(\sqrt{\log N/N})\) 的统一常数，以及精确的条件矩/平移 Taylor 余项；
   - 由实际倾斜密度导出的 \(R_t/f_t\)、\(f_t'/f_t\) 与二阶导数界，而不是将 \(E_A,E_B\) 当作假设；
   - 在整个指定自然管和紧动作区间上对三阶参数导数的可积控制。

   r033 的 Fourier 局部材料可作为这些引理的来源，但 r035 必须明确其向参数 \(a_N\) 的均匀延拓。没有该延拓，(G) 不能被视作已证明。

3. **更正（先前的量级反对不成立）。** (G) 的中间误差有形式
\(f_t\rho_t\)，而非直接的绝对 payoff 误差。桥/Poisson 传播后总量为
\[
 f_N(q)\sum_{t\ge t_0}\mathbb E[\rho_t(S_t)\mid S_N=q].
\]
四阶矩给出的相对系数是
\[
O\!\left(\sum_{t\ge t_0}t^{-2}+Nb^4\right)=O(1),
\]
因为 \(b=O(\sqrt{\log N/N})\) 蕴含 \(Nb^4=O(\log^2N/N)=o(1)\)。所以这部分确实是 \(O(f_N(q))\)，不需要 \(Nb^4=O(f_N(q))\)。此前审计把这个传播后的 \(f_N(q)\) 因子遗漏，现明确撤回该阻断意见。

   与之不同，(G) 中的常数项 \(C_\varepsilon b^2/N\) 不携带 \(f_t\)，逐时间相加为 \(O(b^2)=O(\log N/N)=O(f_N(q))\)；文中的区分是正确的。

   事实上不必使用条件桥四阶矩来做这一步。取 \(t_0>e\) 固定，在
\(t_0\le t\le N\) 上 \(\log t/t\ge\log N/N\)，故
\(b^2\le C\log t/t\)。局部管内还具有 \(\eta^2\le C\log t/t\)，从而
\[
\left[(|b|+|\eta|)^2+t^{-1}\right]^2
\le C\frac{(\log t)^2}{t^2}.\tag{C}
\]
把该**局部**中间源在管外置零（或在拼接中由外侧屏障单独吸收），其全局上包络为
\(Cf_t(\log t)^2/t^2\)。参考半群恒等式
\(P_a^{N-t}f_t=f_N\) 给
\[
\sum_{t=t_0}^N P_a^{N-t}
 \left[f_t C\frac{(\log t)^2}{t^2}\right](q)
\le Cf_N(q)\sum_{t=t_0}^{\infty}\frac{(\log t)^2}{t^2}
=O(f_N(q)).\tag{D}
\]
这比四阶矩桥论证更直接。它不处理管外的实际动作增益；该部分仍必须由已证明的全局屏障或单独源项控制，不能仅把局部 (G) 延拓到管外。

4. (G) 的首项仅对每个固定 \(\varepsilon>0\) 有 \((1+\varepsilon)p_a/8\)，且常数和开始时间可依赖 \(\varepsilon\)。它支持“先固定 \(\varepsilon\)，再令 \(N\to\infty\)”的非一致一阶上界；不能在没有 \(\varepsilon\)-定量余项的情况下令 \(\varepsilon=\varepsilon_N\downarrow0\)，也不能由此宣称锐系数已经以统一 all-action 方式闭合。

结论：在实际桥四阶矩估计已经给出的前提下，reference-summability 的量级链条正确。局部结构、全局根定位的目标和 `d_t` 分母修复都合理；全局 Bellman MAX 校正仍未完成。尚需补写的是参数 \(a_N\) 的统一实际 Taylor/密度导数引理，而不是 \(Nb^4\) 的额外消去。

## tensor lift：分母计算通过，障碍的量词需保持一阶范围

对 \(G_{k,\varepsilon}=\prod_i(1+\varepsilon a_i)\)，高块数
\(J\sim\mathrm{Binomial}(k,15/16)\)，所以 (2) 是实际 \(\Phi_n\) 分母而非错误的二维分母乘积。直接微分给
\[
 D'_{n,k}(A,0)=A(1+\log A)^{n-1}k\alpha
 \left(1+\frac{n-1}{1+\log A}\right),
\]
且对一般 \(H\ge1\)，(4) 等于 \(k\alpha\int H\Phi_n'(H)\)，与写出的式子相同。因而有界基线的**相对一阶**价格确为 \(\Theta(nk\varepsilon)\)。

这足以否定把二维 \(\Phi_2\) 比值相乘当作高维真实分母的做法。它也说明直接 product witness 的可见 \(O(k\varepsilon)\) 体积变化不能在一阶与分母竞争。

不过，`exp(Theta(n k epsilon))` 的叙述应限于已经给出统一导数控制的 epsilon 范围，或改称“在 epsilon=0 的对数一阶导数为 Theta(nk)”。仅凭 (3)--(4) 的点导数，不能对任意随 n,k 变化的 \(\varepsilon\) 推出同样的指数比较。纯立方体的 r024 ceiling 在 \(k\varepsilon=o(1)\) 下提供了独立的严格式界，故这部分范围正确。

该张量障碍仍只针对有界基线、直接积输入与同时积见证；文档已正确保留薄双曲支撑、非积接收并集和大对数振幅为开放机制。
