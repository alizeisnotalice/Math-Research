# 第81轮：平坦软阈值的实际列覆盖与外部halo结算

> 后续校正（第82轮）：本轮提出的统一软列预算目标已被新的必要下界排除；原目标少了d₁−(log N)/3这一远大于log N的项。平滑归约及曲边界模型定理仍有效。见 [第82轮完整校正](../r082/allocation82.md)。

所有对数均为自然对数，阈值为1，源列为实际函数 \(0\le G\le F\)，费用为实际 \(J_N(G)=\int\Phi_N(G)\)。本附件的覆盖结论只涉及明确的曲边界模型。

## 1. 二维曲边界上的有限软列预算

沿用第80轮已证明的源和接收条带
\[
 F(s,t)=2st\mathbf1_{[0,1]^2},\quad
 I=[25/24,13/12],\quad E=\{(x,Y(x)-\delta):x\in I,0<\delta<\delta_0\},
 \quad\delta_0=1/96.
\]
这里 \(u(x)=x-\sqrt{x^2-1}\)、\(u(Y(x))=1/(2u(x))\)。该轮已逐项证明
\[
 0<M_SF(x,Y(x)-\delta)-1\le6\delta,\qquad J_2(F)=(1+5\log2)/8.
\tag{1}
\]
并且对任意 \(\epsilon>0\)，取
\[
 M=\left\lceil\frac1{12\sqrt\epsilon}\right\rceil,\quad
 a_k=2/3+k/(12M),\quad b_k=1/(2a_k),\quad k=0,\ldots,M,
\]
实际完整源 \(F\) 的这些角点列覆盖 \(E_\epsilon=E\cap\{\delta\ge\epsilon\}\)。其严格平均余量的未归一化下界为 \(229\epsilon/768>0\)。这是连续参数上的证明；第81轮不以有限取样替代它。

设 \(\chi\) 为 `smoothing81.md` 中两端无限阶平坦的截断。令
\[
 \delta_j=\delta_0 2^{-j},\quad \epsilon_j=\delta_{j+1},\quad
 n_j=\left\lceil\sqrt{2^{j+2}/3}\right\rceil+1,\quad
 w_j=\chi(6\delta_j/\eta).
\tag{2}
\]
给第 \(j\) 个角点族的每列相同常数权重 \(w_j\)。在条带
\(\delta_{j+1}\le\delta\le\delta_j\) 内，该族至少一列实际成功，且由(1)和单调性
\[
 \chi((M_SF-1)/\eta)\le w_j.
\]
所以可数列的成功指示函数加权和覆盖整个 \(E\) 上的软目标。总费用不超过
\[
 B_{\rm flat}(\eta)=J_2(F)\sum_{j\ge0}n_j\chi(6\delta_j/\eta)<\infty.
\tag{3}
\]
注意同一完整源在不同角点重复收费；这里没有把源复制当成免费。

若 \(\eta=2^{-k}\)，令 \(m=j-k+4\)。在 \(m\ge1\) 时
\[
 \chi(2^{-m})\le e^2e^{-2^m}<8\,2^{-2^m},\qquad n_{j+1}\le2n_j.
\tag{4}
\]
因此从 \(J=k+4\) 开始的全部尾项满足有理上界
\[
 \sum_{j\ge J}n_jw_j\le16n_J2^{-256}.
\tag{5}
\]
指数常数由 \(2<e<11/4\)、\((11/4)^2<8\) 推得；代码用正Taylor和及正几何余项给出有理证明。前缀的每个指数函数均由正Taylor区间、反复平方及向外有理舍入包围，最后乘以实际 \(J_2(F)\) 的独立区间，而非只结算源质量。

(2)--(4)也给出 \(B_{\rm flat}(\eta)=O(\eta^{-1/2})\)（\(\eta\downarrow0\)）。另一方面，固定选足够大的常数 \(A\)，在 \(\delta\asymp A\eta\) 的截面，\(M_SF-1\ge(32/15)\delta\) 使软目标等于1。第80轮的截面下界随即给出任何此类软列覆盖的预算 \(\gtrsim\eta^{-1/2}\)。所以此模型的软预算阶为 \(\Theta(\eta^{-1/2})\)。它既非最优预算常数，也非一般 \(C_N\) 上界。

## 2. 全维乘积源的外部halo是正项级数

取
\[
 F_N(s)=2\prod_{i=1}^N s_i\,\mathbf1_{[0,1]^N},\quad T=\log2.
\]
仅考虑所有 \(x_i>1\) 的一个外部角区。逐坐标区间最优化给
\[
 M_SF_N(x)=2\prod_i u(x_i).
\]
令 \(v_i=-\log u(x_i)>0\)，则 \(x_i=\cosh v_i\)、\(dx_i=\sinh(v_i)dv_i\)。所以这个外部halo的体积为
\[
 V_N(T)=\int_{\sum v_i<T}\prod_i\sinh v_i\,dv
 =\sum_{k\ge0}\binom{N+k-1}{k}\frac{T^{2N+2k}}{(2N+2k)!}.
\tag{6}
\]
证明：展开每个 \(\sinh\) 的正幂级数，Tonelli允许交换积分；单纯形Dirichlet积分恰好约去各 \((2k_i+1)!\)，固定 \(\sum k_i=k\) 有 \(\binom{N+k-1}{k}\) 个组合。

归一化级数 \(R_N=V_N(T)(2N)!/T^{2N}\) 首项为1，项比不超过
\[
 r_N=\frac{T^2}{2(2N+1)}<1,\qquad 1\le R_N\le(1-r_N)^{-1}.
\tag{7}
\]
于是
\[
 V_N(T)=\frac{T^{2N}}{(2N)!}(1+O(N^{-1})),\quad
 V_2(T)=1-\cosh T+\tfrac T2\sinh T=(3\log2-2)/8.
\tag{8}
\]
所有项为正，小尾不会由两个接近大数相减产生。

实际Orlicz分母满足
\[
 2^{1-N}\le J_N(F_N)\le2^{1-N}
 \left(1+\frac{(2T)^N}{(N-1)!}\right).
\tag{9}
\]
确实 \(Z=-\log\prod_i s_i\) 的积分密度是 \(e^{-z}z^{N-1}/(N-1)!\)，故
\[
 J_N(F_N)-2^{1-N}
 =\frac2{(N-1)!}\int_0^T e^{-2z}z^{N-1}
 \big[(1+T-z)^{N-1}-1\big]dz.
\]
由 \(T<1\)、\(z(1+T-z)\le T\) 可得(9)。所以
\[
 \log\frac{V_N(T)}{J_N(F_N)}
 \le (N-1)\log2+2N\log T-\log(2N)!+O(N^{-1}).
\tag{10}
\]
这一外部接收角区在大维数极小。第80轮无限列预算及本轮幂次障碍都不等于实际弱型常数发散；此例尤其显示不能用列数代替接收者体积结算。这里没有计算其余内外混合角区的整个halo。

## 3. 该外部halo的平坦阈值损失

记同一外部角区上
\[
 V_N^{\rm flat}(\eta)=\int\chi((M_SF_N-1)/\eta)\,dx.
\]
令 \(f_N=V_N'\)。以 \(h=M_SF_N-1\) 作换元，精确有
\[
 V_N(T)-V_N^{\rm flat}(\eta)
 =\int_0^\eta [1-\chi(h/\eta)]
 \frac{f_N(T-\log(1+h))}{1+h}\,dh,
\tag{11}
\]
其中先取 \(0<\eta<1\)，\(T-\log(1+\eta)>0\)。

正项级数给出，固定 \(T=\log2\)，
\[
 \frac{f_N(T)}{V_N(T)}=\frac{2N}{T}+O(N^{-1}),\quad
 \frac{f_N'(s)}{f_N(s)}=O(N)
 \quad (T/2\le s\le T).
\tag{12}
\]
例如(7)的同类项比界控制 \(R_N'\) 为 \(O(N^{-1})\)，并且所有次数的加权平均为 \(2N+O(N^{-1})\)；再求导同样控制二阶矩，得到(12)。常数在该固定区间一致。

取 \(\eta=N^{-2}\)，则(11)中密度相对 \(f_N(T)\) 为 \(1+O(N\eta)\)。由 \(\chi(z)+\chi(1-z)=1\)，
\(\int_0^1[1-\chi(z)]dz=1/2\)。从而
\[
 \boxed{\quad
 \frac{V_N^{\rm flat}(N^{-2})}{V_N(T)}
 =1-\frac1{N\log2}+O(N^{-2}).\quad}
\tag{13}
\]
这是此明确外部模型的定理；不是一般输入逐函数的相对halo估计。数值代码用正项密度与有限支撑缺额积分复核(13)，没有对极小体积直接做机器精度相除。

## 4. 与一般尖锐上界之间仍缺的非循环估计

令 \(\mathcal W_{N,\eta}(A)\) 为所有 \(F\le A\) 的最优真实角点列软覆盖费用除以 \(J_N(F)\) 后的上确界，覆盖目标为全空间 \(\chi((M_SF-1)/\eta)\)。此定义不预设其有限。

固定角点张量Hardy结算与 \(t^p\le\Phi_N(t)\)（\(t\le A\)）给
\[
 S^+_{N,A}(\eta)\le K_N(A)\mathcal W_{N,\eta}(A),\quad
 K_N(A)=(p_A/(p_A-1))^{Np_A},\quad
 p_A=1+\frac{\log(1+\log A)}{\log A}(N-1).
\tag{14}
\]
同一角点各方向按源正交象限分割，不多付 \(2^N\)。对 \(2\le A\le3\)，\(\log K_N(A)=N+O(1)\)。结合第73轮的一般振幅归约及本轮平滑比较，若能非循环证明
\[
 \log\mathcal W_{N,N^{-2}}(3)
 \le(\mu-1)N-d_N+\frac{9-\mu}{24}\log N+o(\log N),
\tag{15}
\]
就足以将一般上界匹配现有下界的 \(\log N\) 系数。

本轮没有证明(15)，也没有证明对一般输入 \(\mathcal W\) 统一有限。平坦截断消除了特定曲边界的有限阶障碍；它没有提供全输入、全接收者的维数尖锐预算。第80轮的硬指示函数预算不可能定理仍然有效。
