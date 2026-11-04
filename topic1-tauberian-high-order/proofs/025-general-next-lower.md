# 各维对称剪切构造的下一项

日期 2026-11-05。固定整数 \(n\ge2\)，置 \(u=t^{1/n}\)、\(v=2u\)、\(a=v^n\)、
\(g=(1+2u)^{-1}\)、\(m_*=ag^n\)、\(L_n(t)=(1-m_*)^{-1}\)，记
\[
 B_n=\binom n2,\qquad K_n=B_n n^2/27 .
\]

**命题。** 对所有充分小的实数 \(t>0\)，存在有限正测度输入 \(E_t\)，使其全部半径严格水平集的实际比值满足
\[
 \frac{|\mathcal H_n(1-t;E_t)|}{|E_t|}
 =L_n(t)+K_n\,a v^6 g^{n+4}+O_n(a v^8)
 =L_n(t)+K_n\,a v^6-(n+4)K_n\,a v^7+O_n(a v^8).                     \tag{1}
\]
特别地，存在 \(C_n'<\infty,t_n>0\)，使
\[
 C_n(1-t)\ \ge\ L_n(t)+K_n a v^6-(n+4)K_n a v^7-C_n' a v^8
 \qquad(0<t<t_n),                                                       \tag{2}
\]
其中 \(C_n'\) 是只依赖维数的余项常数。(1) 展开具体输入族，(2) 给真实 Tauberian 常数的一侧界；确定其七阶系数还需相应的全输入估计。

## 1　孔质量和参数

用 \(J\) 表示全一矩阵。定义显式有理函数和维数常数
\[
 \bar\kappa=\frac{8n}{3},\quad \kappa=\bar\kappa g^2,\quad
 Z_n=\frac{3^n}{4}\binom{n-1}{2}\bar\kappa^4,\quad c_n=1+Z_n/n ,
\]
其中 \(Z_2=0\)。取
\[
 e=\kappa u^3,\quad T=(1-e)I+eJ,\quad d=\det T,\qquad
 \Lambda_*=\frac{B_n\kappa^3}{8n^2 g},\quad \Lambda=\Lambda_*+c_nu^2,
\]
\[
 w=ug+\Lambda u^6,\qquad \rho=1-2w=g-2\Lambda u^6,\qquad M=a\rho^n.   \tag{3}
\]
小 \(u\) 时 \(w/u\to1\)，\(\rho,d\to1\)，\(e/w\to0\)。实际孔质量定为 \(M\)，外包半宽为 \(w\)。

## 2　所有邻孔与坏点的三面层

令 \(C_w=[-w,w]^n\)。对 \(x\in C_w\)，写 \(s=w-\|x\|_\infty\)、\(r=\rho+s\)、\(h_{i,\sigma}=w-\sigma x_i\)。见证 \(Q(x,r)\) 包含自身整个外包，且 \(r+w=1-\|x\|_\infty\)。由 \(\|Tk\|_\infty\ge(1-(n-1)e)\|k\|_\infty\)，外包孔两两分离；\(\|k\|_\infty\ge2\) 时见证不与 \(Tk+C_w\) 相交。对 \(0\ne k\in\{-1,0,1\}^n\)，若 \(\sum_i k_i\ne0\)，取与该和同号的坐标，其中心绝对值至少一，故无正交量。若 \(\sum_i k_i=0\)，则 \(Tk=(1-e)k\)，支撑大小 \(q(k)\) 为偶数。支撑轴的交长为 \((e+s-h_{i,k_i})_+\le e\)，其余轴为整宽 \(2w\)。

因此邻孔外包总交量恰为 \((2w)^{n-2}R(x)\)，其中
\[
 R(x)=\sum_{\substack{0\ne k\in\{-1,0,1\}^n\\ \sum k_i=0}}
   (2w)^{2-q(k)}\prod_{k_i\ne0}(e+s-h_{i,k_i})_+ .
\]
记
\[
 K=\frac a{(2w)^{n-2}},\quad F(z)=K[(\rho+z)^n-\rho^n],\quad
 G=\{x\in C_w:s>0,\ R(x)<F(s)\}.
\]
令严格单调的 \(F\) 的正根 \(h_*,h_H\) 分别满足 \(F(h_*)=e^2,\ F(h_H)=3^n e^2\)。由于
\[
 F'(0)=nK\rho^{n-1}\sim4nu^2,\qquad
 h_H\le\frac{3^n e^2}{F'(0)}
      =\left(\frac{3^n\bar\kappa^2}{4n}+o_n(1)\right)u^4,             \tag{4}
\]
有 \(h_*\le h_H=O_n(u^4)\)、\(e+h_H<w\)。每一项 \(R(x)\) 都不超过 \(e^2(e/(2w))^{q-2}\le e^2\)，故任何内点坏点都有 \(s\le h_H\)。

令 \(\mathcal T\) 由满足 \(s\le h_H\) 且至少三个不同坐标各有一个符号面深度 \(\le e+h_H\) 的点组成。全局最小面深 \(s\) 自己就是一个这样的近面。选其坐标、符号，再选另外两个近面，得到 \(n\ge3\) 时
\[
 |\mathcal T|\le 8 n\binom{n-1}{2}\,h_H(e+h_H)^2(2w)^{n-3},
 \qquad \frac{|\mathcal T|}{a}\le (Z_n+o_n(1))u^7 .                    \tag{5}
\]
后式用 (4)、\((e+h_H)/u^3\to\bar\kappa\) 和 \((2w)^{n-3}/a\sim1/(8u^3)\)。\(n=2\) 时取 \(\mathcal T=\varnothing\)。即使达到最小面的坐标不在最初列出的三个近面中，也可改选该坐标并保留另外两个近面，(5) 的覆盖仍然成立。

## 3　单角面积与严格预算

在 \(\mathcal T\) 外的内点坏点，有且仅有两个坐标以相反符号组成进入的有向角。确实，高支撑至少需要四个近面，不同的有向角至少用到三个坐标，同一坐标的两个符号因 \(e+h_H<w\) 互斥。既然 \(R(x)\ge F(s)>0\)，确有一角进入。全局最浅面由该角的两个坐标之一取得，否则又有第三个近面。

将该角的两个面深写为 \(z,z+y,\ y\ge0\)，则 \(s=z,\ R=e(e-y)_+\)。两个排序合计的二维面积为
\[
 \mathscr A=2\int_0^{h_*}\left(e-\frac{F(z)}e\right)\,dz .
\]
每个坐标对有两个相反的有向角。在另外 \(n-2\) 轴取全宽，并计入 (5)，故
\[
 |C_w\setminus G|\le 2B_n(2w)^{n-2}\mathscr A+|\mathcal T|,\qquad
 \mathscr A\le\frac{e^3}{F'(0)} .                                      \tag{6}
\]
面积上界直接将 \(F(z)\) 以其在零点的切线下界代入，再积分到 \(e^2/F'(0)\)。这里 (6) 只需要上界。若要核对单角渐近式，\(0\le z\le h_H\) 上有 \(F'(z)/F'(0)=1+O_n(u^4)\)，故用上下两条直线夹住 \(F\) 得 \(\mathscr A=e^3/F'(0)[1+O_n(u^4)]\)。一个单角模型点必有该角较浅面 \(\le h_*\)，两面均 \(\le e+h_*\)。如出现另外一个更浅面或第二条有正交量的邻孔，必落入 \(\mathcal T\)。因此在 \(\mathcal T\) 外真实内点坏集与所有单角模型区域的并完全相同，且这些区域互不交；双侧比较仅多 \(O_n(|\mathcal T|)\)。下述可行性只用带显式三面系数的上侧 (6)。

将单角的切线上界代入并保留 \(g\)，有精确代数式
\[
 \frac{2(2w)^{n-2}}a\,\frac{e^3}{F'(0)}
 =\frac{\kappa^3}{8n}u^5
  \frac{(g+\Lambda u^5)^{2n-4}}{(g-2\Lambda u^6)^{n-1}}
 =\frac{\kappa^3}{8n}u^5 g^{n-3}[1+O_n(u^5)] .                       \tag{7}
\]
实际可删除量 \(D=(2w)^n-M\) 则满足
\[
 \frac Da=(g+\Lambda u^5)^n-(g-2\Lambda u^6)^n
   =n\Lambda u^5 g^{n-2}+O_n(u^{10}),                                  \tag{8}
\]
因为一阶项为 \(n g^{n-1}\Lambda u^5(1+2u)\)。所有参数由 (3) 显式固定，故余项常数只依赖 \(n\)。将 (7)、(8) 与 (5)、(6) 合并，\(\Lambda_*\) 的项精确抵消，剩下
\[
 \frac{D-|C_w\setminus G|}{a}
 \ge [n c_n g^{n-2}-Z_n-o_n(1)]u^7
 =(n-o_n(1))u^7>0 .                                                     \tag{9}
\]
这一步只用 \(\rho=1-2w\)，未要求 \(w=u\rho\)。

## 4　真实孔形、有限化和下一项

按 (9)，\(|G|>M\)。函数 \(z\mapsto|G\cap\{x_1<z\}|\) 连续，故可用一个坐标切片选出 Borel 集 \(P\subset G\) 且 \(|P|=M\)。记 \(\mathscr D=\bigcup_{k\in\mathbb Z^n}(Tk+P)\)。对每个 \(x\in P\)，
\[
 |\mathscr D\cap Q(x,r)|\le M+(2w)^{n-2}R(x)<a r^n=t(2r)^n .
\]
这是包含全部自身质量和邻孔的严格水平比较。

取 \(N=\lceil u^{-9}\rceil,\ \Omega_N=T([-N,N]^n),\ E_t=\Omega_N\setminus\mathscr D\)。按周期积分，\(|E_t|=(2N)^n(d-M)>0\)。小 \(u\) 时 \(\|T^{-1}\|_{\infty\to\infty}<2\)，且 \(w+r\le1\)。因而 \(\|k\|_\infty\le N-3\) 的每个孔及其所有上述指定见证均在 \(\Omega_N\) 内。连同输入密度点，这给出下界 \(1+(2N-5)^n M/[(2N)^n(d-M)]\)。凸体 \(\Omega_N\) 外点的任一中心立方体至多一半体积在 \(\Omega_N\) 所在的分离半空间内，故在 \(1-t>1/2\) 时整个严格水平集在 \(\Omega_N\) 内。因此真实有限比值为
\[
 \frac{|\mathcal H_n(1-t;E_t)|}{|E_t|}=\frac d{d-M}+O_n(M/N)
                                    =\frac d{d-M}+O_n(a u^9).          \tag{10}
\]

最后 \(d=(1-e)^{n-1}(1+(n-1)e)=1-B_n\kappa^2u^6+O_n(u^9)\)，而 \(M/a=(g-2\Lambda u^6)^n\)。相减可得
\[
 \frac{M-dm_*}{a}
 =B_n\kappa^2g^n\left(1-\frac{\kappa}{4n g^2}\right)u^6+O_n(u^8)
 =64K_n g^{n+4}u^6+O_n(u^8).                                          \tag{11}
\]
其中 \(c_nu^2\) 对质量根的改变量为 \(-2c_nu^8\)。利用
\[
 \frac d{d-M}-L_n(t)
 =\frac{M-dm_*}{(d-M)(1-m_*)},
\]
分母为 \(1+O_n(u^2)\)，因为 \(n\ge2\)。即使在 \(n=2\) 其修正也只从 \(a u^8\) 起。由 \(v=2u,\ g=(1+v)^{-1}\)，(10)、(11) 给出 (1)、(2)。二维的七阶系数为 \(-6(4/27)=-8/9\)，与已审 R024 的另一有限输入族一致。
