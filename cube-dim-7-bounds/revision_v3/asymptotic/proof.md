# 一般维数公式族的渐近、导数与有限参数障碍

本报告只分析已证明上界的数值族，不推导最优算子常数的下界。记

\[
\ell_t(\rho)=1+\frac1{2t}-\frac1{2\rho},\quad
h(t)=\frac{1+1/t}{2},\quad
I_n(t)=\int_1^t\ell_t(\rho)^{n-1}\frac{d\rho}{\rho},
\]
\[
A_n(t)=(1+1/t)^n,\qquad
B_n(t)=1-h(t)^n+n\log t-\frac n2 I_n(t),\qquad
F_n(t)=(\sqrt{A_n(t)}+\sqrt{B_n(t)})^2,\quad
U_n=\inf_{t>1}F_n(t).
\]

所有对数均为自然对数。以下结论为解析证明；脚本输出只作独立代数核验。

## 1. 对所有 t 的一致积分估计

对每个整数 n≥2 及每个 t>1，有

\[
0\le I_n(t)\le\log_+(t/n)+8+\log2.
\tag{1}
\]

证明：在 1≤ρ≤t/2 上，用 1−z≤exp(−z) 和 n−1≥n/2 得
\[
\ell_t(\rho)^{n-1}
\le\exp[-(n-1)/(4\rho)]\le\exp[-n/(8\rho)].
\]
在该区间中 ρ≤n 的部分，由 exp(−s)≤1/s 得
\[
\int_1^{\min(n,t/2)}e^{-n/(8\rho)}\frac{d\rho}{\rho}
\le\int_0^n\frac8n\,d\rho=8.
\]
上式在上端小于下端时按空区间理解。ρ≥n 的部分不超过 log₊(t/n)。余下区间 [max(1,t/2),t] 的对数长度不超过 log2，且被积幂不超过1。此证明同时涵盖 t↓1、t≪n、t≈n 以及 t≫n，不依赖 t/n 的预设范围。

## 2. 主渐近式及显式起始维数

令 K=(8+log2)/2、c=log2、M=4/c。对所有整数 n≥16，

\[
n\log n-n\log\log n-(K+\log M)n
\le U_n
\le n\log n-n\log\log n+(4+\log2)n.
\tag{2}
\]

上界：取 t=2n/log n>1。由 log(1+s)≤s，A_n(t)≤√n。又因 I_n≥0 及 h^n≥0，
\[
B_n(t)\le1+n\log t
=1+n\log n-n\log\log n+n\log2
\le1+n\log n.
\]
最后一步用 log n≥2。n≥16 时 log n≤√n（由导数可知此不等式此后保持），因此
\[
\sqrt{A_nB_n}\le\sqrt{n^2+\sqrt n}\le n+1.
\]
将 F=A+B+2√(AB) 展开，利用 √n+2n+3≤4n，得到(2)右侧。

下界：1−h^n≥0 及(1)给出
\[
B_n(t)\ge n\log t-\frac n2\log_+(t/n)-Kn.
\tag{3}
\]
若 t≥n，则 B≥n log n−Kn，已足够。若 1<t≤n，置 x=n/t∈[1,n)。凹函数 log(1+s) 在 s∈[0,1] 上不低于端点弦 cs，故
\[
A_n(t)\ge e^{cx},\qquad
F_n(t)\ge A_n(t)+B_n(t)
\ge n\log n+e^{cx}-n\log x-Kn.
\]
若 x≤M log n，则 log x≤log log n+log M，舍去正的指数项即可。若 x>M log n，则 e^{cx}>n⁴≥n log n≥n log x，所得下界更强。两种情况均给出(2)左侧。

因此严格成立
\[
U_n=n\log n-n\log\log n+O(n).
\]
由论文已经证明的 C_n≤U_n，只能推出
\[
\limsup_{n\to\infty}\frac{C_n}{n\log n}\le1.
\]
式(2)左侧是公式族下确界 U_n 的下界，绝不是 C_n 的下界。这里未证明 C_n 与 U_n 渐近等价。

### 2a. 强化为 o(n) 余项

还可严格证明
\[
U_n=n\log n-n\log\log n+
O\left(\frac{n\log\log n}{\log n}\right).
\tag{2a}
\]
以下给出显式充分范围 n≥⌈exp(32)⌉（此范围只为简化证明，未作优化）。令 L=log n，T=nL−n log L。首先，n≥2 时
\[
I_n(t)\le\frac{2t}{n-1}.
\tag{2b}
\]
事实上 ell^(n−1)≤exp[−(n−1)(1/ρ−1/t)/2]，置 v=1/ρ−1/t，有
\[
I_n(t)\le\int_0^{1-1/t}\frac{e^{-(n-1)v/2}}{v+1/t}\,dv
\le t\int_0^\infty e^{-(n-1)v/2}\,dv.
\]
由第4节的 B′>0，记 x=n/t、x₀=L+4 log L≤2L，分三段：

* x≤L/2 时，B(t)≥B(2n/L)≥T+n log2−4n/L≥T。
* L/2≤x≤x₀ 时，(2b)及 1−h^n≥0 给出
  B(t)≥nL−n log x−2n/x≥T−4n(log L+1)/L。
  这里用 n/(n−1)≤2 和 log(1+4 log L/L)≤4 log L/L。
* x≥x₀ 时，由 log(1+s)≥s−s²/2（s≥0），
  A(t)≥(1+x₀/n)^n≥exp(x₀−x₀²/(2n))≥nL⁴/e≥nL≥T。
  所用 x₀²≤4L²≤2exp L 在 L≥32 时成立。

所以 U_n≥T−4n(log L+1)/L。上界取 t=n/(L−3 log L)，此时 L−3 log L≥L/2，且
\[
A\le n/L^3,\quad
B\le1+n\log t\le T+6n\log L/L+1,\quad B\le2nL.
\]
其中使用 −log(1−y)≤2y（0≤y≤1/2）。因此
\[
U_n\le F_n(t)\le T+\frac{n(6\log L+5)}L.
\]
以上显式两侧估计证明(2a)，并说明 U_n−T=o(n)。它仍然只描述公式族，不给出 C_n 下界；更精细的主导余项及最优参数渐近未在此证明。

## 3. n=1 的端点退化

I₁(t)=log t，且
\[
A_1(t)=1+1/t,\qquad
B_1(t)=\tfrac12(1-1/t+\log t).
\]
对 t>1，由 log t>1−1/t 得
\[
A_1(t)+B_1(t)=\tfrac32+\frac1{2t}+\tfrac12\log t>2.
\]
所以 F₁(t)>2；而 t↓1 时 A₁→2、B₁→0，故 U₁=2。下确界仅在开区间外的端点极限取得，无 t>1 达到它。

## 4. 解析导数及优化证明的边界

对 n≥2，令 J_n(t)=∫₁ᵗ ell_t(ρ)^(n−2)dρ/ρ。Leibniz 法则给出
\[
I_n'(t)=\frac1t-\frac{n-1}{2t^2}J_n(t),
\]
\[
B_n'(t)=\frac n{2t}+\frac n{2t^2}h(t)^{n-1}
+\frac{n(n-1)}{4t^2}J_n(t)>0.
\tag{4}
\]
n=1 时单独有 B₁′=1/(2t)+1/(2t²)，不引入无用的 J₁。又
\[
A_n'(t)=-\frac n{t^2}(1+1/t)^{n-1},
\]
\[
F_n'(t)=(\sqrt A+\sqrt B)
\left[\frac{B_n'}{\sqrt B}
-\frac n{t^2}(1+1/t)^{n/2-1}\right].
\tag{5}
\]
因此内点极小值须满足方括号为零；B′>0 并不保证 F 单调或驻点唯一。尤其 t↓1 时 B_n(t)=n(t−1)+O((t−1)²)，所以 F_n′(t)→+∞。本报告未证明各 n 的全局极小点存在、唯一或某数值驻点的全局性。对固定 n，t→∞ 时 B_n(t)≥1−h(t)^n+(n/2)log t→∞，故若存在 F_n(t)<2^n 的点，则下确界必在内点取得。

## 5. 有限 a 无法改善同一预算族

在论文参数 a>2t、α=a/(a+1)、c=(a(1+1/t)−1)/(a+1) 下，令 u=1/a、L=log t，并记有限参数预算
\[
A_{n,a}=(1+1/t+u)^n,\qquad
E_{n,a}=u^n+\alpha^{-n}J_a(t).
\]
用 ell_a=α ell_t 及 (1+c−α)/(2α)=h(t) 得精确恒等式
\[
E_{n,a}-B_n
=u^n+[(1+u)^n-1](1+nL)-\frac n2uI_n.
\tag{6}
\]
I_n≤L 及 Bernoulli 不等式 (1+u)^n−1≥nu 给出
\[
E_{n,a}-B_n\ge u^n+nu[1+(n-\tfrac12)L]>0.
\]
又 A_{n,a}>A_n。因此对固定 t，每个有限 a 的 τ 优化值
\[
(\sqrt{A_{n,a}}+\sqrt{E_{n,a}})^2
>(\sqrt{A_n}+\sqrt{B_n})^2.
\]
这排除了在现有几何与预算估计保持不变时，仅选择有限 a 来改善最终上界的路线。它不排除改变选择算法或预算后的改进。联合下确界仍满足 inf_{t>1,a>2t} F_{n,a}(t)=U_n：一侧来自逐点支配，另一侧由每个固定 t 下令 a→∞ 得到；不能把逐点严格不等式误报为两个下确界严格不等。

## 6. 复现与状态

运行 `verify_symbolic.py`，使用已安装的 math-research-tools Python/SymPy 1.14.0；另使用本机 MATLAB R2026a Symbolic Math Toolbox 作独立检查。`verification_receipt.json` 记录实际运行结果及必要源文件 SHA-256。所有符号检查均针对有限维数范围的恒等式；任意维数的证明是上述解析论证。
