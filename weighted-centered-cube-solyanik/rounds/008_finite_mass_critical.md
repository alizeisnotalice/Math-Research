# 第八轮独立核验：有限总质量的内禀临界发散

## 结论（通过）

令 \(m\ge2\)、\(a_1,\ldots,a_m>0\)，并定义

\[
P_{\boldsymbol a}(s)=\max\left\{
\prod_{i=1}^m(1-d_i)^{a_i}:d_i\in[0,1),\ \sum_i d_i=s
\right\}.
\]

若 \(\mu=w\,dx\) 为非零有限的 Radon 测度，则对每个
\(0<\lambda\le P_{\boldsymbol a}(1)\)，内禀的共同中心多线性常数满足

\[
C_{c,\mu}^{\boldsymbol a}(\lambda)=\infty.
\]

这里的常数和平均均以同一个 \(\mu\) 定义，零质量方体仍略去。这一命题**只**适用于总质量有限的 \(\mu\)；本轮论证不推出无限总质量权的内禀临界发散。

## 严格证明核对

设 \(M=\mu(\mathbb R^n)\in(0,\infty)\)。把目标
\(\prod_i(1-d_i)^{a_i}\) 连续延拓到紧单纯形
\(\{d_i\ge0,\sum_i d_i=1\}\)。取 \(d_i=1/m\) 得到正值，故

\[
P_{\boldsymbol a}(1)>0.
\]

边界上若有 \(d_i=1\)，乘积为零，因而一个最大点
\(\boldsymbol d^*=(d_i^*)\) 位于原定义域：

\[
0\le d_i^*<1,\qquad \sum_i d_i^*=1,\qquad
P_{\boldsymbol a}(1)=\prod_i(1-d_i^*)^{a_i}>0.
\]

特别地，存在 \(j\) 使 \(d_j^*>0\)。非原子有限测度可把整个
\(\mathbb R^n\) 分割为可测的 \(A_1\dot\cup\cdots\dot\cup A_m\)，使
\(\mu(A_i)=d_i^*M\)；零比例部分可包含任意剩余的零测集。对每个
\(0<t<d_j^*M\)，非原子性再给出 \(F\subset A_j\)、\(\mu(F)=t\)。令

\[
E_i=A_i^c\cup F\quad (1\le i\le m).
\]

全空间分割是必需的：它确保 \(\bigcap_iE_i=F\)，而不是仅在局部盒内成立。
当 \(i\ne j\) 时 \(F\subset A_i^c\)，所以

\[
\frac{\mu(E_i)}M=1-d_i^*\quad(i\ne j),\qquad
\frac{\mu(E_j)}M=1-d_j^*+\frac tM.
\]

从而

\[
G_t:=\prod_i\left(\frac{\mu(E_i)}M\right)^{a_i}
=\left(1-d_j^*+\frac tM\right)^{a_j}
 \prod_{i\ne j}(1-d_i^*)^{a_i}
>P_{\boldsymbol a}(1)\ge\lambda.
\tag{8.1}
\]

严格性不依赖于 \(t\) 的统一余量，只需每个选定 \(t>0\)。因为
\(d_j^*<1\) 且 \(a_j>0\)，增加第 \(j\) 个因子会严格增加乘积。

对任何给定的 \(x\in\mathbb R^n\)，中心方体满足

\[
Q(x,R)\uparrow\mathbb R^n,\qquad
\mu(Q(x,R))\longrightarrow M>0,
\]

且对每一个 \(E_i\)，
\(\mu(E_i\cap Q(x,R))/\mu(Q(x,R))\to\mu(E_i)/M\)。故共同中心乘积收敛到
\(G_t\)。由 (8.1)，某个有限 \(R=R(x,t)\) 已给严格大于 \(\lambda\) 的见证；该方体也自动有正 \(\mu\)-质量。因此光晕为整个 \(\mathbb R^n\)。所取测试集给出商

\[
\frac{\mu(\mathbb R^n)}{\mu(F)}=\frac Mt.
\]

令 \(t\downarrow0\) 即发散。这里半径可以依赖于 \(x\)，正是中心极大算子的定义所允许的；不需要任何统一半径或覆盖选择。

## 边界与量词

- 采用 \(0<\lambda\le P_{\boldsymbol a}(1)\)。阈值通常定义在 \((0,1)\)；若将算子扩展到非正阈值，结论是平凡的，未在本命题中另作声称。
- \(P_{\boldsymbol a}(1)>0\) 与全部 \(d_i^*<1\) 已在紧化论证中明确，故零比例坐标不会产生 \(0^0\) 或零的基底问题。
- \(d_i^*=0\) 时对应 \(A_i\) 可为零测集，\(E_i=A_i^c\cup F\) 仍定义良好；只选择某个 \(d_j^*>0\) 的坐标放置 \(F\)。
- 每一个 \(x\) 都有自己足够大的 \(R(x,t)\)。由于总质量 \(M>0\)，这类方体最终分母为正，符合内禀模型排除零质量方体的约定。
- 仅用到 \(\mu\) 有限且非原子。对本题的 \(w\,dx\)，非原子性自动成立；所附 TeX 片段已直接保留 \(d\mu=w\,dx\) 的假设。若以后要改成一般 Radon 测度，必须显式另加“非原子”。

## 本机精确计算（仅作独立检查）

执行：

~~~text
python3 round08/finite_mass_critical_check.py
~~~

实际输出：

~~~text
N=10: P_a(1)=4/27; global witness product=26/135; excess=2/45; halo quotient=1/f=10
N=100: P_a(1)=4/27; global witness product=103/675; excess=1/225; halo quotient=1/f=100
N=1000: P_a(1)=4/27; global witness product=1003/6750; excess=1/2250; halo quotient=1/f=1000
~~~

该例取 \(m=3\)、\(\boldsymbol a=(1,2,3)\)、
\(\boldsymbol d^*=(2/3,1/3,0)\)，故 \(P_{\boldsymbol a}(1)=4/27\)。在
\([0,1]\) 的 Lebesgue 概率模型中，可取
\(A_1=[0,2/3]\)、\(A_2=(2/3,1]\)、\(A_3=\varnothing\)、
\(F=[0,1/N]\)。充分大的中心区间精确包含整个支撑，程序用有理数验证了
\(G_{1/N}>4/27\) 和商 \(N\to\infty\)。这只是对符号证明的独立实例检验，不替代上述一般证明。

## 可输入 TeX 与审查意见

可直接纳入的章节片段为 [finite_mass_critical.tex](finite_mass_critical.tex)。它适合放在不等指数传递节之后、原先“外层临界发散”说明之前或之后，并应将主文的“临界以下内禀多线性发散仍开放”改为“无限总质量权的情形仍开放”，以免与本命题冲突。

独立复核结论：该候选在上述有限总质量范围内完整成立，适合保留为命题。它不使用错误的 BIP/覆盖选择推理，也不需要把有限点见证错误地当作全体光晕；全空间光晕来自逐点的大尺度极限。
