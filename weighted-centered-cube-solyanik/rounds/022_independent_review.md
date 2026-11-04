# 第 22 轮独立审查：角向稳定不等式

## 结论

三个候选不等式及给出的普适尖锐性机制成立。唯一需补的技术措辞是 (i) 中
若 \(f\) 未被假设远离零，二阶微分应先在 \(0\le t<1\) 进行，再令
\(t\uparrow1\)；不能无说明地在端点直接使用含 \(u_t^{p-2}\) 的式子。

令 \(g=f-1\)，\(\mathbb Ef=1\)，
\[
S=\|f\|_p,\qquad p=\frac n{n-\tau}>1.
\]
这里 \(S\) 确为 (H.2) 中的端点系数：按 cone probability 写
\(a=(\mathbb Ea)f\)，则第 21 轮的极坐标计算给
\(S(h)=\|a\|_p/\mathbb Ea=\|f\|_p\)。

## \(1<p\le2\)

对 \(u_t=1+tg\)，在 \(0\le t<1\) 定义
\[
A=\mathbb Eu_t^p,\quad B=\mathbb Eu_t^{p-1}g,\quad
C=\mathbb Eu_t^{p-2}g^2,\quad F(t)=A^{2/p}.
\]
直接微分得
\[
F''=2(2-p)A^{2/p-2}B^2+
2(p-1)A^{2/p-1}C.
\]
前项非负。Hölder 给
\[
\mathbb E|g|^p\le C^{p/2}A^{(2-p)/2},
\quad\text{即}\quad
A^{2/p-1}C\ge\|g\|_p^2.
\]
故 \(F''\ge2(p-1)\|g\|_p^2\)。又 \(F(0)=1,F'(0)=2\mathbb Eg=0\)。
先积分至任意 \(r<1\)，再以有界 \(f\) 的 \(L^p\) 连续性令 \(r\uparrow1\)，得
\[
\|f-1\|_p^2\le\frac{S^2-1}{p-1}.\tag{R22.1}
\]
取任意 cone-balanced 半集上 \(f=1+\delta\)、补集上 \(f=1-\delta\)，
令 \(\delta\downarrow0\)，有
\(S^2-1=(p-1)\delta^2+o(\delta^2)\)。故常数 \(1/(p-1)\) 普适尖锐。

## \(p\ge2\)

标量不等式
\[
x^p-1-p(x-1)\ge|x-1|^p\quad(x\ge0)
\]
成立。对 \(x\ge1\) 可由导数和
\((1+t)^{p-1}\ge1+t^{p-1}\) 验证；对 \(x=1-s\in[0,1]\) 则导数为
\(p[1-(1-s)^{p-1}-s^{p-1}]\ge0\)。取期望并用 \(\mathbb Ef=1\)，得到
\[
\|f-1\|_p^p\le S^p-1.\tag{R22.2}
\]

此外 Bregman 余量有
\[
x^p-1-p(x-1)\ge(p-1)(x-1)^2.\tag{R22.3}
\]
Taylor 积分式为
\[
p(p-1)(x-1)^2\int_0^1(1-s)[1+s(x-1)]^{p-2}\,ds.
\]
括号底数至少为 \(1-s\)，所以积分至少
\(\int_0^1(1-s)^{p-1}ds=1/p\)。取期望后得
\[
\|f-1\|_2^2\le\frac{S^p-1}{p-1}.\tag{R22.4}
\]

对非原子 cone probability，令
\(\mathbb P(A)=\varepsilon=M^{-(p+1)}\)，令 \(f=M\) 于 \(A\)，
\(f=(1-\varepsilon M)/(1-\varepsilon)\) 于补集。令 \(M\to\infty\)，
(R22.2) 左侧和右侧之比趋于一。再令 \(f=q\) 于概率为 \(q\) 的集合，
并令 \(f=1+q\) 于补集，则
\[
\mathbb E(f-1)^2=q(1-q),\qquad
\mathbb Ef^p-1=(p-1)q+O(q^2),
\]
故 (R22.4) 的常数也尖锐。此两个尖锐性构造需要可取任意小质量集合；
因此可安全宣称为 \(n\ge2\) cone measure 的尖锐性，而不应声称固定一维
两原子 cone measure 同样尖锐。\(p=2\) 时两个不等式本身为恒等式。

## 空间 \(L^1\) 应用

令 \(\bar a=\mathbb Ea\)。相对空间误差满足
\[
\frac{\int_K|h(x)-\bar a\rho(x)^{\tau-n}|\,dx}{\int_Kh(x)\,dx}
=\mathbb E|f-1|.
\]
由 (R22.1) 或 (R22.4) 并 Cauchy--Schwarz/Hölder，
\[
\mathbb E|f-1|=O_p(\sqrt{S-1})\qquad(S\downarrow1).
\]
balanced \(1\pm\delta\) 例有 \(\mathbb E|f-1|=\delta\) 而
\(S-1\asymp\delta^2\)，故指数 \(1/2\) 不能提高。这里不需要、也不应声称
最佳 \(L^1\) 常数。

本机 Python 已对方体 cone 例的 \(\mathbb E\theta_1^2\) 及上述端点代数作
有理算术复核；数值仅用于算术检查，不替代证明。


## 附录：相对 \(L^1\) 的尖锐渐近常数

该强化结论成立。写
\[
D=S^p-1=\mathbb E\phi_p(f),\qquad
\phi_p(x)=x^p-1-p(x-1),\qquad \mathbb Ef=1.
\]
对任意固定 \(0<\delta<1\)，在 \(|x-1|\le\delta\) 上，积分 Taylor 式给
\[
\phi_p(x)\ge\frac{p(p-1)}2m_\delta(x-1)^2,\qquad
m_\delta=\min\{(1-\delta)^{p-2},(1+\delta)^{p-2}\}.
\]
而在闭尾部 \(\{x\ge0:|x-1|\ge\delta\}\)，
\[
c_{p,\delta}:=\inf\frac{\phi_p(x)}{|x-1|}>0.
\]
确实，商在有限区间上连续且严格正；当 \(x\to\infty\) 时该商
\(\sim x^{p-1}\to\infty\)，所以不会在无穷远丢失正下界。

分割期望、Cauchy--Schwarz 和尾部线性界给统一估计
\[
\mathbb E|f-1|
\le\sqrt{\frac{2D}{p(p-1)m_\delta}}+\frac{D}{c_{p,\delta}}.\tag{R22.5}
\]
先令 \(D\downarrow0\)，再令 \(\delta\downarrow0\)，并用
\(D/(S-1)\to p\)，得到
\[
\limsup_{S\downarrow1}\frac{\mathbb E|f-1|}{\sqrt{S-1}}
\le\sqrt{\frac2{p-1}}.\tag{R22.6}
\]
(R22.5) 的常数只依赖 \(p,\delta\)，故同一论证正是对所有 (H.1)
角向剖面的统一上界：
\[
\lim_{\eta\downarrow0}\sup_{1<S(h)\le1+\eta}
\frac{d(h)}{\sqrt{S(h)-1}}
\le\sqrt{\frac2{p-1}},
\quad
d(h)=\frac{\int_K|h-(\mathbb Ea)\rho^{\tau-n}|}{\int_Kh}.
\]
前一轮的 cone 坐标式给 \(d(h)=\mathbb E|f-1|\)。

下界取 cone measure 的互补等质量两半，并令 \(f=1+\varepsilon\)、
\(1-\varepsilon\)。中心对称的任意范数球都可取这种两半：
\(n\ge2\) 时可用 \(\theta_1\) 的正负号（零超平面的锥测度为零），
\(n=1\) 时就是两个端点。于是
\[
d=\varepsilon,\qquad
S-1=\frac{p-1}{2}\varepsilon^2+o(\varepsilon^2),
\]
给出相反不等式并证明
\[
\boxed{\displaystyle
\lim_{\eta\downarrow0}\sup_{1<S(h)\le1+\eta}
\frac{d(h)}{\sqrt{S(h)-1}}=\sqrt{\frac2{p-1}}.}\tag{R22.7}
\]
故尖锐性对所有维数成立，不需要 cone measure 非原子性。

这在 \(p>2\) 时严格改善上一轮由 \(L^2\) 界所得的渐近常数
\(\sqrt{p/(p-1)}\)；在 \(1<p\le2\) 时与 (R22.1) 加 Hölder 给出的
一阶渐近常数相同，但现在确定了全类的最优性。没有必要声称有限偏差范围内的
最佳 \(L^1\) 常数。

本机 Python 对 \(p=1.25,1.75,2,3,5\) 的 balanced 两值剖面复核了
比值向 \(\sqrt{2/(p-1)}\) 收敛；对 \(p=3,5\) 亦复核该值严格小于
旧 \(L^2\) 渐近常数。

主代理裁决：已落实t<1微分后趋近端点、固定K/n/τ的统一上确界、p=2恒等与p>2尖锐实现范围；将含糊的“逐点稳定性”精确改为“不推出L∞控制”，不以本质上确界发散声称否定a.e.点态收敛。经典Lp来源与本文角向应用分别说明。
