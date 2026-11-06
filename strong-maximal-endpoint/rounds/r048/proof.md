# 第48轮：下端 KKT 带的平方根稳定性与精细比例抵消

日期：2026-11-07。以下是固定独立核的局部结论，不是一般 C_N 的新上界。
统一使用 Φ_N(h)=h(1+log⁺h)^(N−1)，p=(mu−1)/mu，L=logN。

## 1. 有限带平方根引理

设 dν=f(t)dt=R(t)dξ，a=λR，λ>0，输入0≤H≤A，A>1。
Φ_A 是 Φ_N 限制于[0,A]的凸函数，在区间外置+∞。
记 D_a(H)=Φ_N(H)−aH+Φ_A*(a)。
在 I_r=[q,q+r] 假定 a(q)=1、1≤a≤N、
(loga)'≥κ>0、f≤F、|s|≤S。令
\[
 u=(1-H)_+,\quad
 \Delta_r=\lambda^{-1}\int_{I_r}D_a(H)d\xi,\quad
 \kappa_r=\frac{1-e^{-\kappa r}}r.
\]
第47轮的尖点不等式给
\[
 \Delta_r\ge\int_{I_r}(1-a^{-1})u\,d\nu.
\]
写 x=t−q，则1−a(q+x)^−1≥1−e^(−κx)≥κ_r x。
若 M=∫_I u dν，因0≤uf≤F，层蛋糕积分给
\[
 \int_0^r x u(q+x)f(q+x)dx
 \ge\int_0^{M/F}(M-Fy)dy=M^2/(2F).
\]
因此对所有这些 H，严格有
\[
 \boxed{\left|\int_{I_r}u s\,d\nu\right|^2
 \le\frac{2S^2F}{\kappa_r}\Delta_r.}\tag{1}
\]
这是不含端点距离发散因子的平方根控制，而非第47轮已经排除的统一线性控制。

若 a≤a_max<N，则超出1的部分另由
D_a(H)≥(N−a)(H−1)_+ 控制，得到
\[
 \left|\int_{I_r}(H-1)s\,d\nu\right|
 \le S\sqrt{2F\Delta_r/\kappa_r}
 +S\frac{a_{\max}}{N-a_{\max}}\Delta_r.\tag{2}
\]
这只包含 I_r 上的贡献；积分号不能直接扩为全空间。

## 2. 局部常数是尖锐的，两个极限有先后

先固定 N、q，并假定 f、s、k=(logR)' 在q连续，f(q)>0、k(q)>0。
取完整 capped-KKT 参考输入，在[q,q+ell]将其高度1平台降为0，其他位置不变。
在a=1的单点取KKT允许值0；低于该层本已为0，带后跳到1。
因此这是非负有界输入，缺口仅出现在该带。若a在全轴非减（本项目标量块核满足），
该修改也保持非减；一般有限带引理本身无需全局单调。
当 ell↓0 时，
\[
 \int_q^{q+\ell}u s\,d\nu=s(q)f(q)\ell+o(\ell),\qquad
 \Delta_\ell=\tfrac12 k(q)f(q)\ell^2+o(\ell^2).
\]
于是
\[
 \boxed{\lim_{\ell\downarrow0}
 \frac{|\int u s\,d\nu|^2}{\Delta_\ell}
 =\frac{2f(q)s(q)^2}{k(q)}.}\tag{3}
\]
该式先取固定维数的局部极限。以下渐近式再对其系数取 N→∞；
没有未经证明地交换两个极限。

## 3. 整个临界阈值带上的精确 score 展开

Z_c=Exp(1)+c Bern((1+e^c)^−1)，f_(N,c) 为N重卷积。
R_c(t)=e^t f_(N,c)(t)，s_(N,c)=∂c logf_(N,c)。
参数偏导固定空间自变量；随后仅在取值时令t=q。
令 barPhi(z_N)=sqrt(mu/N)z_N，d_N=sqrt(mu N)z_N，
取
\[
 c=\mu+h,\quad |h|\le K\sqrt{L/N},\quad
 q=N\mu-d_N+O(L).
\]
所有O常数有固定上界。第17、19轮统一平台根及根附近平坦密度估计给
\[
 f_{N,c}(q)\sim L/N.\tag{4}
\]
这不要求q恰好等于平台根：平台根本身在该临界带，两个阈值只差O(L)。

取负鞍点ell>0使 Λ_c'(−ell)=q/N，设
ell_N=z_N/sqrt(mu N)，d=c(1+ell)。均值/方差Taylor给
\[
 \ell=\ell_N+O(L/N),\qquad d=\mu+h+\mu\ell_N+O(L/N).
\]
第45轮精确桥公式及其紧参数展开为
\[
 s_{N,c}(q)=Np_c\left[1-\frac{u(d)}{u(c)}B_N(d)\right],\quad
 u(x)=\frac{x}{1+e^x},
\]
\[
 B_N(d)=1+\frac{1/2+O(\sqrt{L/N})}{N}+O(N^{-2}).
\]
用u'(mu)=0、(logu)''(mu)=−1/mu，且
N(L/N)^(3/2)→0，展开得到一致的**绝对误差o(1)**公式
\[
 \boxed{s_{N,c}(q)=p\left\{
 \frac{h z_N\sqrt N}{\sqrt\mu}+\frac{z_N^2}2-\frac12\right\}+o(1).}\tag{5}
\]
例如logu(d)−logu(c)=−h ell_N−mu ell_N²/2+O((L/N)^(3/2))，
Np_c的替换误差也是O(L^(3/2)/sqrtN)=o(1)，保留桥项−p/2。
同一倾斜密度在中心的g'/g=O(1/N)，给
\[
 k(q)=1+f'_{N,c}(q)/f_{N,c}(q)=1+\ell+O(N^{-1})\to1.\tag{6}
\]
这是空间logR导数，不是平台平均对q的导数。

## 4. 精细比例使下端常数降低两个logN因子

定义
\[
 \bar c_N=\mu-\frac{\sqrt\mu\,z_N}{2\sqrt N}.
\]
代入(5)，两个z_N²项精确抵消，得到
\[
 \boxed{s_{N,\bar c_N}(q)=-p/2+o(1).}\tag{7}
\]
也可直接看c、d关于mu的二次对称性：
c=mu−mu ell_N/2，d=mu+mu ell_N/2+O(L/N)，
N[u(d)−u(c)]=O(L^(3/2)/sqrtN)=o(1)，只留下桥项。
该结论对任意固定宽度[q,q+r]一致成立。

由(3)—(7)，三种比例的尖锐局部系数分别为
\[
 \frac{2f(q)s(q)^2}{k(q)}\sim
 \begin{cases}
 (p^2/2)L^3/N,&c=\mu,\\
 (p^2/2)L(\log L)^2/N,&c=\mu-\sqrt\mu\sqrt L/(2\sqrt N),\\
 (p^2/2)L/N,&c=\bar c_N.
 \end{cases}\tag{8}
\]
第二行用Mills根方程：z_N²=L−2logL−log(2pi mu)+o(1)，
故粗sqrtL参数留下score=−(p/2)logL+O(1)。
这是精细比例对局部稳定性常数的改善，不是C_N下界系数的改善。

更直接地，在精细比例和任意固定r>0下，(2)给
\[
 \left|\int_{I_r}(H-1)s\,d\nu\right|
 \le\left(\frac p2+o(1)\right)
 \sqrt{\frac{2r}{1-e^{-r}}\frac LN\Delta_r}
 +O(\Delta_r/N).\tag{9}
\]
这里F~L/N、κ→1、S=p/2+o(1)，a_max=e^(r+o(1))，
所有估计均来自同一固定宽度带。H可任意非负且有界，无需输入单调。

## 范围与来源

(9)没有处理a<1一侧、a约N的上端、soft区或多块角向coarea。
精细比例也没有使score在全空间一致有界；第45轮全状态下确界在该窗口仍可为−O(L)。
因此本轮不把局部平方根估计当作全局近极值稳定性，不改变C_N账本。

来源：[第17轮](../r017/report.md)、[第19轮](../r019/report.md)的临界根与密度，
[第45轮](../r045/proof.md)的精确缩放桥及紧参数展开，
[第47轮](../r047/deficit.md)的完整capped-KKT预算缺口。
