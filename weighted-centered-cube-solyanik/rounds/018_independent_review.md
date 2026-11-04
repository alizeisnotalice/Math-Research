# 第十八轮独立审查：可比齐次权

## 结论

候选的一阶结论成立。设 \(K\) 为范数单位球、\(\rho=\|\cdot\|_K\)、
\(V=|K|\)，\(0<\tau<n\)，且正 Borel 权满足逐点齐次性及可比性

\[
w(tx)=t^{\tau-n}w(x),\qquad
a\rho(x)^{\tau-n}\le w(x)\le b\rho(x)^{\tau-n}.
\tag{18.1}
\]

令

\[
D=|\{w>1\}|,\quad \kappa=\int_Kw\,dx,\quad
\beta=\tau/n,\quad
S=\frac{n}{\tau}\frac{D^{1-\beta}V^\beta}{\kappa}.
\tag{18.2}
\]

则在第十七轮已经闭合的范数球统一弱型界下，

\[
C^{\rm out}_{K,w}(1-t)-1\sim S\,t^\beta.
\tag{18.3}
\]

若额外 \(w\) 为偶函数且在 \(\mathbb R^n\setminus\{0\}\) 为 \(C^2\)，则对
\[
0<\tau<T_n=\min\left\{1,\frac{n+4-\sqrt{n^2+16}}2\right\}
\]
有
\[
C^{\rm out}_{K,w}(1-t)-1
=S\,t^\beta+S^2t^{2\beta}+o(t^{2\beta}).
\tag{18.4}
\]

对粗糙角向齐次剖面，仅保留 (18.3)，不应声称 (18.4)。

## 齐次分布与精确重排常数

设 \(d=n-\tau>0\)。齐次性给

\[
|\{w>s\}|=D\,s^{-n/d}\qquad(s>0).
\tag{18.5}
\]

这里超水平集的边界具有零 Lebesgue 测度：沿每条范数极射线，方程
\(w(r\theta)=s\) 至多有一个半径。故递减重排积分给任意可测 \(A\)

\[
\mu(A)\le\frac{n}{\tau}D^{1-\beta}|A|^\beta.
\tag{18.6}
\]

取唯一的缩放阈值 \(c\)，使
\[
A_0=\{w>c\},\qquad |A_0|=V.
\]

由 (18.5)，\(A_0\) 有界、包含于某个 \(D_0K\)，且

\[
\mu(A_0)=\frac n\tau D^{1-\beta}V^\beta=\kappa S.
\tag{18.7}
\]

这一点是下界常数恰为 \(S\) 的关键；不能只取任意体积为 \(V\) 的集合。

## 上界

短类环带估计为

\[
O\!\left(\frac ba\,tC^{n-\tau}\right).
\tag{18.8}
\]

因 \(w\in L^1_{\rm loc}\)，平移在固定有界盒上的 \(L^1\) 连续性给

\[
A_C=\inf_{\rho(u)\le C^{-1}}\frac{\mu(u+K)}{\kappa}\longrightarrow1.
\tag{18.9}
\]

对长见证，\(|Q\setminus E|<tVr^n\) 和 (18.6) 给
\[
\mu(Q\setminus E)\le\kappa S r^\tau t^\beta.
\]
先在 \(A_C>St^\beta\) 下由
\[
\mu(E)\ge\kappa r^\tau(A_C-St^\beta)
\]
界定所有实际长见证半径，再取其上确界及共同核心，得到

\[
\frac{\mu(H_l)}{\mu(E)}
\le\frac{S t^\beta}{A_C-St^\beta}.
\tag{18.10}
\]

这与第十七轮的 \(R\) 顺序一致；不需要最大见证球或覆盖选择。
取 \(C=t^{-a}\)、\(0<a<1/n\)，(18.8) 为 \(o(t^\beta)\)，而 (18.10)
给出首项 \(S t^\beta\)。

## 下界与严格性

令 \(q=t^{1/n}\)，取
\[
E=RK\setminus hA_0,\qquad \frac hR=\frac q{1+D_0q}.
\tag{18.11}
\]

当 \(t\) 充分小（例如 \(3D_0h<R\)）时，对
\(\rho(x)<D_0h\) 取 \(r=R-\rho(x)\)。三角不等式给
\[
K(x,r)\subset RK,\qquad hA_0\subset K(x,r).
\]

且
\[
\frac{|hA_0|}{|K(x,r)|}
=\left(\frac hr\right)^n<t.
\]

故孔内点有严格见证。另一方面，Lebesgue 微分定理对范数球适用，所以 \(E\) 的
Lebesgue 密度点（从而 \(\mu\)-几乎处处的 \(E\)）具有极大值 1。于是

\[
C^{\rm out}_{K,w}(1-t)-1
\ge\frac{S t^\beta}{(1+D_0q)^\tau-St^\beta}.
\tag{18.12}
\]

这收敛到 \(S t^\beta\)，证明 (18.3) 的下界。输入有有限正 Lebesgue 体积及
有限正 \(\mu\)-质量，满足外层常数的定义。

## 偶 \(C^2\) 剖面的二阶平移

取偶 cutoff \(\chi\)，其支撑包含于 \(\frac12K\)，并令
\(g=(1-\chi)w\)。当 \(\rho(u)\) 足够小时，\(\chi w\) 的支撑包含于
\(u+K\)，所以
\[
\int_{u+K}\chi w\,dx
\]
不依赖于 \(u\)。函数 \(g\) 在一个包含全部 \(u+K\) 的紧邻域上 \(C^2\)，且为偶；
因 \(K\) 对称，
\[
\int_K\nabla g(y)\,dy=0.
\]
Taylor 公式（Hessian 在该紧邻域有界）遂给
\[
\mu(u+K)=\kappa+O(\rho(u)^2).
\tag{18.13}
\]

由于 \(u=0\) 属于下确界的参数集合，\(A_C\le1\)，故 (18.13) 同时给
\[
A_C=1+O(C^{-2}).
\tag{18.14}
\]

取第十六轮的参数窗口
\[
\frac{\beta}{2}<a<\frac{1-2\beta}{n-\tau}.
\]

则 (18.8) 是 \(o(t^{2\beta})\)，(18.10) 展开为
\(St^\beta+S^2t^{2\beta}+o(t^{2\beta})\)。在 (18.12) 中，多出的
\(t^\beta q\) 为 \(o(t^{2\beta})\) 当且仅当 \(\tau<1\)，故得到 (18.4)。

## 欧氏径向权配中心方体

若 \(n\ge2\)、\(w(x)=|x|_2^{\tau-n}\)，但平均球仍是中心方体 \(K=[-1,1]^n\)，
则

\[
S=\frac n\tau\,
\frac{v_n^{1-\beta}(2^n)^\beta}
{\int_{[-1,1]^n}|x|_2^{\tau-n}\,dx}>1.
\]

严格不等式源于等体积 Euclidean 球对径向递减权的严格重排优于非球方体。这说明
首项 1 是“权的径向几何与平均球匹配”的结论，不能把 Euclidean 径向权无条件
代进中心方体公式。

## 审查意见

一阶、具有适当偶局部正则性的二阶以及欧氏/方体的常数差异均已闭合。主文须保留
超水平最优孔 \(A_0\) 与先界定实际长见证半径再取 \(R\)；二阶正则性已在下节降为
“偶且原点外局部 Lipschitz”，不能再降为仅有有界可测角向剖面。

## 追加审查：二阶条件可降为偶且局部 Lipschitz

该改进成立。可把上节“偶 \(C^2\)”统一替换为：

\[
w\ \text{为偶函数，且在 }\mathbb R^n\setminus\{0\}\text{ 局部 Lipschitz}.
\tag{18.15}
\]

令 \(\rho(u)=\delta<1/2\)，
\[
P=(K+u)\setminus K,\qquad L=K\setminus(K+u).
\]
由于 \(K=-K\)，映射 \(x\mapsto u-x\) 是从 \(P\) 到 \(L\) 的保体积双射：
若 \(x-u\in K\) 而 \(x\notin K\)，则 \(u-x=-(x-u)\in K\)，同时
\((u-x)-u=-x\notin K\)。故

\[
\begin{aligned}
\mu(K+u)-\mu(K)
&=\int_Pw(x)\,dx-\int_Lw(y)\,dy\\
&=\int_P\bigl[w(x)-w(u-x)\bigr]\,dx\\
&=\int_P\bigl[w(x)-w(x-u)\bigr]\,dx,
\end{aligned}
\tag{18.16}
\]

最后一步只用偶性。

又
\[
P\subset(1+\delta)K\setminus K,\qquad |P|=O_K(\delta).
\]
对 \(x\in P\)，以及连接 \(x\) 与 \(x-u\) 的线段，有
\[
\frac12<\rho(\cdot)<\frac32
\]
当 \(\delta\) 充分小时。因此全部相关点处于固定紧环带。局部 Lipschitz 性在该
紧环带上给出统一的小距离 Lipschitz 常数，故
\[
|w(x)-w(x-u)|=O_w(\delta).
\]
代入 (18.16) 得
\[
|\mu(K+u)-\mu(K)|=O_{K,w}(\delta^2).
\tag{18.17}
\]

因为 \(u=0\) 属于 \(A_C\) 的下确界参数域，\(A_C\le1\)；(18.17) 同时给
\[
A_C=1+O(C^{-2}).
\]
所以第十八轮的二阶上界、参数窗口和下界不变。

这一改进可覆盖许多范数不匹配的齐次权（例如角向不光滑但在原点外局部 Lipschitz
的范数幂权）。但**不能**由仅有的有界可测角向剖面推出：那种剖面可能没有
\(O(\delta)\) 的局部平移模量，故仍只应声称一阶结论。

本机 Wolfram Language 15 实际检查一维偶且原点外 Lipschitz 的
\(w(x)=|x|^{-1/2}\)、\(K=[-1,1]\)：在
\(u=1/10,1/20,1/40,1/80\) 时，
\((\mu(K+u)-\mu(K))/u^2\) 约为
\(-0.50157,-0.50039,-0.50010,-0.50002\)，与 (18.17) 的二阶平移损失一致。

## 追加审查：粗糙角向剖面的较小二阶窗口与 \(S=1\) 刻画

这项进一步改进也成立。仅由 (18.1)，在固定范数环带
\(\{1/2\le\rho\le3/2\}\) 上 \(w\) 有统一上界。若 \(\rho(u)=\delta<1/2\)，则
\[
(K+u)\mathbin\triangle K
\subset(1+\delta)K\setminus(1-\delta)K,
\]
其 Lebesgue 体积为 \(O_K(\delta)\)。所以无需偶性或正则性已有
\[
|\mu(K+u)-\mu(K)|=O_{K,a,b}(\delta),
\qquad
A_C\ge1-O(C^{-1}).
\tag{18.18}
\]

令 \(q\) 表示平移损失阶：粗糙情形 \(q=1\)，偶且局部 Lipschitz 情形
\(q=2\)。长类余项为 \(t^\beta C^{-q}\)，小类余项仍为
\(tC^{n-\tau}\)。以 \(C=t^{-a}\) 要求
\[
a>\beta/q,\qquad
a<\frac{1-2\beta}{n-\tau}.
\tag{18.19}
\]
这两个条件可同时满足当且仅当
\[
\tau^2-(n+2q)\tau+qn>0,
\]
即在相关小根区间
\[
0<\tau<
\frac{n+2q-\sqrt{n^2+4q^2}}2.
\tag{18.20}
\]
再与下界展开所需 \(\tau<1\) 相交。故仅有粗糙角向剖面也在
\[
0<\tau<
\min\left\{1,\frac{n+2-\sqrt{n^2+4}}2\right\}
\]
得到 \(S t^\beta+S^2t^{2\beta}+o(t^{2\beta})\)；偶局部 Lipschitz 对应
\(q=2\)，恢复先前较大的 \(T_n\) 窗口。

本机 Wolfram Language 15 已实际核对 (18.19) 与正二次条件的等价性；结果为
True。此粗糙情形的结论来自环带有界性和对称差体积，**不是**把有界角向条件误说成
局部 Lipschitz。

最后，\(S\ge1\)，且
\[
S=1\quad\Longleftrightarrow\quad
w(x)=c\,\rho(x)^{\tau-n}\quad\text{a.e.}
\tag{18.21}
\]
的刻画正确。事实上 \(S\) 是体积 \(V\) 的集合可取得的最大 \(\mu\)-质量与
\(\mu(K)\) 之比。由于 (18.5) 严格递减且其等值面为零 Lebesgue 集，等号迫使
\(K\) 模零等于相应超水平最优集。沿范数极射线使用齐次性，便得到该超水平半径在
几乎每个方向均为常数，亦即 \(w=c\rho^{\tau-n}\) 几乎处处。反向蕴含显然。


---

# 第十八轮统一主稿终审

审查对象为 homogeneous_extension.tex 与 homogeneous_applications.tex。结论：
**数学主线、重新接入顺序和应用常数均通过；仅建议一个严格包含的表述补强。**

1. **定义、依赖与循环。** 齐次权定义 (H.1)--(H.2) 在外层常数之前给出；
   Aldaz/BIP 弱型引理在端点证明前独立完成；多指数回传只在
   weighted-transfer 已输入后使用。新的多指数反向证明明确给出有界
   \(D_1\)、Boolean 原子分割、盒外延拓和 \(E_i=\mathbb R^n\setminus A_i\)，
   因而没有原先的公共交集或循环依赖问题。

2. **分布、\(S=1\) 与可测孔。** 逐点齐次 Borel 代表使超水平集的缩放为集合恒等式；
   比较条件使它们有有限正体积且有界。层析重排给 (H.3)，最优孔
   \(A_0=\{w>\lambda_0\}\) 满足 (H.4)。等号情形也正确：严格径向分布与零
   Lebesgue 等值面迫使 \(K=A_0\) 模零；沿极射线和齐次性即得
   \(w=c\rho^{\tau-n}\) a.e.。孔输入 \(RK\setminus hA_0\) 是可测且合乎
   外层输入条件，Lebesgue 密度点定理保证环部分模 \(\mu_w\)-零在严格光晕中。

   唯一建议：第 117 行应把“\(A_0\subset D_0K\)”写成
   \(A_0\subset\{\rho<D_0\}\)（可任选严格更大的 \(D_0\)）。这样第 121--123
   行对全部孔点的严格见证无需处理可能的边界等号；这只是文字补强，不改变下界。

3. **上界和二阶窗口。** 短类 (H.6) 使用 \(b/a\) 的 Tonelli 比较正确；长类先以
   (H.7) 控制实际见证半径、后取 \(R\) 和共同核心，顺序正确。粗糙角向剖面
   \(q=1\) 与偶局部 Lipschitz 剖面 \(q=2\) 分别给 \(A_C\ge1-O(C^{-q})\)，
   (H.12) 的窗口及额外 \(\tau<1\) 的下界限制均正确。数值检查只能展示，不被用作
   一般证明。

4. **应用常数。** 方体平均/有限 \(p\) 的体积、面锥分解
   \(\kappa=n2^n I_{n,p,\tau}/\tau\) 与 (A.2) 一致；(A.4) 的两个闭式正确。
   一维非对称式 (A.5) 由 \(D=c_+^{1/(1-\tau)}+c_-^{1/(1-\tau)}\) 和
   \(\kappa=(c_++c_-)/\tau\) 直接得到。Wolfram Language 15 实际返回
   \(S_{\tau=1/2,(1,2)}=\sqrt{10}/3\)，与文中例子一致。
