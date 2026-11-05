# 真实第三反馈与三维微孔：可实现族之间的九阶改进

2026-11-05。固定维数。只写本独立目录，旧论文、冻结源及他人文件只读。本稿给出实际构造的比较定理，不将已有 \(O(av^9)\) 余项拟合成一个绝对九阶系数。

## 1. 可检验目标和结论范围

记
\[
 u=t^{1/n},\quad v=2u,\quad a=v^n,\quad
 g=(1+2u)^{-1},\quad B_n=\binom n2,\quad K_n=B_nn^2/27.
\]
沿用已证宏观斜格孔的参数，但所有对象在下文重新定义。取
\(\kappa(u)\) 在固定正紧区间内，最终令 \(\kappa(u)\to8n/3\)。
外半宽写为 \(w=ug+\Lambda(u)u^6\)，其中 \(\Lambda\) 有界。

**实际比较命题 A。** 若原参数对两次反馈层 \(U_2\) 安全，即
\((2w)^n-a(1-2w)^n>|U_2|\)，则对每个固定 \(0<\theta<1\)，可以改用明确的四步层 \(U_4\)，并将参数改为
\[
 \Lambda^-=\Lambda-
 \theta\,\frac{c_n\kappa^6}{1024n^5g^4}\,u^3,\qquad
 c_n=2^{-19-12(n-2)}.
 \tag{A}
\]
新参数仍安全，原安全余量不需要重新购买。原、新实际有限输入的完整 halo 比值分别记为 \(\mathscr R_2(u),\mathscr R_4(u)\)，则可使
\[
 \mathscr R_4(u)-\mathscr R_2(u)
 =\left(\frac{\theta c_n n^2}{729}+o_n(1)\right)av^9>0.
 \tag{B}
\]
这里的保守常数来自一个真正中央角块，不假定一个平面常数能自动迁移到高维全部源层。

**实际比较命题 B。** 对任一这样的安全 \(U_4\) 宏观族，可以在 \(U_4\) 内真实造微孔，再用一次显式保护操作维持所有宏观见证。适当减小 \(\Lambda\) 后，实际比值增加
\[
 (2K_n+o_n(1))av^{n+6}.
 \tag{C}
\]
特别 \(n=3\) 时是 \(2av^9+o(av^9)\)，可与 (B) 顺序合并。\(n\ge4\) 时该微孔项晚于九阶。

另证明真实平面 min 帽迭代的第三反馈常数 \(241/60\)。二维纯宏观族可据此把 (B) 的保守增益提升为任意固定 \(\theta<1\) 倍的
\[
 \frac{482}{10935}av^9+o(av^9).
 \tag{D}
\]
这尚不计算二维已经加入微孔的最佳族的九阶项。微孔保护在二维恰可能影响该阶，旧安全常数也在该阶，不能省略它们。

以上均是明确实际族的差值结论。原 \(U_2\) 族只展开到八阶且其九阶安全项尚未识别，因此本稿不据 (B) 或 (C) 宣布 \(C_n-P_8\) 的绝对九阶系数。

## 2. 固定实际质量的坏集算子

令
\[
 e=\kappa u^3,\quad T=(1-e)I+eJ,\quad d=\det T,\quad
 C_w=[-w,w]^n,\quad \rho=1-2w,\quad M=a\rho^n .
\]
\(J\) 为全一矩阵。每个源点取同一个规则半径
\[
 s=w-\|x\|_\infty,\qquad r(x)=\rho+s.
\]
自身总孔质量始终是 \(M\)，不按测试删除集重新计算。置
\[
 F_w(s)=\frac{a}{(2w)^{n-2}}\big[(\rho+s)^n-\rho^n\big],
 \quad f=F_w'(0),\quad \delta=\frac e f,\quad
 A_w=(2w)^{n-2}e^2 .
\]
均匀地
\[
 f=4nu^2g[1+O_n(u^5)],\quad \delta\asymp_n u,\quad
 e/w=O_n(u^2).
 \tag{1}
\]
对 Borel \(V\subset C_w\) 定义
\[
 R_{w,V}(x)=\frac1{(2w)^{n-2}}
   \sum_{k\in\mathbb Z^n\setminus\{0\}}
       |(Tk+(C_w\setminus V))\cap Q(x,\rho+s)|,
\]
\[
 \mathcal B_w(V)=\{s=0\}\cup\{R_{w,V}(x)\ge F_w(s)\}.
 \tag{2}
\]
这是原阈值的真实坏集，不是只扣一个名义帽的安全测试。

有限邻孔清单直接由几何得到。因为
\(r+w=1-\|x\|_\infty\)、
\(\|Tk\|_\infty\ge[1-(n-1)e]\|k\|_\infty\)，指标模至少二者不能进入。
对 \(k\in\{-1,0,1\}^n\)，若 \(\sum k_i\ne0\)，取与和同号的非零坐标，其晶格坐标模至少为一，至多边界接触。只有平衡非零指标进入，支撑 \(j\ge2\) 为偶数；此时 \(Tk=(1-e)k\)。支撑轴交长均不超过 \(e\)，其余方向整宽进入；自身全孔由 \(r\ge\rho>2w\) 捕获。

所有坏点除零测外边界外满足
\[
 s\le h_H:=F_w^{-1}(3^ne^2)\le C_ne\delta.
 \tag{3}
\]
在最近面符号和全局最浅轴 \(q\) 的坐标中，写
\(d_i=s+y_i,\ y_q=0\)。Jacobian 为一。活动平衡目标的各支撑穿透为 \((e-y_i)_+\)，不依赖 \(s\)。故 (2) 在每张图上是
\[
 0\le s\le
 \min\{w-\max y_i,\ F_w^{-1}(R_{w,V}(y))\}.
 \tag{4}
\]
所有图都保留，没有删去三近面。有限帽交量连续，所以各 \(\mathcal B_w(V)\) 闭、图高度连续；边界由有限图分界和竖纤维高度图组成，零测。

反单调性是精确集合关系。定义
\[
 U_0=\mathcal B_w(\varnothing),\quad U_{j+1}=\mathcal B_w(U_j)
 \quad(0\le j\le3).
\]
则
\[
 U_1\subset U_3\subset U_4\subset U_2\subset U_0,\qquad
 \mathcal B_w(U_4)\subset U_4.
 \tag{5}
\]
因此只要可删量超过 \(|U_4|\)，在 \(C_w\setminus U_4\) 切取质量 \(M\) 的孔即可得到所有孔点的原严格见证。

## 3. 整个真实算子的小常数

对任意 \(V,W\subset C_w\)，
\[
 |\mathcal B_w(V)\triangle\mathcal B_w(W)|
 \le L_u|V\triangle W|,\qquad L_u\le C_n\delta=O_n(u).
 \tag{6}
\]
证明使用 (4)、\(F_w'\ge f\) 和真实顶端截断的 1-Lipschitz 性。打开 Tonelli 固定一个目标点。若目标支撑为 \(j\)，最浅轴在支撑内，可出现的源底面至多
\(C_ne^{j-1}w^{n-j}\)；否则至多 \(C_ne^jw^{n-j-1}\)。
除以 \(f(2w)^{n-2}\)，分别为
\[
 C_n\delta(e/w)^{j-2},\qquad
 C_n\delta(e/w)^{j-1}.
\]
有限平衡清单求和得到 (6)，任意薄帽和跨角输入均在其中。
同样积分完整名义帽，
\[
 |U_0|\le C_n A_w\delta.
 \tag{7}
\]
于是
\[
 |U_2\setminus U_1|\le L_u^2|U_0|
       \le C_nA_w\delta^3,\qquad
 |U_2\setminus U_4|\le L_u^3|U_0|
       \le C_nA_w\delta^4 .
 \tag{8}
\]
这些界只比较真实集合，不是形式 Neumann 展开。

## 4. 单个中央短角释放的实际正体积

本节给可直接审核的保守下界
\[
 |U_2\setminus U_4|\ge c_n A_w\delta^4,
 \qquad c_n=2^{-19}8^{-4(n-2)}.
 \tag{9}
\]
固定一个有向短角，例如平衡指标 \((1,-1,0,\ldots,0)\)，并固定其他轴的最近面符号。记两个支撑面深为 \(X,Y\)，先取排序 \(X=s,Y=s+eD\)，其中
\[
 D\in[1/8,1/4],\qquad q=1-D\in[3/4,7/8].
 \tag{10}
\]
其他底变量取 \(y_h\in[2e,w/3]\)。小 \(u\) 后每个区间长至少 \(w/4\)，物理面深小于 \(w/2\)，所以不触顶端截断。这些点的同一个原见证只有该短角目标产生正交量，所有其他平衡支撑均因额外坐标深度大于 \(e\) 而不进入。

记此目标帽为
\(\Pi_q=[0,e]\times[0,eq]\times[-w,w]^{n-2}\)，
所有坐标均从目标相应面对端解释。由于全反射和轴置换对称，目标的反向短角与源的短角具有相同坏层。以下每个正质量下界因此可反馈给原角。

先在目标中取
\[
 0\le X\le e\delta/8,\qquad e/8\le Y\le e/4,
\]
其余面深取 \([2e,w/3]\)。这里 \(F_w(X)\le2fX\le e^2/4\)，而完整短角交量至少 \(3e^2/4\)，故整块在 \(U_0\) 中；对所有 \(q\ge1/2\) 都在 \(\Pi_q\) 内。因此
\[
 |U_0\cap\Pi_q|\ge b_0 A_w\delta,\qquad
 b_0=8^{-(n-2)}/64.
 \tag{11}
\]
高度窗 (3) 上 \(F_w'\le2f\)。在 (10) 的中央源图内，以 \(h_j\) 记 \(U_j\) 的真实物理高度，则
\[
 h_0-h_1\ge
 \frac{|U_0\cap\Pi_q|}{2f(2w)^{n-2}}
 \ge (b_0/2)e\delta^2 .
 \tag{12}
\]
此差层的两个支撑面深都小于 \(e/2\)，其余横向区间至少长 \(w/4\)。于是对所有 \(q\ge1/2\)，
\[
 |(U_0\setminus U_1)\cap\Pi_q|
 \ge b_1 A_w\delta^2,\qquad b_1=2^{-10}8^{-2(n-2)}.
 \tag{13}
\]
这里 \(1/8\) 来自 (10) 的底变量长度，\(1/2\) 来自逆导数，另一个 \(8^{-(n-2)}\) 是真实横向比例。

由 (8)，取小 \(u\) 后可减去整个 \(U_2\setminus U_1\)，从而
\[
 |(U_0\setminus U_2)\cap\Pi_q|
 \ge (b_1/2)A_w\delta^2\quad(q\ge1/2).
 \tag{14}
\]
现在 \(U_3\setminus U_1=\mathcal B_w(U_2)\setminus\mathcal B_w(U_0)\)，同一中央图的高度差至少 \((b_1/4)e\delta^3\)。再作 (10) 的积分，得到
\[
 |(U_3\setminus U_1)\cap\Pi_q|
 \ge b_2A_w\delta^3,\qquad b_2=2^{-15}8^{-3(n-2)}.
 \tag{15}
\]
最后 \(U_2\setminus U_4=\mathcal B_w(U_1)\setminus\mathcal B_w(U_3)\)，再次用 (15)、逆导数和同一底面，给
\[
 |U_2\setminus U_4|\ge (b_2/16)8^{-(n-2)}A_w\delta^4,
\]
恰为 (9)。全程使用实际层差、同一帽和同一原半径；没有把三个不交假设或正名义体积当成真实正差。

## 5. 改变外半宽的体积稳定性

后面会改变 \(\Lambda\)。为避免用旧孔宽直接代新孔，先证明对有界 \(\Lambda,\Lambda'\) 和 \(0\le j\le4\)，
\[
 \big||U_j(\Lambda')|-|U_j(\Lambda)|\big|
 \le C_na u^{10}|\Lambda'-\Lambda|.
 \tag{16}
\]
这里 \(\kappa,u\) 固定，所有 \(U_j\) 均按自己的新 \(w,\rho,F_w\) 重新构造。

用内禀最近面坐标 \((\sigma,d_1,\ldots,d_n)\)，其中
\(x_h=\sigma_h(w-d_h),\ 0\le d_h\le w\)。
将不同 \(w\) 的集合嵌入同一个有限符号直和的 \(\mathbb R_+^n\)，在各自域外补零。这个比较不把物理边界移动错当作体积变化。

在这些坐标中，所有支撑前缀区间仍是 \([0,(e-y_h)_+]\)，不随 \(w\) 改变；非支撑目标坐标仅在 \([0,w]\) 内积分并求最近面符号和。变化的只有这些非临界积分限、归一化 \((2w)^{n-2}\)、\(F_w\) 及源中分面上界。
固定高度窗上
\[
 |\partial_wF_w(s)|\le C_nF_w(s)/w,\qquad F_w'(s)\asymp f.
\]
所以在固定目标测试数据下，高度变化由名义高度乘
\(C_n|w'-w|/w\) 上包。
目标非临界积分限的差条同样只占完整名义帽的这个比例。
源域或顶端截断的变化用层积分处理：在 \(s\le h_H\) 时，活动支撑坐标都远离中分面，只有自由底变量的区间端随 \(w\) 移动；Fubini 给相同的名义体积比例界。整个名义体积为 \(O_n(au^5)\)。

因此对两组内禀测试数据 \(V,V'\)，真实算子满足
\[
 \|\mathcal I_w\mathcal B_w(V)-\mathcal I_{w'}\mathcal B_{w'}(V')\|_1
 \le C_nu\|\mathcal I_wV-\mathcal I_{w'}V'\|_1
       +C_nau^5|w'-w|/w.
 \tag{17}
\]
第一个常数由 (6) 的同一核计数取得；第二项的边界差条用完整帽上包，故不要求测试缺集有规则边界。对从空集开始的固定四次迭代应用 (17)，再用
\(w'-w=u^6(\Lambda'-\Lambda)\)、\(w\asymp u\)，得到 (16)。

## 6. 比较命题 A 的实际有限化

令 \(\Lambda^-=\Lambda-\theta c_n\kappa^6u^3/(1024n^5g^4)\)，
所有新对象加上上标 \(-\)。精确可删量
\[
 D_\Lambda=(2w)^n-a\rho^n,\qquad
 D_\Lambda/a=n\Lambda u^5g^{n-2}+O_n(u^{10})
 \tag{18}
\]
在有界参数区间可求导；其关于 \(\Lambda\) 的导数为
\(n u^5g^{n-2}+O_n(u^{10})\)。
又由 (1)，
\[
 \frac{A_w\delta^4}{a}
 =\frac{\kappa^6g^{n-6}}{1024n^4}u^8[1+O_n(u^5)].
 \tag{19}
\]
所以 \(D_{\Lambda^-}-D_\Lambda=-\theta c_nA_w\delta^4+o(au^8)\)。
另一方面 (16) 给
\(|\,|U_2^-|-|U_2|\,|=O_n(au^{13})\)，(9) 对新参数仍有效。
因此
\[
 D_{\Lambda^-}-|U_4^-|
 \ge D_\Lambda-|U_2|+(1-\theta)c_nA_w\delta^4+o(au^8)>0.
 \tag{20}
\]
原安全余量未被消耗；严格正的新增空间来自第三反馈。

这个输入类非空，且无需预先知道其九阶系数。例如既有宏孔证明供应
\[
 \Lambda=
 \frac{B_n\kappa^3}{8n^2g}
 -\frac{5B_n\kappa^4}{96n^3g^2}u+
 \left[\frac{13B_n\kappa^5}{512n^4g^3}
       +\frac{B_n(n-2)\kappa^4}{32n^2g^2}\right]u^2+cu^3,
\]
其中先选充分大的固定 c 再取小 u 即对 \(U_2\) 安全。本稿也适用于比此安全选择更紧的任何有界安全参数；(20) 保留其原有余量，因而新增空间不是仅由任意减小一个粗安全常数所得。

在 \(C_{w^-}\setminus U_4^-\) 中连续分位切取质量
\(M^-=a(\rho^-)^n\) 的紧宏孔。由 (5)，所有孔点的规则半径严格满足原缺量条件。原、新周期的 determinant 同为 \(d\)，而
\[
 M^--M=
 2n a g^{n-1}
 \frac{\theta c_n\kappa^6}{1024n^5g^4}u^9+o(au^9).
 \tag{21}
\]
比值 \(d/(d-M)\) 的质量导数为 \(1+O_n(a+u^6)\)，因为 \(d=1+O_n(u^6)\)。取
\(N=\lceil u^{-12}\rceil\)，外盒 \(\Omega_N=T([-N,N]^n)\)，保留全部边界切孔。
周期输入质量精确为 \((2N)^n(d-M)\)。内层全部孔点及其原见证保留，边界孔质量比误差 \(O_n(a/N)=o(av^9)\)；凸外盒之外所有半径由分离半空间的半体积界排除，输入密度点进入 halo。因此这是完整极大算子比值而非指定半径子族。
将 (21) 除 \(av^9\)，再令 \(\kappa\to8n/3\)，即得 (B)。

## 7. 真实平面 min 帽的第三常数

这一节的尖常数先只对真实平面迭代证明，不能仅凭量纲把它作为全部高维层差的系数。

在单位短角中，用 \(d=|X-Y|\) 记两排序，设薄层高度为 \(H(d)\)，满足 \(0\le H(d)\le\delta(1-d)\)。真实帽面积为
\[
 \mathcal M(H)(q)=
 \int_0^q\min\{H(d),q-d\}\,dd+
 \int_0^1\min\{H(d),q\}\,dd .
 \tag{22}
\]
第一项是 \(X\le Y\) 的排序，第二项是 \(Y\le X\) 的排序；\(H(d)\le1-d\) 已保证单位外帽另一个轴的裁剪。没有丢掉 \(q\) 小于薄高的帽。

令
\[
 H_0(d)=\delta(1-d),\qquad
 H_{j+1}(d)=\delta[(1-d)-\mathcal M(H_j)(1-d)].
 \tag{23}
\]
因为真实帽面积不超过 \(q\)，高度始终非负，且不超过 \(\delta(1-d)\)。这恰为两方向相等的原平面反馈。

定义正线性算子
\[
 (\mathcal A\phi)(y)=\int_0^{1-y}\phi(z)\,dz+\int_0^1\phi(z)\,dz,
 \qquad \phi_0(y)=1-y.
\]
以下是处理任意近零帽的统一估计。若 \(0\le G\le H\le\delta\)，\(W=H-G\)，则
\[
 0\le
 \left(\int_0^qW+\int_0^1W\right)
       -[\mathcal M(H)(q)-\mathcal M(G)(q)]
 \le\int_{(q-\delta)_+}^qW(d)\,dd+
           \mathbf1_{\{q<\delta\}}\int_0^1W.
 \tag{24}
\]
因为当截断阈值至少为 \(\delta\) 时 min 的差恰为 \(W\)，其余处差仍在 \([0,W]\)。对 q 积分，两项 Fubini 出现长度均不超过 \(\delta\)，所以 (24) 的 \(L^1\) 误差至多 \(2\delta\|W\|_1\)。这个常数不依帽的最小宽度。

由 (24)，正高度差每反馈一次为
\[
 \delta\mathcal A W+O_{L^1}(\delta^2\|W\|_1).
 \tag{25}
\]
于是
\[
 H_0-H_1=\delta^2\mathcal A\phi_0+O_{L^1}(\delta^3),\qquad
 H_0-H_2=\delta^2\mathcal A\phi_0+O_{L^1}(\delta^3),
\]
\[
 H_3-H_1=\delta^3\mathcal A^2\phi_0+O_{L^1}(\delta^4),\qquad
 H_2-H_4=\delta^4\mathcal A^3\phi_0+O_{L^1}(\delta^5).
 \tag{26}
\]
这些式只对真实正层差应用小算子，绝非将集合直接作几何级数相减。

本机精确多项式核给出
\[
 \mathcal A\phi_0=1-y^2/2,\quad
 \mathcal A^2\phi_0=5/3-y/2-y^2/2+y^3/6,
\]
\[
 \mathcal A^3\phi_0=31/12-5y/6-y^2/2+y^4/24,\quad
 2\int_0^1\mathcal A^3\phi_0=\frac{241}{60}.
 \tag{27}
\]
前两次积分分别为 \(5/3,31/12\)，与已证二次反馈回归一致。
平面总面积为两排序 \(2\int H\)，所以
\[
 |U_2^{\rm plane}\setminus U_4^{\rm plane}|
 =\frac{241}{60}\delta^4+O(\delta^5).
 \tag{28}
\]

二维实际宏孔没有横向坐标或更高平衡支撑，两个相反短角就是两份上述模型。由 \(\sup_{s\le h_H}|F_w'(s)/f-1|=O(u^4)\)，真实逆容量相对线性高度误差为 \(O(u^4)\)；每个高度的绝对误差 \(O(\delta u^4)=O(\delta^5)\)，有限迭代的小常数保持此界。顶端截断不活跃。因此
\[
 |U_2\setminus U_4|
 =2e^2\left[\frac{241}{60}\delta^4+O(\delta^5)\right]
 \quad(n=2).
 \tag{29}
\]
以 (29) 代替保守 (9)，同样减少可删预算并保留 \(1-\theta\) 份额，得到 (D)。
若将全部高维短角也完整迁移，形式上会给 \((241/810)K_n\) 的族间九阶增益；本稿没有证明该尖高维迁移，实际全维结论只使用 (9)。

## 8. 一次显式保护下的全维微孔

本节从一个安全 \(U_4\) 宏观族出发；其参数可以是命题 A 的输出。记
\[
 \Lambda_0=\frac{B_n\kappa^3}{8n^2g}.
\]
由完整名义两支撑积分、(6) 和 (7)，
\[
 |U_4|/a=n\Lambda_0u^5g^{n-2}+O_n(u^6).
 \tag{30}
\]
三近面名义项更晚，且 \(U_0-U_4\) 总量为 \(O_n(au^6)\)，故本式不需新的九阶展开。

取新参数
\[
 \Lambda^\mu=\Lambda-a\Lambda_0+Cau ,
 \tag{31}
\]
其中 C 是最后选定的大固定安全常数；所有新 \(U_j^\mu,w^\mu,\rho^\mu,M^\mu\) 均据此重新定义。以下省略上标仅指这组新参数。

因为 \(U_4\) 是有限连续高度图的闭集，边界零测且内部体积等于总体积。固定当前 \(u\)，取
\[
 t'=t(1-u),\quad v'=2(t')^{1/n},\quad a'=(v')^n.
\]
只保留 \(q\in h\mathbb Z^n\) 中 host \(Q(q,h)\subset\operatorname{int}U_4\) 的微格，孔半宽为
\[
 w_\mu=\frac{hv'}{2(1+v')}.
\]
用不交基本胞 \(q+[-h/2,h/2)^n\) 计算有效体积。随着 dyadic \(h\to0\)，它趋于 \(|U_4|\)；可选第一个满足 \(h\le u^{20}\) 且有效体积至少 \((1-u)|U_4|\) 的间距。有限并集 H 因而满足
\[
 |H|=a|U_4|[1+O_n(u)].
 \tag{32}
\]
不要求 H 在每个宏观帽中均匀分布。

对每个微孔点 z，半径
\[
 r_\mu(z)=h-w_\mu-\|z-q\|_\infty
\]
含自身全孔，与其他格孔至多边界接触，并在 host 内。其缺陷率不超过
\((w_\mu/(h-2w_\mu))^n=t'<t\)。这些是原阈值的所有点严格见证。

固定已构造的 \(H\subset U_4\)，只作一次显式运算
\[
 \widetilde U=U_4\cup\mathcal B_w(U_4\setminus H).
 \tag{33}
\]
这就是宏孔保护。因为 \(H\subset\widetilde U\) 且
\(\widetilde U\setminus H\supset U_4\setminus H\)，反单调性给
\[
 \mathcal B_w(\widetilde U\setminus H)
 \subset\mathcal B_w(U_4\setminus H)\subset\widetilde U.
 \tag{34}
\]
实际其他周期的总目标正是 \((C_w\setminus\widetilde U)\cup H\) 的子集，
所以 (34) 是所需真实超解，而非只保护某个指定角点。
由 \(\mathcal B_w(U_4)\subset U_4\) 及 (6)，
\[
 |\widetilde U|\le|U_4|+C_nu|H|.
 \tag{35}
\]
这是对任意实际 H 的物理 \(L^1\) 核界，不含密度替代。

比较新旧预算。由 (16)、(18)、(30)–(32)，
\[
 \frac{D_{\Lambda^\mu}-D_\Lambda}{a}
 =-na\Lambda_0u^5g^{n-2}+nCau^6g^{n-2}+O_n(au^{10}),
\]
\[
 |H|/a=na\Lambda_0u^5g^{n-2}+O_n(au^6),\quad
 u|H|/a=O_n(au^6),\quad
 (|U_4^\mu|-|U_4|)/a=O_n(au^{10}).
\]
故
\[
 D_{\Lambda^\mu}+|H|-|\widetilde U|
 \ge D_\Lambda-|U_4|+a^2u^6[nCg^{n-2}-C_0]>0
 \tag{36}
\]
只需先在固定参数区间取得统一 \(C_0\)，再选 C，最后小 u 使 \(Cau\le1\)。没有误差常数反过来依未知 C 的选参循环。

在 \(C_w\setminus\widetilde U\) 内切取宏孔 A 的质量恰为 \(M-|H|>0\)，总孔为 \(A\cup H\)，总实际质量仍为 M。宏孔点用原长半径，(34) 保证自身 M 加外国实际质量严格小于 \(ar^n\)；微孔 host 在 \(U_4\subset\widetilde U\) 内，没有宏孔。不同宏胞外包孔分离，所以微见证不受外国孔污染。

相对原安全 \(U_4\) 族，(31) 给
\[
 M^\mu-M=2n a^2g^{n-1}\Lambda_0u^6+O_n(a^2u^7).
 \tag{37}
\]
固定实际总质量的周期比值因而增加
\[
 (2K_n+o_n(1))av^{n+6}.
 \tag{38}
\]
若 \(\kappa=8ng^2/3+O_n(u)\)，余项可写 \(O_n(av^{n+7})\)。
有限化取 \(N=\lceil u^{-(n+12)}\rceil\)，保留边界切孔，误差 \(O_n(a/N)\) 比 (38) 更晚；宏／微全部指定见证和外盒外全半径排除与第 6 节相同。
这证明命题 B。三维的 (38) 恰为 \(2av^9+o(av^9)\)，并可在命题 A 后再作 (31)，其交叉影响已经包含于上述参数稳定误差。

## 9. 已证增益与尚未识别的绝对九阶项

1. 对所有固定 \(n\ge2\)，真实 \(U_4\) 超解严格优于 \(U_2\)；全维保守正增益为 (B)，没有拟合旧误差。
2. 真实平面第三常数为 \(241/60\)，二维纯宏观族有 (D)。高维恢复全部平面尖系数还需额外证明，不在本稿签署。
3. 任意 H 的显式保护 (33) 供给全维真实微孔；三维产生新的九阶 \(+2\) 族间收益。
4. 二维保护量 \(O(u|H|)\) 对最终 halo 正好可能是九阶，且旧 \(cu^3\) 安全项同阶。不能把 (D) 直接加到已有最优八阶微孔族后便宣布其九阶系数，也不能把 (38) 的 \(O(av^9)\) 余项当作已知常数。

若要识别绝对九阶，仍须计算宏观 \(U_4\) 的全部 \(u^8\) 项（包括源中分面及高维非中央修正），并在二维控制实际微孔保护的九阶积分。本稿停在已证明的实际比较定理，不宣布 \(C_n\) 归一化九阶极限或某个九阶最优系数。

## 10. 实际本机证据

已实际读取 Terra 本机执行的：

* r060-computation/third_feedback.wl，SHA 605b30abe1d6faec96589ad322107f62d2145165701abfd8042919a1d6521fa8；
* third_feedback.log，SHA 881ad0c4e02ec336781ed82b01fe9e963a74544ffb40f219b62b511cfe31db67；
* fraction_polynomial.json，SHA 84de2e01c543432afdb91707a12901321cd47e4c5a29ab7ddb801d44b6200d40；
* cap_l1.wl，SHA d4882286b940e80fa77537b3cae7862fb136b140b894026ad2d38a47d1b51ee9；
* cap_l1.log，SHA 4f7b7e8a57148824297c957763ba5ba7552d9eaea7743865730f14650fccb264；
* cap_fraction.json，SHA cf57e2879fa084e229652706d2e001ecd0678683c213ce930373453de780a98d。

成功 Wolfram 15 日志及独立 Fraction 精确核对了 (27)、c1/c2 回归、平面尖系数的量纲换算、general-n 微孔增益 \(2K_n\) 及 n3=2；min 帽三个违反集为 False，Fubini 两条带合常数为 2。第一次替换次序未完全化简的 partial 日志不作认证。

另实际读取 single_corner.wl/log，SHA 分别为 30162e6f2bb12801b805735efc1825f0ed3099a567f3a9a4f1617e3584f3a709、a1592fc295d2476080aa6da2c6f82cc511440929b6941fab7abf478e5e4237c1；single_corner_fraction.py/json 的 SHA 分别为 501bde2fc529936336fc5a7aba45c0141b612b45aa99ddfb3931606b77617261、8a56978930a1af583881d779b812c9ab6888a0f1549b69ac1d35b03b3727d698。自由体积归一化、参数预算比例、比值增益和最优 \(\kappa\) 极限的四个通式残差均为零；\(\theta=1/2\) 给 \(\epsilon_n=c_nn^2/1458\)。两组有理参数只核代数及正预算比例，不认证该有限尺度满足全部几何阈值。

这些程序只核已明确定义的积分算子与有限代数，不证明邻孔清单、中央图、参数变动、宏／微实际集合、严格阈值或全半径有限化；对应解析证明均在正文中。

