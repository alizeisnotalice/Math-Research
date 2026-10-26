# audit_frame075：全 frame 曲率与 body 方向比值界的独立对抗审计

独立理论结论：两个提议均成立。下面从正产品Cauchy密度和固定纤维坐标重新证明，不把主代理的新推导当作已认证输入。主代理保留最终证明责任。

固定五维16条未单位化body方向、等权共同Cauchy时间，令
\[
B=[v_1\ \cdots\ v_{16}],\quad
V=\{(1,\pm1,\pm1,\pm1,\pm1)\},\quad
A=B^\dagger=B^T/16.
\]
对任意 \(\epsilon>0\)、\(x,u\in\mathbb R^5\)，截时势
\(H_\epsilon(x)=\int_\epsilon^\infty p_t(x)dt\) 正、有限、二次连续可微，并满足
\[
\boxed{D^2\log H_\epsilon(x)[u,u]\ge-\frac{|u|^2}{8\epsilon^2},}
\tag{1}
\]
\[
\boxed{\log H_\epsilon(x+u)\ge\log H_\epsilon(x)
+\nabla\log H_\epsilon(x)\cdot u-\frac{|u|^2}{16\epsilon^2}.}
\tag{2}
\]
对每个 \(t>0\)，另有
\[
\boxed{\|\nabla\log p_t(x)\|_1\le\frac1t,\qquad
e^{-\|u\|_\infty/t}\le\frac{p_t(x+u)}{p_t(x)}\le e^{\|u\|_\infty/t}.}
\tag{3}
\]
类似地 \(\|\nabla\log H_\epsilon\|_1\le1/\epsilon\)，因而H的位移比也有相同形式、把t替换为epsilon的界。这些结论覆盖所有实空间点，包括零/重复锥节点和完整Green的S1/S2奇异集合；不依赖扇形表达式分母的非退化条件。

实际有限资格只运行一次标准库Fraction矩阵检查，PID64058，exit0，记录耗时0.03653487502015196秒。没有H、梯度、Green、log、exp、Gamma或Arb求值，没有重放162父盒。它验证代数身份而不替代下述解析证明。这些是供后续认证使用的组件，不构成新的全域g、完整I、C5或Cn突破。

## 1. 非方阵推前与纤维测度的正确常数

16个符号方向的坐标交叉和为零，平方和为16，故
\[
BB^T=16I_5,\quad BA=I_5,\quad A^TA=I_5/16.
\tag{4}
\]
令K=ker B，维数11，使用K上Euclidean Hausdorff/Lebesgue测度 \(dw\)。它由任意正交基识别为 \(\mathbb R^{11}\) 的普通Lebesgue测度，与基的选择无关。

映射 \((x,w)\mapsto Ax+w\) 是从 \(\mathbb R^5\oplus K\) 到 \(\mathbb R^{16}\) 的线性同构。因为im A与K正交，且A的五个奇异值均为1/4，其16维Jacobian为 \(4^{-5}=1/1024\)。等价地，非方阵B的coarea Jacobian是
\(J_B=\sqrt{\det(BB^T)}=16^{5/2}=1024\)。

写
\[
s_t(z)=\frac{t}{\pi(t^2+z^2)},\qquad
f_t(z)=\prod_{j=1}^{16}s_t(z_j).
\]
因此正推前密度的纤维公式是
\[
\boxed{p_t(x)=\frac1{1024}\int_Kf_t(Ax+w)dw.}
\tag{5}
\]
该常数不是 \(1/16\)、\(1/\det B\) 或单位常数。B非方阵，不能对其直接写普通行列式。

式(5)先由变量替换或coarea给出a.e.版本。下一节将证明右侧对每个x有限且连续，从而它是推前分布的canonical连续密度；不能只凭a.e.版本未经解释就在异常纤维上求导。

## 2. 每条纤维有限、处处严格正及完整时间尾

从V中固定选择五列C：全+向量及四个只翻转一个尾坐标的向量，另外11列记为D。列减法给 \(\det C=16\)。以剩余坐标 \(y\in\mathbb R^{11}\) 参数化纤维：
\[
z_C=C^{-1}(x-Dy),\qquad z_D=y.
\]
切向矩阵为把 \(-C^{-1}D\) 与 \(I_{11}\) 竖直堆叠所得的16×11矩阵M。其Gram行列式由Sylvester恒等式或直接消元给
\[
\det(M^TM)=\frac{\det(BB^T)}{(\det C)^2}=4096,
\]
故 \(dw=64\,dy\)。结合(5)，得到对每个实x有效的正积分
\[
p_t(x)=\frac1{16}\int_{\mathbb R^{11}}
\prod_{j=1}^{5}s_t\big((C^{-1}(x-Dy))_j\big)
\prod_{j=1}^{11}s_t(y_j)dy.
\tag{6}
\]
被积函数严格正；余下11个变量为概率测度，而每个选中密度至多 \(1/(\pi t)\)。所以
\[
0<p_t(x)\le\frac1{16\pi^5t^5}\quad(t>0,x\in\mathbb R^5).
\tag{7}
\]
同一个支配函数也给(6)的连续性。式(5)因而处处有限且等于该连续密度，消除了异常纤维问题。

对于任意固定 \(\epsilon>0\)，(7)给
\[
0<H_\epsilon(x)\le\int_\epsilon^\infty\frac{dt}{16\pi^5t^5}
=\frac1{64\pi^5\epsilon^4}<\infty.
\tag{8}
\]
时间下端是严格正epsilon，不涉及0处奇性；上端∞由完整 \(t^{-5}\) 尾控制。积分域端点不随x移动，因此空间求导不会产生时间端点项。此处只需一个粗解析峰值，不消费或重算057、旧Gamma或Green常数。

## 3. 固定右逆平移与微分支配

对任意u，取固定提升 \(a=Au\)，于是Ba=u且
\(\|a\|_2^2=|u|^2/16\)。沿同一纤维坐标w，
\[
p_t(x+su)=\frac1{1024}\int_Kf_t(Ax+w+sa)dw.
\tag{9}
\]
这是同一K与同一测度dw，不是随x改变的纤维参数化或归一常数。

单坐标log-score及二阶log导数为
\[
\ell_t'(z)=-\frac{2z}{t^2+z^2},\qquad
\ell_t''(z)=\frac{2(z^2-t^2)}{(t^2+z^2)^2}.
\]
它们满足全局界
\[
|\ell_t'(z)|\le\frac1t,\qquad
|\ell_t''(z)|\le\frac2{t^2},\qquad
\ell_t''(z)\ge-\frac2{t^2}.
\tag{10}
\]
第一式来自 \(2|z|t\le z^2+t^2\)；最后一式还可准确写成
\[
\ell_t''(z)+2/t^2
=\frac{2z^2(z^2+3t^2)}{t^2(t^2+z^2)^2}\ge0.
\]

在任意固定实z上，定义
\[
\sigma_a(t,z)=\sum_ja_j\ell_t'(z_j),\quad
\tau_a(t,z)=\sum_ja_j^2\ell_t''(z_j).
\]
于是
\[
|\sigma_a|\le\|a\|_1/t,\quad
|\tau_a|\le2\|a\|_2^2/t^2,
\quad \tau_a\ge-2\|a\|_2^2/t^2.
\tag{11}
\]
沿任何有限区间 \(|s|\le R\) 积分第一阶score界，得到
\[
f_t(z+sa)\le e^{R\|a\|_1/t}f_t(z)
\le e^{R\|a\|_1/\epsilon}f_t(z)\quad(t\ge\epsilon).
\tag{12}
\]
而
\[
\partial_sf_t(z+sa)=f_t(z+sa)\sigma_a(t,z+sa),
\]
\[
\partial_s^2f_t(z+sa)
=f_t(z+sa)(\sigma_a^2+\tau_a)(t,z+sa).
\]
故零、一、二阶导都被常数乘原 \(f_t(z)\) 支配，常数只依赖a、R、epsilon；比如二阶可取
\(e^{R\|a\|_1/\epsilon}(\|a\|_1^2+2\|a\|_2^2)/\epsilon^2\)。该支配同时适用于整个K和完整 \([\epsilon,\infty)\)。由(8)原被积函数总质量有限，因此DCT允许两次微分。

多变量的局部支配同样成立：对h在一个有界邻域，\(\|Ah\|_1\le4\|Ah\|_2=|h|\)；混合二阶导由相应有界scores乘积和单坐标二阶项控制。故H确实为C²，导数连续。固定t的p亦为C²；其实可继续证明更高光滑性，但本审计不需要。

所有支配使用的是有界score，完全没有对重尾变量要求 \(\mathbb EZ_j\)、\(\mathbb EZ_j^2\) 或任意绝对矩。尤其不能把(11)替换为关于z的多项式矩支配。

## 4. 联合纤维—时间概率与全 frame 曲率

由(5)、(8)，对固定x定义正概率测度
\[
d\mu_x(t,w)=\frac{f_t(Ax+w)}{1024H_\epsilon(x)}dt\,dw,
\quad(t,w)\in[\epsilon,\infty)\times K.
\tag{13}
\]
归一化恰为1。上一节合法微分给
\[
\frac{D_uH}{H}=\mathbb E_{\mu_x}\sigma_a,\qquad
\frac{D_u^2H}{H}=\mathbb E_{\mu_x}(\sigma_a^2+\tau_a).
\]
所有scores有界，方差确实有限。于是
\[
D_u^2\log H
=\mathbb E\tau_a+\operatorname{Var}(\sigma_a)
\ge-2\|a\|_2^2\mathbb E(t^{-2})
\ge-\frac{2\|a\|_2^2}{\epsilon^2}
=-\frac{|u|^2}{8\epsilon^2}.
\]
混合带来的是非负方差项，方向不能反用；无需假设归一后的坐标彼此独立，也无需假设t与纤维变量独立。

对 \(h(s)=\log H(x+su)\)，积分有界二阶导：
\[
h(1)-h(0)-h'(0)=\int_0^1(1-s)h''(s)ds
\ge-\frac{|u|^2}{16\epsilon^2}.
\]
因此(1)、(2)中的1/8与1/16分别是Hessian常数和积分罚项，二者不能混淆。

补充：若采用任意其它固定提升a满足Ba=u，则 \(a=Au+k\)、k∈K，且Au⊥k，故
\(\|a\|^2=|u|^2/16+\|k\|^2\)。因此A是当前“逐坐标统一二阶下界+固定提升”方法的最小范数选择。这不证明(1)是真实H的最优曲率常数；条件方差和更精细的score分布仍可能给进一步改进。

## 5. 独立的log-convexity交叉证明及其陷阱

对每个固定t≥epsilon与w∈K，考虑
\[
F_{t,w}(x)=e^{\|Ax\|^2/\epsilon^2}f_t(Ax+w).
\]
其log在方向u的Hessian为
\[
2\|Au\|^2/\epsilon^2+\tau_{Au}(t,Ax+w)
\ge2\|Au\|^2(\epsilon^{-2}-t^{-2})\ge0.
\]
所以每个 \(F_{t,w}\) 都log-convex。对0<lambda<1逐点用该性质，然后用Hölder：
\[
\int F_{t,w}(\lambda x+(1-\lambda)y)dt\,dw
\le\int F_{t,w}(x)^\lambda F_{t,w}(y)^{1-\lambda}dt\,dw
\le\left(\int F_{t,w}(x)\right)^\lambda
\left(\int F_{t,w}(y)\right)^{1-\lambda}.
\]
所有积分对固定x,y有限，因共同修正因子是有限常数且H有限。于是
\[
e^{|x|^2/(16\epsilon^2)}H_\epsilon(x)
\quad\text{是log-convex函数。}
\tag{14}
\]
结合已证C²即与(1)等价。这给方差推导的独立交叉，使用log-convex函数积分的Hölder方向，不是把log-concavity的边际定理用反。

必须保留共同的epsilon修正。对单个t，\(e^{\|Ax\|^2/t^2}p_t(x)\)的log-convexity也成立；但直接对t积分所得不是(14)，因为指数因子依赖t，不能从积分号中取出。先将每个t修正为共同epsilon尺度才闭合完整H的证明。

## 6. 单body方向微分、全部cube顶点与比值界

固定v_j=Be_j。其最小提升Av_j与坐标提升e_j只差K内向量：
\(e_j-Av_j\in K\)。在(5)内对w作这个常向量的平移，Hausdorff测度不变，所以对每个s有
\[
p_t(x+sv_j)=\frac1{1024}\int_Kf_t(Ax+w+se_j)dw.
\tag{15}
\]
因此可以只微分第j个Cauchy坐标。第3节的支配同样适用，给
\[
|D_{v_j}p_t(x)|\le\frac1{1024}\int_K f_t(Ax+w)\frac1t dw
=\frac{p_t(x)}t.
\tag{16}
\]
注意这里v_j未单位化；若先把v_j单位化，常数必须相应调整。p_t处处正，所以
\(|v_j\cdot\nabla\log p_t|\le1/t\)。

\(V\cup(-V)=\{-1,1\}^5\)，恰为32个cube顶点。对任意g∈R5，
\[
\|g\|_1=\max_{s\in\{-1,1\}^5}s\cdot g
=\max_{v\in V}|v\cdot g|.
\]
由此得到(3)的梯度界。再沿线段积分并用对偶范数，
\[
|\log p_t(x+u)-\log p_t(x)|
\le\int_0^1\|\nabla\log p_t(x+su)\|_1\|u\|_\infty ds
\le\|u\|_\infty/t.
\]
指数化即得到两侧比值。没有额外的5、√5、10或16因子。

完整截时势可直接积分(16)：
\(|D_{v_j}H_\epsilon|\le\int_\epsilon^\infty p_t/t\,dt\le H_\epsilon/\epsilon\)，从而得到H的同型梯度及比值界。对任何固定正时间带[a,b]，同理可用a；无需截掉空间奇异集合。

比值界对所有位移成立，但它相对旧asinh型比值包络并非全位移统一更小：本界指数随 \(\|u\|_\infty/t\) 线性增长，asinh型在大参数时可只对数增长。后续band计算应在实际参数上比较，或取两种合法上界的最小值。本审计未运行band优化，也未断言067所有参数都得到同等改进。

## 7. 尾坐标曲率、旧056比较及16角点闭盒证书

限制u=(0,v)，完整frame罚项为
\(v^TQ_{\rm full}v/\epsilon^2\)，其中
\[
Q_{\rm full}=I_4/16.
\]
056使用 \(Q_s=I_4/4-ss^T/36\)，s∈{−1,1}4。新旧之差准确为
\[
Q_s-Q_{\rm full}
=\frac{11}{144}I_4+\frac{4I_4-ss^T}{36}
\succeq\frac{11}{144}I_4>0,
\tag{17}
\]
因为 \(|s|^2=4\)，Cauchy给第二项半正定。故新罚项严格小于每个旧Qs的罚项，非零位移的统一改善至少 \(11|v|^2/(144\epsilon^2)\)。这也严格改进其逐点最小psi。数值检查独立验证了16个旧基的逆矩阵、Qs身份以及(17)的精确矩阵分解，没有新基搜索。

若真实中心梯度g已知，(2)的右侧是单一凹二次函数。若只持有保存的梯度区间 \(g_j\in[g_{j,L},g_{j,U}]\)，安全下界为
\[
\Phi(v)=\Delta_L+\sum_{j=1}^4\min(g_{j,L}v_j,g_{j,U}v_j)
-\frac{|v|^2}{16\epsilon^2},
\quad\Delta_L\le\log(H_L/g_0),\quad0<H_L\le H_\epsilon(x_0).
\tag{18}
\]
这是全局凹的分段二次函数：两个仿射函数的最小值凹，有限和仍凹，再减正二次保凹。使用梯度区间时不能字面称它为单一二次，但不影响角点结论。

对任意四维闭矩形盒，每个点都是其最多16个角点的凸组合。凹性给 \(\Phi(v)\ge\sum\lambda_i\Phi(v_i)\ge\min_i\Phi(v_i)\)。因此
\[
\min_{16\text{ corners}}\Phi>0
\quad\Longrightarrow\quad H_\epsilon>g_0\text{ 在整个闭盒上成立。}
\tag{19}
\]
盒内穿越坐标符号面、投影重合面、S1/S2都不影响这一推导。anchor甚至无需位于盒内；仅需保存的H/梯度确实属于该anchor且位移端点正确。退化宽度和重复角点可去重。与056的逐象限psi不同，这里不需要每坐标加零后的81点覆盖。

本轮只证明该新判据，未消费162盒的HL、gradient、Delta，未声称那些盒的新余量、树计数或全局覆盖已经计算。

## 8. 对抗检查与适用边界

以下潜在失败点均在证明中显式排除：

- 非方阵Jacobian：使用1024的coarea因子，另用固定C给出64的纤维参数化体积，二者组合为1/16，与普通卷积分组一致。
- 异常纤维：先以固定可逆5列证明每个x的有限性、严格正性与连续性，再求log和微分。
- 重尾支配：只用score及二阶log导的统一界，不使用任何Cauchy矩。
- 归一混合：对联合(t,w)构造一个概率测度，方差符号为正；不把posterior分布误当原产品独立概率。
- 时间端点：epsilon必须严格正，上尾由t^-5完整支付；本曲率常数不能令epsilon=0用于可能无穷的G。
- 可移动纤维：所有积分在固定K上，任何两种lift只差K内平移，测度不变。
- log-convex积分：使用共同epsilon修正与Hölder；直接混合t依赖修正项不足以证明H结论。
- 空间退化：零、重复节点和S1/S2无须剔除；它们只影响其它核表示或完整G，不影响正时间产品卷积。
- 闭盒角点：区间梯度下界仍全局凹，因此16角点充分；这不是未经检查地套用056的81点逻辑。

没有发现反例或遗漏导致(1)–(3)失效。以上说明限于同一16方向、共同时间t、等权Cauchy模型；方向尺度、权重、不同时间或低维矩阵改变时，应重新计算frame矩阵、lift范数与单坐标score常数。

## 9. 唯一有限资格、文件锁与未完成事项

登记SHA `9b9e181158d082a7e3430e08f99dbb8d13f63184d43fa9263ea5248489f25250`；源码SHA `0abba4a9d759ef2e86780b449d543dbb63960ec3e1146ff4d4e31fa8fe1e2c87`。仅一次标准库Python执行，真实PID64058，exit0，记录时间0.03653487502015196秒，整个工具执行完成于一次短调用内，远低于180秒总wall限制。内部160秒alarm包含输入读取、计算与结果落盘阶段，未触发；无等待启动的armed阶段。

固定检查为：BBᵀ=16I、BA=I、AᵀA=I/16、AB为正交投影、全部16个 \(e_j-Av_j\) 属于ker B、±V覆盖全部32cube顶点、固定det C=16、coarea常数1024、11维参数化Gram行列式4096、纤维体积因子64，以及16个旧056基与新曲率的准确比较。完整矩阵与有理输出保存在results；不把这些有限代数检查当作测度论形式化验证。

读取过本clone AGENTS，只新增`audit_frame075/`；旧068的全部6文件hash与056只读报告/登记在计算前后均未变化。未修改旧052/056/068，未启动或干预070核，未使用MATLAB，未提交或推送。最小白名单为登记、检查源码、完整结果、本报告、manifest；没有复制旧树或核缓存。

后续仍需另登记实际消费：在现有HL/梯度/Delta上重做16角点余量，或将(3)代入实际band积分误差。它们未在本轮运行，不能从理论改善直接写成更少树节点、新全域g或更优Cn。主证明责任及是否接入后续批仍由主代理保留。
