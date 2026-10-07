# 原 winner 的有界 receiver test-energy：准确归约、两类障碍与局部配对接口

2026-10-07。独立研究稿。固定物理尺度为1、固定入口参数 \(c>0\)，任意非负 \(f\in L^1(\mathbb R^n)\)，\(W=\int f\)。本文保留真实有序核与 receiver winner；不增加 cube 主账费用，也不声称满足 FIRST/CP/GP/历史门。使用 L03 的“精确常数链与浮点诊断分开”工作流；技能不提供本文定理。耗散恒等式引用已证明并独审的 [连续有界熵稿](ordered_continuous_bounded_entropy_20261007.md)，这里另核 receiver 配对的常数。

## 1. 有限未来与确切 flow

令
\[
v_z=T_{1-e^{-z}}^c f,\qquad
T_r^c=\prod_{i=1}^n[(1-r)I+rG_{c,i}],\qquad
J_z=\sum_i(G_{ce^{-z},i}-I).
\tag{1}
\]
这是原过程，\(v'_z=J_zv_z\)。每个轴核正、对称、保质量，Fourier 乘子为
\[
g_d(\xi^2)=\frac1{1+dB(\xi^2)},\qquad
B(v)=\frac{v}{\log(1+v)}-1,\quad B(0)=0.
\tag{2}
\]
给定有限 \(Z>0\)，\(M_Z(x)=\max_{0\le z\le Z}v_z(x)\)，\(E=\{M_Z>\lambda\}\)。下文写 \(M=M_Z\)。有限 mask 展开保证存在同一满测集，其上轨道是连续函数，\(M\le\sum_A G_Af<\infty\)；故 \(|E|<\infty\)。选择最早最大点 \(\tau(x)\in[0,Z]\) 是可测的：连续紧参数最大值可用可数有理点描述，最早 argmax 的子水平集亦可由紧区间最大值描述。任何可测最大点都满足以下公式。

对任意可测 \(S\subset E\)，设
\[
\theta=\mathbf1_S\lambda/M,\qquad
H_z=\theta\mathbf1_{\{z<\tau\}},\qquad
t_z=H_z(1+v_z/\lambda).
\tag{3}
\]
由于 \(v_z\le M\)，有 \(0\le t_z\le1+\lambda/M<2\)，且 \(t_z=0\) 于 \(S^c\)。这里没有把 receiver 选择改为独立时间选择。点态微积分及 Fubini 给
\[
\lambda|S|-\int\theta f
=\int_0^Z\!\int H_zJ_zv_z.
\tag{4}
\]
绝对可积由 \(H\le1\)、\(\int_0^Z\|J_zv_z\|_1dz\le2nZW\) 保证。时间端点 gate 的值不影响积分；若 \(\tau=0\)，该点没有 flow，全部费用在初始项。

## 2. 正 Bregman 及 commutator：常数核对

记 \(h(v)=\lambda/(\lambda+v)\)、\(L(v)=\log(1+v/\lambda)\)。点态恒等式
\[
h(v_z)J_zv_z=\lambda J_zL(v_z)+b_z,\qquad b_z\ge0
\tag{5}
\]
的每个 jump 被积式是
\[
\lambda\left[a-1-\log a\right],\qquad
a=\frac{\lambda+v_z(x+h)}{\lambda+v_z(x)}>0.
\]
令 \(D=\int_0^Z\int b_z\)。保质量及原有界熵恒等式给
\[
D=\int\Phi_\lambda(f)-\int\Phi_\lambda(v_Z)
=\frac{\lambda}{2}\int_0^Z\sum_i\int dx\,G_{ce^{-z},i}(dh)\,
\frac{(v_z(x+h)-v_z(x))^2}{(\lambda+v_z(x+h))(\lambda+v_z(x))}
\le W,
\tag{6}
\]
其中 \(0\le\Phi_\lambda(v)=v-\lambda\log(1+v/\lambda)\le v\)。没有使用 \(L\log L\) 假设。

乘以 \(t_z\) 后，式(4)等于 \(\int t_zb_z+C_S\)，其中
\[
C_S=\lambda\int_0^Z\!\int t_zJ_zL(v_z)
=-\frac{\lambda}{2}\int_0^Z\sum_i\int dx\,G_{ce^{-z},i}(dh)\,
\Delta t_z\,\Delta L(v_z).
\tag{7}
\]
所有项绝对可积：\(t_z\le2\)，\(\lambda L(v_z)\le v_z\)。故
\[
\lambda|S|\le\int\theta f+2D+C_S.
\tag{8}
\]
保留 signed \(C_S\) 比先取绝对值更强。

对 \(a,b>0\)，令 \(s=\log(a/b)\)，则
\[
(\log(a/b))^2\le\frac{(a-b)^2}{ab},
\tag{9}
\]
因为右边是 \(4\sinh^2(s/2)\ge s^2\)。定义未缩放的 test-energy
\[
\mathcal T_S=\int_0^Z\sum_i\int dx\,G_{ce^{-z},i}(dh)\,(\Delta t_z)^2.
\tag{10}
\]
对称 Cauchy 及(6),(9)准确给
\[
|C_S|\le\sqrt{D\cdot\frac{\lambda}{2}\mathcal T_S}
\le\sqrt{W\cdot\frac{\lambda}{2}\mathcal T_S}.
\tag{11}
\]
因此 \(S=E\) 时
\[
\lambda|E|\le3W+\sqrt{W(\lambda/2)\mathcal T_E}.
\tag{12}
\]
常数中没有缺失因子2。粗界
\[
\mathcal T_S\le8nZ|S|
\tag{13}
\]
来自 \((\Delta t)^2\le4(\mathbf1_S(x)+\mathbf1_S(x+h))\)。它只复得 \(O(nZ)W\) 的水平集费用。

可先付 \(E_0=\{f>\lambda/2\}\)，\(\lambda|E_0|\le2W\)。余集
\[
E_{\rm new}=\{f\le\lambda/2,\ M>\lambda\}
\]
满足 \(\int\theta f\le(\lambda/2)|E_{\rm new}|\)，故
\[
\frac{\lambda}{2}|E_{\rm new}|\le2D+C_{E_{\rm new}}.
\tag{14}
\]
以下反例说明这个初始 gap 仍不足以得到 test-energy 单独的 polylog 界。

## 3. 全水平集的近并列真实 winner 障碍

取 \(c=1,Z=1,\lambda=1/2,K\) 为整数且 \(K\ge64\)。在完整 \(n\)-torus 的归一化 Haar 测度上，令
\[
f(x)=1+\varepsilon\cos\vartheta,\qquad
\vartheta=K\sum_i x_i,\qquad0<\varepsilon\le10^{-3}.
\tag{15}
\]
来源高度相关；独立的是算子轴向 jump，不是来源坐标。原响应准确为
\[
v_z=1+\varepsilon m(z)\cos\vartheta,\qquad
m(z)=[e^{-z}+(1-e^{-z})g_1(K^2)]^n.
\tag{16}
\]
其中 \(m'(z)=-n[1-g_{e^{-z}}(K^2)]m(z)<0\)。因此除零测的 \(\cos\vartheta=0\) 外，原唯一 winner 为
\[
\tau=0\quad(\cos\vartheta>0),\qquad
\tau=1\quad(\cos\vartheta<0).
\tag{17}
\]
所有点属于 \(E\)，\(W=1,|E|=1\)。在 \(0<z<1\)，
\[
t_z=\mathbf1_{\{\cos\vartheta<0\}}
\frac{3/2+\varepsilon m(z)\cos\vartheta}
{1+\varepsilon m(1)\cos\vartheta}.
\tag{18}
\]
即便 \(\varepsilon\downarrow0\)，此 test 仍趋 \((3/2)\mathbf1_{\{\cos\vartheta<0\}}\)，不是趋常数。

这里所有常数只用 **\(c=1\)**。\(z\le1\) 保证 \(e^{-z}>1/4\)：\(e<4\)。又 \(e>8/3\) 且 \((8/3)^9>4097\)，故
\[
B(4096)>4096/9-1>450,\qquad
g_{e^{-z}}(j^2K^2)<1/100\quad(j\ne0).
\tag{19}
\]
原 \(G_d\) 的谱位于 \([0,1]\)。对相位函数定义
\[
d_d(q)^2=\int dx\,G_{d,i}(dh)(q(x+h)-q(x))^2
=2\langle q,(I-G_d)q\rangle.
\]
每个坐标的相位响应相同。由(19)，\(d_d(q)^2\ge2(99/100)\operatorname{Var}(q)\)。

对 \(q_0=(3/2)\mathbf1_{\{\cos<0\}}\)，
\[
d_d(q_0)\ge(3/2)\sqrt{99/200}>(3/2)(7/10)=21/20.
\]
式(18)与 \(q_0\) 的差只支撑半圆，其幅度至多
\[
b=\frac{(5/2)\varepsilon}{1-\varepsilon}\le5/1998<1/200.
\]
由于 \(0\le G_d\le I\) 的谱，\(d_d(t-q_0)\le b\)。Dirichlet 半范数三角不等式给
\[
d_d(t_z)^2>(21/20-1/200)^2>1,\qquad
\boxed{\mathcal T_E\ge n}.
\tag{20}
\]
它否定所有来源上的 \(\mathcal T_E\le C\,\mathrm{polylog}(n)|E|\)。但这个例子 \(\lambda|E|=W/2\)，原弱费用已经由初始项支付。它不是弱端点反例。

## 4. 初始深 gap 也不能单独支付 test-energy

仍取同一原算子、\(c=1,Z=1,\lambda=1/2,K\ge64\)，改为
\[
f=1+\cos\vartheta,\qquad W=1,\qquad
S=E_{\rm new}=\{\cos\vartheta\le-3/4\}.
\tag{21}
\]
以下对所有 \(n\ge8\) 成立。因为 \(e^{-1/2}<2/3\)，(19)给
\[
m(z)\le(7/10)^n<1/16\quad(1/2\le z\le1).
\tag{22}
\]
故 \(M=v_1\ge15/16>\lambda\)，并且在 \(S\) 上唯一 winner 是1、\(f\le1/4=\lambda/2\)。集合比例
\[
p=|S|=\arccos(3/4)/\pi,\qquad1/5<p<1/4.
\tag{23}
\]
两界来自 \(\cos(\pi/5)=(1+\sqrt5)/4>3/4\) 及 \(\cos(\pi/4)=\sqrt2/2<3/4\)。

对 \(q_0=(3/2)\mathbf1_S\)，有
\[
d_d(q_0)^2\ge2(99/100)(9/4)p(1-p)
>2(99/100)(9/4)(3/20)>(4/5)^2.
\]
在 \(1/2\le z<1\)，真实 \(t_z=\mathbf1_S(3/2+m(z)\cos\vartheta)/(1+m(1)\cos\vartheta)\)，其与 \(q_0\) 的幅度差小于
\[
\frac{(5/2)(1/16)}{1-1/16}=1/6.
\]
其支撑比例 \(p<1/4\)，故 \(d_d(t-q_0)\le\sqrt{2p}/6<1/6\)。于是
\[
d_d(t_z)^2>(4/5-1/6)^2>1/3,\qquad
\boxed{\mathcal T_{E_{\rm new}}\ge n/6}.
\tag{24}
\]
这里 \(|E_{\rm new}|<1/4\)、\(W/\lambda=2\)。因此
\[
\mathcal T_{E_{\rm new}}\le C\,\mathrm{polylog}(n)
\bigl(|E_{\rm new}|+W/\lambda\bigr)
\tag{25}
\]
这一单独能量合同同样不成立。没有让小扰动 selector 冒充任意抽象 gate；输入仅有一个真实有限 Fourier mode，winner 由(16)严格决定。

## 5. 两个障碍的有限 \(L^1(\mathbb R^n)\) 提升

torus 本身不是 \(L^1(\mathbb R^n)\) 来源。以下提升针对上面 test-energy 合同，不声称实际 geom 门。

令 \(Q_L=[-L/2,L/2]^n\)，\(L\) 取 \(2\pi\) 的整数倍，\(f_L=\mathbf1_{Q_L}(1+A\cos\vartheta)\)，其中 \(A=\varepsilon\) 或1。则 \(W_L=L^n\)、\(0\le f_L\le2\)。在入口 \(c=1\)，\(G_1=w\)，其一维 Fourier 乘子是 \(\log(1+\xi^2)/\xi^2=\int_0^1(1+t\xi^2)^{-1}dt\)。因此 \(w\) 是尺度 \(\sqrt t\le1\) 的 Laplace 混合，\(\Pr(|Y|>R)\le e^{-R}\)。

原 \(T_{1-e^{-z}}\) 的每轴是 holding 与 \(w\) 的混合。若 \(x\) 距 \(Q_L\) 边界至少 \(R\)，在所有 \(0\le z\le1\) 一致有
\[
0\le v_z^{\rm per}(x)-v_z^L(x)\le2ne^{-R}.
\tag{26}
\]
若某一坐标在 \(Q_L\) 外的距离大于2，进入来源盒必须在该坐标跳过至少2，故所有 \(z\) 有 \(v_z^L(x)\le2e^{-2}<1/2\)。于是
\[
E_L\subset Q_L+[-2,2]^n,\qquad |E_L|\le(L+4)^n.
\tag{27}
\]
flat 来源在 \(Q_L\) 内处处大于 \(\lambda\)，故 \(|E_L|/L^n\to1\)。gap 来源在内部的初始低集恰是(21)的相位集合；由(26)及 \(v_1^{\rm per}\ge15/16\)，其 bulk 最终都越过 \(\lambda\)，而外部只占(27)的边界条带，故 \(|E_{{\rm new},L}|/L^n\to p\)。

还需检查 winner/test，而不只检查输入质量。固定 \(n\) 后，bulk 轨道由(26)在紧时间区间一致趋向周期轨道。周期最大点在正/负 cosine 区分别唯一为0/1；紧参数 argmax 稳定性给 \(\tau_L\to\tau\)，除零测相位外成立。对 gap 集，cosine 至少离零 \(3/4\)，稳定性更直接。因 \(t\le2\)，对固定有界 jump、\(0<z<1\) 及避开相位边界的 bulk 周期单元，可用支配收敛。先丢弃宽度固定的空间边界条带，按周期单元取相位平均，再令 bulk 误差趋零；最后以正 Tonelli 放开 jump 截断及时间/相位小排除集，得到
\[
\liminf_{L\to\infty}\frac{\mathcal T_L}{L^n}
\ge\mathcal T_{\rm per}.
\tag{28}
\]
这个证明只需 jump 核正与总质量1，不假定轴核有 Lebesgue 全维密度。固定 \(n\) 时先让 \(L\to\infty\)，之后才让 \(n\) 增大，避免把极小 winner gap 的尺寸偷当统一常数。

由(20),(24),(27),(28)，可在每个 \(n\ge8\) 选择有限大盒，得到真正有限 \(L^1\cap L^\infty\) 来源族，使 full test-energy/receiver-volume 为 \(\Omega(n)\)，并且 gap test-energy/\((|E_{\rm new}|+W/\lambda)\) 亦为 \(\Omega(n)\)。例如极限下后者至少 \((n/6)/(p+2)>2n/27\)。这是解析提升；没有把周期上界直接当作欧氏上界。

## 6. 仍可用的局部配对，而非强 test-energy 合同

定义逐时间真实熵速率 \(d(z)\) 与 scaled test 速率 \(a_S(z)\)：
\[
d(z)=\frac{\lambda}{2}\sum_i\int G_{ce^{-z},i}(dh)dx\,
\frac{(\Delta v_z)^2}{(\lambda+v_z(x+h))(\lambda+v_z(x))},
\qquad
a_S(z)=\frac{\lambda}{2}\sum_i\int G_{ce^{-z},i}(dh)dx\,(\Delta t_z)^2.
\]
式(7),(9)给更准确的
\[
|C_S|\le\int_0^Z\sqrt{d(z)a_S(z)}\,dz.
\tag{29}
\]
任何可测的时间/坐标/jump/空间划分 \(\mathcal P\) 还给
\[
|C_S|\le\sum_{P\in\mathcal P}\sqrt{D_P A_P},
\qquad
\sum_PD_P=D\le W.
\tag{30}
\]
细化划分只会改善此界；来源熵预算仍仅一份 \(W\)。这不是对 conditional weak 常数做平均。

更强可先保留 \((-\Delta t\,\Delta L)_+\)，因为其他方向的 signed commutator 是负费用。特别，若 \(H\) 在一条边不变，则 \(\Delta t=H\Delta v/\lambda\)，该边 commutator 非正。未付的是 winner/level 选择发生变化且与响应差方向相反的边；没有证明这些边的一般空间支付。

两个反例恰说明必须保留配对。对 \(v_z=1+A m(z)\cos\vartheta\)，\(0<A\le1,\lambda=1/2\)，真实正谱性及 \(\lambda+v\ge\lambda\) 给
\[
d(z)\le\frac{nA^2m(z)^2}{2\lambda},\qquad
a_S(z)\le4\lambda n|S|\le4\lambda n.
\]
由(19)及原微分式，\(n\int_0^1m(z)dz\le100/99\)。于是
\[
\int_0^1\sqrt{d(z)a_S(z)}dz
\le\frac{100\sqrt2}{99}A
<\frac32 A\le\frac32 W.
\tag{31}
\]
full near-tie 版本更得到 \(O(\varepsilon W)\)。gap 版本的熵消耗集中在早期 \(O(1/n)\) 时窗，而测试能量在晚期仍粗糙；全时间 Cauchy 把它们错误配对，数值中上界会随 \(n\) 增大，实际(29)却保持常数。

式(31)仅是这两个原单相位族的严格守卫，**不是一般 \(f\) 的新费用**。下一步真正最小待证接口是对实际 winner 的 oriented commutator，或其时间/边局部配对(29)/(30)，给 \(\mathrm{polylog}(n)W\) 的来源一次预算。单独 test-energy 体积界已经被排除。当前没有证明该配对合同；也没有由此闭合 \(A_{\rm ord}\) 或原时变 cube 账。

## 7. 三轮专属登记、执行与误差范围

先写 [预登记](ordered_capped_test_energy_registration_20261007.json)，再运行新 [脚本](ordered_capped_test_energy_guard_20261007.py) 一次；[结果](ordered_capped_test_energy_results_20261007.json) 保存18组诊断及18个 NPZ/hash。固定 \(c=1,Z=1,\lambda=1/2\)，三轮 \((n,K,\varepsilon)=(8,64,10^{-3}),(32,128,10^{-4}),(128,256,10^{-5})\)，每轮同时 flat/gap。相位1024/2048/4096、时间 Simpson 256/512/1024，原 \(g_{e^{-z}}(j^2K^2)\) FFT；不是 reset kernel、不是任意 mask 系数。来源的单相位来自完整 \(n\)-torus 的 Haar 推前，不假设源坐标独立。

终态 exit0，1项注册状态、23项 Fraction 常数链、180项浮点 screen，共204项通过。最细结果：

|n|flat \(\mathcal T\)|flat 局部配对|gap \(\mathcal T\)|gap 局部配对|gap 全时间 Cauchy|
|---:|---:|---:|---:|---:|---:|
|8|8.973156|0.000176714|5.542164|0.109008|0.307160|
|32|35.969678|0.000017678|24.662385|0.108986|0.647952|
|128|143.965449|0.000001768|101.163011|0.108976|1.312310|

解析下界负责所有 \(n\ge8\)；不拟合维数阶数。原符号导数身份与正 Bregman/flow 分解经过数值 screen。gap 的 sharp 相位边界是可见误差来源：最细相位占用为0.22998046875，而连续比例是 \(\arccos(3/4)/\pi\)；2048→4096 的 gap test-energy 差分别0.00850、0.03694、0.15077，不能称网格误差已认证很小。最细 \(n=128\) 熵时间积分与解析 drop 差 \(1.755\times10^{-5}\)，winner flow 与相位网格上闭式 flow 差 \(1.424\times10^{-7}\)。脚本 terminal 时间点用左极限的 test 值做 Simpson，真实 strict gate 在终点为0，但单点不改变 Lebesgue 时间积分。

FFT/非线性 log/不连续 test 的采样及 Simpson 都是双精度非区间诊断；Fraction 只认证稿中明确列出的有理常数链。\(L^1(\mathbb R^n)\) 提升是§5解析证明，不是这些 torus 数值的有限盒实验。两个例子都没有原 actual geom 全门资格；它们否定新一般强 test-energy 合同，并保留更准确的配对路线。
