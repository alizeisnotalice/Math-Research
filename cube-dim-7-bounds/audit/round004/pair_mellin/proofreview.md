# R004 pair-Mellin 独立 Sol 理论审计

2026-10-27。独立 clone 起点 `b2139fc1fe70fd55e7d409031256195c7d10e768`，仅写 `audit/round004/pair_mellin/`。本审不调用 MATLAB，不重复任何数值 mesh/panel 系数，不提交推送，不改变其他工作者或 R001–R003 工件。主集成人仍保留最终证明和采用结果的责任。

**结论：接受 pair-Mellin 解析推导、严格改善及 Terra 端点到完整 I 的接口。** Mellin/Fourier 半线系数、绝对 Fubini、复解析域、Cauchy 包络、Taylor 总误差和无限尾均有正确的自足证明。给定 Terra 交付的512面板系数区间，精确账本汇总和本审768-bit标量重组证明

`I < 8179950070308/10^18 = 0.000008179950070308`，

严格优于 R002 的完整解析预算表达式。此为完整 I 上预算，没有认证全域 g，没有新 C7 上界。数值面板由 Terra 负责；本审没有独立重算面板，因此不把本审称为第二套 quadrature 实现。

两项不影响主证明的记录偏差必须保留：实际 pilot 只核 F 的零阶系数与直接点计算以及点态 Gamma recurrence，没有执行建议的 f(0) 或 series-level recurrence；Terra 的 postrun compact-proof hash 是较早 b512…，现含独立 Laplace/Plancherel §8 的最终 pair_proof 为 f12c…，而预运行完整 theory.md 的1c6f…保持一致并已含全部积分/误差公式。这不应表述为 compact proof 在整个运行期间哈希未变。

## 1. 实际源与审计边界

完整读取 Sol `strategy/theory.md`、`pair_proof.md`、`pair_registration.md`，及 Terra `pair_mellin/registration.md`、run.py、assemble.py、check.py、certificate、report、sources、receipt、manifest 和前后源哈希记录。准确路径/哈希全列于 results.json，最终 compact proof 锁定 SHA256 `f12ce577b2fac6f9c845a9664895bd9ba39136c86c8eaea222795888b98c97f6`。原 proof 的与本轮 pair 无关的诊断/条件化建议不属于本审接受范围；只接受下面的 pair 路线和完整 I 链。

本审 minchecker 未导入或执行作者代码；仅读取512条面板终态记录并按 Fraction 精确汇总二进制上下端。没有新增面板、节点、随机数、旧 baseline science、优化或指数变体。标准 Gamma/Fourier/Plancherel、普通AM–GM、Tonelli/MCT、Cauchy估计和球算术被当作明确工具，不由有限 PASS 替代。运行中还核对主源前后相同、Terra manifest全项相同，以及本 clone 201个旧封存文件前后相同。

## 2. 八独立块、负矩与配对系数

R002 的64方向被八个不交7×8 Hadamard tight frame 分割；每组 Gram8I、全部七阶minor模512。原 iid S_i=1/(4Y_i)、Y_i~Gamma(1/2,rate1)，共同 Cauchy 时间未变。每块 Σ_g SPD，T_g=detΣ_g^(1/7)>0；不交变量保证八个T_g独立同分布，非随机frame对称保持determinant law。所有这些均是既有被审结构，没有另换生成元或高维尺度。

由 `detΣ_g=16(ΣY_i)/(ΠY_i)`，令 R=ΣY_i、p_i=Y_i/R，有 R~Gamma4 与 Dirichlet(1/2)^8独立，`T=16^(1/7)R^-1(Πp_i)^(-1/7)`。对 Re z>0，绝对 Gamma/beta 积分给

`m(z)=E T^-z=2^(-4z/7) Γ(4+z)/Γ(4+8z/7) [Γ(1/2+z/7)/Γ(1/2)]^8`。

Re z>0 下所有绝对积分正有限，且在紧子条带可由邻近实指数矩控制导数/差商，故 m holomorphic，`|m(v+it)|≤m(v)`。Hölder 给 log m(v) 凸。T非退化连续：固定内部p后，R有连续Gamma密度，T=正系数/R，所以 T1=T2 的概率为0。

SPD determinant Minkowski 给 `M≤E(Σ_1^8T_g)^(-1/2)`。把八个T配成四个独立同分布 U_j=T_(2j−1)+T_(2j)，四项AM–GM使

`M≤4^(-1/2) E Π_j U_j^(-1/8)=(1/2)q_(1/8)^4`。

因子1/2与pair负指数1/8正确；a=s/2=1/16只属于后续 Mellin分裂。pair AM–GM本身又给 `(T1+T2)^(-s)≤2^-s(T1T2)^(-s/2)`，以概率1严格，正期望差有限而非零，故

`q_s<2^-s m(s/2)^2`，

`M_pair=q_(1/8)^4/2<2^(-3/2)m(1/16)^8=M_H(R002)`。

这是上界表达式的解析严格改善；不是 actual M 的下界，也没有量化其与 actual M 的间隙。

## 3. Mellin-Barnes 系数、Tonelli/Fubini 和独立检验

a=s/2>0。beta积分代 y=e^r，得 `[2cosh(r/2)]^-s` 的 Fourier变换是 `Γ(a+it)Γ(a−it)/Γ(s)`（偶性使 Fourier正负号选择不影响值）。因此 Fourier inversion 的全线系数是1/(2πΓs)，而半线因偶性为1/(πΓs)。代

`(x+y)^-s=(xy)^-a[2cosh((logx−logy)/2)]^-s`

并对两个独立T取期望，得到

`q_s=[πΓ(s)]^-1∫_0∞ Γ(a+it)Γ(a−it)m(a+it)m(a−it)dt`。

该积分不是形式 contour操作。取绝对值时 `E(T1T2)^-a=m(a)^2`；Gamma乘积可用正实Gamma积分与下述整数recurrence给可积尾，且t在紧区间无极点，所以绝对Fubini合法。配对后的 integrand 是 `|Γ(a+it)m(a+it)|²≥0`，不会遗漏复部或符号。

独立系数核对可完全绕开 beta inversion：令 L(u)=Ee^-uT，Laplace恒等式与Tonelli给 `q_s=Γ(s)^-1∫u^(s−1)L(u)^2du`。g(r)=e^(ar)L(e^r) 的L1范数为Γ(a)m(a)，L2范数平方为Γ(s)q_s，后者由pair AM–GM有限。其 ordinary Fourier变换是Γ(a−it)m(a−it)，合法绝对Fubini；Plancherel便给同一全线1/(2πΓs)及半线1/(πΓs)。Sol最终§8的新增证明正是这条独立核对，常数没有误差。

## 4. Holomorphic strip 与全局 complex Cauchy 包络

解析延拓必须在积分变量w进行：

`f(w)=Γ(a+iw)Γ(a−iw)m(a+iw)m(a−iw)`。

|Imw|<a时，两参数实部分别a−Imw与a+Imw，均正，Gamma无极点且m的积分表示合法，所以f holomorphic。复域不存在将conjugation误当holomorphic的步骤；实轴上conjugacy仅用于生产实系数。

置 h(v)=Γ(v)m(v)，两个因子均正log-convex，因此h也log-convex。Gamma的定义积分和m的moment积分给 `|f(t+iy)|≤h(a−y)h(a+y)`。对 |y|≤r，以θ=(y+r)/(2r)分别对h(a+y)及h(a−y)应用Hölder/凸性，两条指数互补，相乘恰得 `h(a−r)h(a+r)`。不需要Gamma在某区间单调，且此界与t无关。

r=3/64<a=1/16，左右端为1/64和7/64。每个以实panel中心为圆心的闭Cauchy圆都严格在strip内部，故统一 `M*=Γ(1/64)Γ(7/64)m(1/64)m(7/64)` 确可用于全部512圆，包含最靠近t=0圆的负实部区域。

## 5. Taylor 余项和无限 Gamma 尾

512 panels覆盖[0,32]，每宽1/16、h=1/32、中心(2j+1)/32。Cauchy估计 `|b_k|≤M*/r^k`，degree128的第一遗漏项k129，|u|≤h时剩余≤`M*(2/3)^129/(1−2/3)`。每panel积分乘2h，512项总长度正好32，因此

`E_raw=32 M*(2/3)^129/(1−2/3)`。

没有多乘512或少乘2，不能把128误作遗漏指数。多项式偶数系数积分 `Σ_(k=0)^64 2h^(2k+1)b_(2k)/(2k+1)` 正确，奇数项在对称区间为0。所有负系数和interval半径均需保留，Terra代码确实这样相加。

另对任意整数N≥1，Gamma recurrence与 `|Γ(a+N+it)|≤Γ(a+N)` 给 `|Γ(a+it)|≤Γ(a+N)|t|^-N`。乘 `|m(a+it)|²≤m(a)²` 后，N16、T32的整个无限尾积分≤`m(a)²Γ(a+16)²/(31·32^31)`，再除πΓs。此为对全部t≥32的上界，无非严格asymptotic或隐形截尾。finite-segment误差先除normalization，正tail只加上端，下端无需扣tail。Sol和Terra的单位安排一致。

## 6. Terra 系数生成实现与接受条件

run.py中全局ctx.prec512、ctx.cap129，z-series常数为准确a+i·center，一阶系数为i。`F=Γ(z)m(z)` 使用gamma、reciprocal gamma、exp与准确有理参数；2的幂按real log2计算。固定strip内所有Gamma参数实部正；rgamma entire，未跨未处理branch。已读本机python-flint的acb_series文档确认ctx.cap控制truncation阶，gamma/rgamma是所用series接口；没有运行新scientific series。

令F_j(u)=F(a+i(t_j+u))。因为Gamma和m满足实参数conjugation，其解析配对F(a−i(t_j+u))的Taylor系数正是F_j系数的conjugate。因此coefficientwise conjugate series相乘给f在该中心的正确系数；它不等于对复u直接取非解析absolute-square。检查imaginary coefficients含0正确，取其real enclosure并积分安全。FLINT的O(u^129)只保证有限系数；函数Taylor余项由独立Cauchy界支付，不能把形式O项当数值误差。

所有512终态独立编号完整、正dyadic积分端点、原source/code哈希均通过读取核对。面板系数的正确性在数值层面依赖Terra generator和所用certified ball-series库；本审只检查其数学规格及源码，不声称第二次重算验证。生产report已明确该责任分界。

建议pilot的f(0)和series-level recurrence在run.py没有实现；实际j<4检查的是零阶F系数与direct point overlap及点Gamma recurrence。这是建议诊断未全部执行，需要如实记录，未破坏上述series算术和解析remainder的充分证明。若主要求严格逐项履行这项工程合同，可另登记诊断；本审不为此擅自重跑任何pilot。

## 7. 本审实际三组固定标量、修复记录与端点

已先登记，独立读取ledger端点并Fraction准确求L、U；未用certificate中已算好的polynomial sum当输入。768bit三组标量为：M*/E_raw/tail；由准确L,U到q/M/I上下端；固定R002表达式比较。m用log-Gamma真实数参数计算，I系数由原始C^(7/6)H/896路径算，未执行作者checker。

第一次自身检查直接用Arb union(qlo,qhi)重建球，并要求该新球严格包含在Terra旧球内，断言失败。Arb有限精度radius重新外包不保证两个有效球嵌套；这不是公式或证书失效。改为保留准确qlo、qhi并分别通过单调正幂传播后，第二次全部通过。此实现修复未增加数学范围、mesh或科学采样，收据保留一次失败断言/修复。成功计时0.06721304100938141秒，Python3.15.0、python-flint0.9.0，共享Arb库，非独立算术内核。

M*/raw error/tail/q/M/I的独立768区间均被Terra certificate确切dyadic区间包含。完整I新上端严格低于RI=8179950070308/10^18，RI严格低于本审R002_I_H768lower；同样M_pair768upper<R002_M_H768lower。完整有理上下端和RI减本审上端的正差保存在results.json。这里的lower endpoint是预算表达式的lower，不能充当actual I或M的lower。

## 8. 溯源偏差、原一般L1/全r接口及结论限制

Terra预运行锁完整theory.md为1c6f392f…；该文件今仍相同，§§4–6含全部pair公式、tail和Cauchy参数。Terra postrun compact proof记录b512…，现final pair_proof=f12c…并含独立Laplace/Plancherel§8，registration参数ef054…始终一致。本审未取得旧compact字节，因此没有逐字确认两版本是否仅相差该附录。本审同时记录两hash，接受最终证明但不追认final compact在运行前存在。参数与数值算式未因此改变；主集成应保留这条历史说明。

从M上端到I使用原共同时间Cauchy的全boundary Gaussian/Minkowski链：先对完整boundary做shape积分与Minkowski，再对期望后的原law用14face×64orthant=896，不对单shape正face强加对称。积分支付全部时间和奇集邻域，无需另加R001短时预算。此前接受的正部势全空间source、同一障碍D对全部common r>0、一般L1单调逼近链仍适用；本审不从Mellin策略单独推这条已封存比较。

因此可以采用新的完整I上预算。要得到C7上界仍需另外证明 `0<g0≤essinf_Q G`，再用 `C7≤7^(7/6)I0/g0^(7/6)`；这一步尚未闭合。没有新g证书，没有新C7数值上界，也没有model exclusion结论。此审在固定范围结束，不扩展到其他策略诊断、MC、merge quadruple积分或额外数值网格。

复现本审接口检查：

```
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-7-bounds/audit/round004/pair_mellin/minchecker.py
```

该命令只重读交付ledger并重组标量，不重算面板系数。
