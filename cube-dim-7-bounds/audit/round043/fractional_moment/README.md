# R043 实p矩、熵与uniform-cut模板障碍独立审计

**ACCEPT最终修订稿** `upper/round043/fractional_moment_cut/PROOF.md`，实际SHA-256：`d8fe1bb71c9072590f65af1c78db2500e8fd5e832a815d1f96dbd535a595be5e`。

初审源身份为`84d202360a932dd952d22389bc9df8382eec28fb1d4fa656a38b4c8197756ae8`。初审发现的envelope/真实column能量桥接问题已修正：最终§4只对冻结slice的乘积envelope Y自身的倾斜能量给诊断；最终§5只限制具体moment上包支付模板，或显式以Y-energy为保留目标的模板。**不接受从C_A≤Y反推真实kernel C_A的能量截断必要性。** 新版还补齐bounded-A时固定retention η强制H≥η。

实际read-only diff及独立脚本确认两版本全部21个displayed math blocks相同，改变的是上述scope与bounded-A解释。作者十项检查未重跑，本审计未改作者源文件。

## 1. 实参数矩与固定p展开

条件g≤min(r/2,ρ/2)确保X-law的连续尾为u∈(0,g)上密度1/ρ，且ℓ∈[3/4,1]。其余质量为ℓ=1处的atom。对任何实p>1积分(1−u/(2r))^p，得到式(1)，没有整数p限制或从有限数值点外推。

独立符号核验改用a=1−g/(2r)、k=2r/ρ。矩为

M_p=1−k(1−a)+k(1−a^(p+1))/(p+1)。

对a求导得k(1−a^p)，且a=1时M_p=1，确认整个实p公式。均值与方差分别为

m=1−k(1−a)²/2，Varℓ=k(1−a)³/3−k²(1−a)⁴/4，

等价于原文式(2)。由g/ρ≤1/2，实际可得Varℓ≥(5/96)g³/(r²ρ)，因此作者使用的1/24更保守且正确。

令s=1−a，独立符号展开精确给

log(M_p/m^p)=k p(p−1)s³/6−[k p(p−1)(p−2)/24+k²p(p−1)/8]s⁴+更高项。

三次项就是p(p−1)g³/(24r²ρ)，四次项分别落入作者所列g⁴/(r³ρ)、g⁴/(r²ρ²)误差尺度。固定p下ℓ始终远离0，Taylor常数合法；g/r与g/ρ同时趋零时误差相对leading项消失。该展开本身不用于p↓1的uniform推论。

## 2. 熵下界及所有p>1

U=ℓ/m满足EU=1且U∈[3/4,4/3]。u logu的二阶导≥3/4，使其关于1的强凸下界在取期望后给Entropy≥(3/8)VarU；又VarU=Varℓ/m²≥Varℓ。结合保守variance bound便得g³/(64r²ρ)。均值实际更接近1，但不需要优化常数。

对M_p在p=1求导等于Eℓlogℓ；独立SymPy检查其恰为k[(a²−1)/4−a²log(a)/2]。除以m再减log m即原文精确entropy式(5)，未遗漏atom（atom在ℓ=1时对ℓlogℓ贡献0）。

ψ(p)=log EU^p是凸函数，ψ(1)=0、ψ'(1)=EUlogU；因U正且有上下界，求导及凸性对所有实p合法。切线下界给ψ(p)/(p−1)≥Entropy。固定face只有N=n−1个独立free coordinates，乘积的log矩加N倍；固定face坐标及任何正的deterministic prefactor恰在均值归一化中抵消。式(6)因而覆盖任意接近1或更大的p，不除未知O_p余项。

## 3. Markov支付成本与倾斜方差

式(7)是选择某个p>1并使用该特定moment上包付尾的证书。要求此上包≤θM，就必有(p−1)logH≥log EwidehatZ^p+log(1/θ)，故式(8)正确且uniform于p。即使Y确定，第二项也可能是证书artifact；原文对此有正确区分。

dQ=widehatZ dP确是概率，且独立乘积density使每个坐标仍独立、单坐标density为U。log widehatZ的Q均值为N·Entropy。其方差等于ΣVar_Q(logℓ)，确定的−logm没有variance贡献。

由−logℓ≤(4/3)(1−ℓ)、U≤4/3，

Var_Q(logℓ)≤(64/27)E(1−ℓ)²=(16/81)g³/(r²ρ)≤g³/(r²ρ)。

故总variance≤A。阈值A/128距离均值至少A/128，Chebyshev给A/(A/128)²=16384/A。A≥65536时≤1/4，原常数全部正确。对一般固定Y-energy retention η>0，只需A大于依赖η的固定常数便能强制H>exp(A/128)；小A时该概率界可能平凡，不应拿它当有效尾估计。

## 4. 必须保留的Y与C_A区别

R042证明的是实际列C_A≤Y，而Y为此处deterministic prefactor乘独立坐标乘积。精确恒等式

E[Y1_(Y≤HM)]/EY=Q{widehatZ≤H}

只表示**Y-envelope能量**。一般非负C_A≤Y可以主要集中在低Y事件上；Y损失大部分能量时C_A仍可能保留全部。因此仅上比较不能给反向能量保留结论，也不能把(10)转成真实C_A的必要cap。cap条件、半径混合、actual μ或coordinate-dependent prefactor还会改变归一化law。

新版已经明确：式(10)是冻结full-face乘积law的诊断，不是原kernel或C_n障碍；若将来需要转移，必须另证reverse energy comparison或重新建立对应law。此scope修复后没有剩余阻断缺口。

## 5. 统一列截断模板的order限制

在固定比较常数、充分大n、g~d、r~1、ρ~n/L（本上下文L=log n）下，N=n−1≥n/2给A≥a_0Ld³。

若证书使用式(7)付尾，式(8)直接强制H≥exp(A/64)。若独立指定其目标是保留固定η比例的Y-energy，则在大A由式(10)强制H>exp(A/128)。在bounded-A部分，A≥a_0Ld³已给d=O(L^(−1/3))；再由E[widehatZ1_(widehatZ≤H)]≤H知H≥η，补齐reported deficit≤Cnd/H的同阶上界。

因此在明确的reported deficit≤Cnd/H架构中，最大可报告量受d exp(−a_1Ld³)控制。其导数符号为1−3a_1Ld³，唯一正极大点d=(3a_1L)^(−1/3)，最大值(3e a_1L)^(−1/3)。式(11)–(12)及O(nL^(−1/3))模板障碍正确，允许p依赖n也不能绕过all-p entropy下界。

这不是实际missing mass Δ的上界，更不是C_n下界或所有新证明的障碍。变换kernel、row-dependent付款、radius-dependent cap/retention、真实共同mass、actual-history几何及对C_A更精准的尾估计均未被排除。d必须处于前述几何law的合法范围；无须把形式上的d>0优化理解成所有d均可几何实现。

## 6. 本机新核验及复现

实际使用本机SymPy 1.14.0完成上述实p导数/边界、精确mean/variance、entropy的p导数、三次和四次log系数、16/81 variance系数、16384 Chebyshev常数及单变量极值导数检查。没有运行作者的整数p=1,…,8及两项fractional-case脚本。

另外使用512bit Arb在新点(r,ρ,g)=(5,13,14885/16384)检查精确entropy式与下界的strict margin。该点满足g≤min(r/2,ρ/2)，entropy约为9.3657074871×10^(−5)，下界为10147595705/281474976710656，Arb定向包络严格证明此点margin为正。它只是新实现核验；一般实p与所有n结论来自解析证明。

```sh
PYTHONDONTWRITEBYTECODE=1 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-dim-7-bounds/audit/round043/fractional_moment/verify.py cube-dim-7-bounds/upper/round043/fractional_moment_cut/PROOF.md
```

可选第二路径传入原84d版本，则还核对全部display math blocks保持一致。原stage文件不属于归档必要源；`verification.json`已记录原、新身份及21个block比较结果。

本审计只写本目录三个文件；未改其他代理文件，未用共享MATLAB、n=7、生产扫描、网络或Git修改。
