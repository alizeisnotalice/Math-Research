# R066正式源的新first-profile局部接入锁

**Local interface ACCEPT。** 源：`/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round066/target_profile/PROOF.md`，实际SHA-256为`4b7b3210f063312527c08adab5b7f81dc4b2c0ecf879e466dc715eb56ea6ebd1`。

实际读取§2 active guards、§3新profile、§4中baseline profile概率/匹配transport、§5每source face均值方差及mixing。以下是局部调用范围的确认，不代审count、radial column、actual-history/all-scale或完整C_n结论。

## 1. Active域和真实旧J输入

源已明确q_max=.917而非旧.909，target上端2751/2000−δ<3/2，h<δ/2使匹配q区间完整落入[2/3,.917]。q_b=(r+h)/s给R=r，整个matching上R∈[r−2h,r]且active R>1。εq_max≤2h/3给correct g=(s−1+ε)q−h>0、g/R≤1/3。inactive fullband R<1点不调用这些需要R≥1的矩。

R064冻结c60引理的真实抽象域没有把.909作为数学上界，故以上guards足以接入J及logJ接口。源的说法不是把新J_phi当作旧结果按名字继承；§3重新写出它的原函数、闭式、derivative与Fubini身份。

## 2. 新闭式和导数的逐项匹配

源H(R,u)=R Q(R,u)²的两尾公式正是冻结新引理；m>g时中央phi=0，确无−m项。它给出的H_u≤6/25、|H_R|≤3/25、H≤R/25合法；两项q导数均≤27/50，再计入ρ'=ρ/q后足以推出8/ρ。我们的更强Q-chain界<1/(2ρ)支持该宽界，并不要求作者改用强常数。

源正确使用a.e.导数与绝对连续性，未把onset kink冒称为classical可微。x=917L/1000≤L、q≥2/3、q宽4h/3使ΣJ_phi shift≤16xh≤16Lh=2^(−76)。源进一步仅需16xh≤x/100000，为h→0的eventual guard，方向正确。

## 3. Target一阶profile概率空间

源的E_W J_phi=E_Zphi是在baseline target W~Uniform[-r,r]的Fubini身份。E_Zphi的精确式与我们的新原函数检查相同。R/q实际为s−h/q≤s，而非严格相等；源明确负h/q项，不影响μ_0上界。

μ_0=1/2−(5/2)log(6/5)<221/5000已独立用even-eight lower-log准确核验。源作者用positive-atanh lower-log证明同一方向也合法；这里不能改用log上界。

源的target free-profile variance≤n/ρ_b²、pinned≤1/ρ_b，比冻结引理中9n/(2500ρ_b²)、3/(100ρ_b)更宽且安全。mean≤μ_0x+1/ρ_b与新阈值的固定正gap给uniform full-target-face o(1)失败概率。这不是actual μ或D-cap后target独立性。source仅在这些概率建立后把target good事件与D相交。

baseline good profile系数221/5000+1/100000再付transport额外1/100000，准确成为2211/50000；源式(14)是对同一target所有连续matching q的确定性上界，不是数值样本的推断。

## 4. 对每source face的mean接口和概率混合

源精确使用free coordinate均值J_phi/J，未在均值这里丢掉capture。固定good target后，target profile condition只约束w，不会对source角坐标另作条件化。每source face可省略其pin对应的一个非负profile，再以全n profile上界支付。因此

mean≤Q_pin+[2211/(50000(1−2^(−20)))]x
<Q_pin+(223/5000)x

与本次新Fraction预算完全一致。pinned source phi为确定Q_pin，未当作另一条free law。

二阶矩只在此variance上界步骤丢capture，使用fullsource 2R Q_pin³/(3ρ)并除以J下界，固定face内独立相加，得var<x/150。源明确每face均值可以不同；χ−223/5000−1/5000=1/5000，pin/x≤1/5000后每face Chebyshev≤5000²/(150x)。最后混合这些统一概率上界，未把mixture variance等同于条件variance和，也未在D-cap/soft b条件化后继续假设独立。

## 5. 冻结依赖和本次身份核验

实际核对以下三个数学文件保持首次发布字节，不作任何增量编辑：

* `LEMMA.md`：`67bdf01d9aa65c70bf058df08b527b904e954d19f3327ecbb43e5b5c316b5030`。
* `verify.py`：`28629535461b3460e09b53f458d1f0e70ef1f813b7d95e3c25be3a6082476abb`。
* `verification.json`：`01077db2d4b54aff9cbc1296df374b19af6b1d5e3a344a9ed3c47797760f0f98`。

本文件是纯源码读取、数学接口对照及身份锁，不是新的计算run。没有重跑作者16个gates、旧J/baseline或任何生产grid，未用MATLAB、n=7、下界、Git，也未修改原source或冻结c60文件。
