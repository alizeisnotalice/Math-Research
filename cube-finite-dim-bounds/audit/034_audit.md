# 034 整face短time Lq与11O均值/凸弦接口审计

2026-10-25。结论通过：主 `proofs/round011_face_time_budget.md` 的整face短time Lq证明，以及锁定11O的degree32连续插值、准确mean、stable chord、面积/分part接口均无关键数学缺口。已提前向主反馈可继续全I计算。

本批不是完整8256个patch节点/记录的验收。没有读取live part或节点cache，没有调用Green或新node evaluator。当前严格g≥131/25000已由033单独认证；I仍未在034接受，未宣称新U3/U4。

## 锁源与实际有限核验

16项源锁在 `034_source_snapshot.json`。主短time proof SHA `336244fcffc2c69b8b3af4c213cb1b7a8fc249ec63f2160296b776a971f617cd`，Python `8840d561b60cc9cb4579b5791708fea1fd61b3c70a23a72782eb7544a3ce0096`，结果 `beab285a80fff53f00f4fb4b8d48442fbc0312736b81a0a9c26618d7227e5d4b`。11O checker为 `4f89a377fbdf263061c2d6525d2101a773f7a6ba4127bb367bd3b7328ed57eb5`，report为 `9a837612b6e21f0de4dd0f73b7b82ee130975e541055c9b5cd250c0df8de11c3`。

独立Python第一次通过，内部2.825秒，无实现修补。重建n3/n4全部active方向补基、四个主ε预算与当前ε=9/5000预算；33个完整Bernstein多项式/mean恒等式、1089 DCT entries、8256面积1及12part准确分区均过。四个已锁pilot的L/U/m标量以512bit独立重算chord，其上端不超过作者384bit准确contribution端点。这不把四pilot的节点/cache视作已独立全算，也不以四pilot代替全I。

共享SymPy1.14.0、Fraction、python-flint/Arb及032/033分析基础。256/384/512bit分别仅用于有限常数、矩阵及标量重算，没有改正在运行批的384bit精度。未执行作者checker或point函数。

## 一大跳事件必须在空间density层面处理

令active集合为vi1≠0的方向、k为其个数，A=Σactive ci；这里vi1∈{−1,1}。若第一坐标=1，至少一个active参数|Zi|≥1/k。选该列及n−1条权1的body补列为basis，det绝对di。受限Cauchy一维密度≤ci k²t/π，补列各≤1/(πt)，再除di；其余变量作概率卷积，不增此sup。事件重叠只作上界求和，得

\[
p_t(1,y)\le C_f t^{2-n},\qquad
C_f=k^2\pi^{-n}\sum_{active}c_i/d_i.
\]

固定face是零概率slice，不能直接把其“事件概率”当密度。严谨升级是先在x1>1−η邻域取阈值(1−η)/k，得到a.e.密度上界C_f(1−η)⁻²t^(2−n)，再由真实pt连续延拓至该face点并令η↓0。权已在Cauchy scale ci t中，补基用未加权方向列，不能再在det重复乘ci。

精确枚举核：n3 τ1/3有4条body、4条active face，k8、A16/3，每条可补det4，故Cf=(256/3)/π³。n4 purebody有8条active、A8、det16，Cf=32/π⁴。主最初把face最大det手算为2的较松猜测已经在锁稿前更正，本审核的是准确最终4，没有沿用旧320/3。

## 真plane marginal与Minkowski指数

第一坐标是独立Cauchy和，scale为At，其canonical marginal密度是At/[π(1+A²t²)]。这可以直接从basis-product表示作正Tonelli：用一条|vi1|=1参数解normal坐标，其余n−1个参数作free平面坐标，Jacobian与basis det消去，得到一维Cauchy卷积。因此该恒等式在固定slice逐点成立，不需假设任意L1密度的slice代表自动连续。即便仅用连续pt和Fatou得到该slice的marginal上界，也已足够本不等式。

用**同一真实pt**的sup与marginal，p^q≤sup(p)^(q−1)p，q=n/(n−1)，再舍去1+A²t²≥1，得

\[
\|p_t(1,\cdot)\|_q\le
C_f^{1/n}(A/\pi)^{(n-1)/n}t^{1/n}.
\]

时间幂准确为[(2−n)(q−1)+1]/q=1/n；不是从各坐标独立平均或不同过程相乘。先有限time/空间截断，再Tonelli/Minkowski及单调极限给wholeplane Tε=∫0εptdt的norm≤n/(n+1) times该常数 times ε^(1+1/n)。该证明对任意ε>0成立，不限于主程序原展示的1/1000、1/2000，故11O的9/5000合法。

真实过程在各free坐标分别反射不变，Tε也even。orthant q次积分是wholeplane的2^−(n−1)倍，norm因子是2^−((n−1)²/n)，n3为2^−4/3、n4为2^−9/4。单位positive face是orthant子集，因此主Eε常数正确。它包括整个plane的奇线邻域及空间无穷尾，不删帽或截远处。

Tε先按非负time积分定义，G=Hε+Tε由Tonelli成立；不是作∞−∞。若I_H≤B，则||G||q≤B^(1/q)+Eε，故I≤(B^(1/q)+Eε)^q。所有量正时，目标的允许B为(I_target^(1/q)−Eε)^q；前括号须严格正。使用上界时取E包围上端，使用允许阈值时取B包围下端。不得将确定表达式的下端冒作未知实际T norm下界。

独立四主预算区间均严格落在主192bit端点之间。当前ε=9/5000、I_target=29/50000：

\[
E_\epsilon\approx0.0001299691295463827922367652480088,
\quad B_{allowed}\approx0.0005638179328648305750722316975048.
\]

完整准确有理端点见 `034_results.json`；这里显示小数不作为严格阈值。

## 11O全域插值、mean与凸弦

Hε正则化real Hermite DD与032/033同一已审time公式，ε>0时零/重节点都正常，源signed项完整保留。新参数h=1/256、rho5/2、degree32给ellipse高度21/5120、ε−s/3=277/640000>0。由已审ψ≥3(|ξy|+|ξz|)，Fourier tube bound给M；Cheb系数与roots alias、Lebesgue≤65、tensor误差如033逐项替换，得到Epoly=176 M(2/5)^32<7.706293505929e−8。real Cauchy radiusR=ε给导数和2Dh2⁻128位移误差<3.938477182864e−36。代码把该误差在node ball中先传播，未在高放大变换后少付。Arb radius为向上mag舍入，实际node radius可略大于这个解析表达式，不应错误要求两者的展示数逐位相等。

degree32精确Bernstein转换的33个完整多项式恒等式由独立SymPy重建。准确normalized mean为μ0=1、奇数μk=0、偶数μk=−1/(k²−1)，它等于∫0¹T_k(2u−1)du，也等于Bernstein列系数和/33。独立全部33个等式及1089 DCT inverse entries通过。mean(P)=μ A V Aᵀ μᵀ对应**平均值**，面积尚未乘入，不能再漏/多乘1/2。

Bernstein上下端及Epoly给0≤L≤H≤U，U>0；H的mean≤mean(P)_upper+Epoly且≤U，所以m取这两个上界的min合法。程序要求m≥L，不能向上修补一个不合格mean。对凸f(h)=h^(3/2)，整个[L,U]的上弦给

\[
\operatorname{mean}(H^{3/2})
\le L^{3/2}+S(m-L),\quad
S=\frac{U+\sqrt{UL}+L}{\sqrt U+\sqrt L}.
\]

设a=√U,b=√L，则S=(a³−b³)/(a²−b²)，证明准确。L=0时S=√U；U=L>0时S=3√L/2且m=L，结果L^(3/2)正确；源码不除U−L。U=0本模型不发生，源码assert排除它，没有将失败静默当零贡献。S≥0使mean上界m可安全代入。四pilot仅对给定L/U/m和面积的标量运算被本批高精度复核；其point/node范围的全批接受仍是另一个步骤。

## 面积、分part及最终验收边界

全部indices 0≤j≤i<128共8256。j<i方块与转置方块积分相等，权2；i=j的整方块已经含diagonal两侧，权1。权和16384，乘面积1/128²恰1。闭边界重叠为Leb零测，不遗漏三角半或重复付diagonal。不能把033的extra closed覆盖论证误当成这里的积分乘重。

serial index i(i+1)/2+j从0到8255，mod12各688个、两两不交且并集完整。代码aggregate要求全部8256唯一patch并准确Fraction求和contribution，没有先舍入double。完整接收还应核冻结checker/config、各part正确index集、缓存来源/执行记录和每条range/mean/chord的准确端点一致；作者selfcheck提供这些必要检查，不等同本批独立重算所有节点。

本034没有打开live part或cache，未验收全8256实际数值证书。数学上，只要全部patch的node/系数包围正确，接口给严格Bsum，再以wholeface Etime闭合I。源代码/四pilot/主正在运行的part不能未经最终覆盖和聚合核对就声称I目标通过。

033已独立给g≥131/25000；若以后I≤29/50000获接受，则64g³−27I²=979537/7812500000000>0，可按031得到<8的成本。这个准确条件比较不是本批已经更新U3。

## 实际工具与终止

主导入后最小有限核验（stdout-only；可选--save）：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-finite-dim-bounds/audit/034_check.py
```

审计clone加 `--project /Users/zhengzhihao/Desktop/T/Math-Research-finite-dim/cube-finite-dim-bounds`。16来源开头/收尾锁定，未改任何正在运行的参数。MATLAB仅主串行核指数、orthant、stable chord、mean、面积和当前time预算，实际见 `034_matlab_receipt.json`；不宣称它重算全I。本理论有限批锁manifest后停止，后续n4表示/g/I按新035登记接续。


主串行MATLAB已实际exit0（session79257），同SHA dbe4c5f5e43694014b77cde5f1eeee97279eaeadb199c54268d229fb8cf39b0b。审计已读stdout，指数/orthant/chord/means/面积assert通过，Etime=.00012996912954638281、allowedIH=.00056381793286483062；log SHA5a414efd58f172206dd9d23166cf8382dd7715d77044edefe690d3c565b4f681。浮点展示不冒称位于256bit窄包围内。收尾16来源未变，034理论有限批已停止。
