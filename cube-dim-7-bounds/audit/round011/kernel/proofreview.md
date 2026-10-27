# R011 kernel/DD与phase接口独立Sol审计

2026-10-27。新独立clone起点 `75eedb49b0e95dff748081d658d5a9f53a593974`，仅写 `audit/round011/kernel/`。全文读取kernel_proof.md、extensions.md、finite_contract.md、testcases，以及Terra generator/checker/counts/report/receipts和主线程sharp-floor receipt。本审不执行DD、phase、kernel或旧21orbit functional，不重复R012已release的完整generic单点。

**结论：未发现R011数学或程序接口阻断。** positive DD6、confluent derivatives、B^6缩放、真实phase fibres、完整高对称H与generic partial区分、near-distinct余项和sharp ψ floor链均成立。冻结generator与完整phase coverage支持四个真实H_epsilon point intervals；本审不冒称独立重算它们。没有全域g、I或新Cn界，R012另属主线程和Terra。

一项证据层级需保留：checker中重组aggregate的overlap是consistency diagnostic，不能单独证明某个更窄stored interval的外包性；本轮严格外包仍由已冻结generator的完整signed Ball arithmetic、实际source/coverage和原分析承担。本审没有把overlap单独当enclosure proof。

## 1. 结构输入、alpha与DD6正号

R010的135 cells/212 top simplices覆盖ordered cone、内部不交，46 primitive rays来自完整R009输入。ψ在每simplex同一weak-sign chamber上线性，原weight-one64body尺度保持不变。primitive ξ_i=Mw_i、q_i=ψ(ξ_i)>0，detM1，所以

`alpha_S=|detW_S|/Πq_i`，`b_i=x·ξ_i/q_i`。

lambda分母N_i、u-volume的7!不属于primitive alpha。正scale各ray时det和q同步相消。Terra compiler核source复制hash、全64prefix q、212tuple唯一性和source q/alpha，直接读取已审detW；它不是重新证明R010覆盖，完整结构前提来自前轮及其checker。

cone/time积分是∫epsilon∞Π(t−ib_i)^−1dt。distinct knots的partial fraction coefficient i^-6/Π(b_j−b_k)=−1/Π，upper logT系数相消，lower boundary再取负，故恰为**正**DD6[log(epsilon−i·)]。没有额外6!。epsilon>0时right-half-plane log解析，∫e^-epsilonψ/ψ有限（ψ≥c|ξ|），Fourier/time/finite-cone交换绝对合法。

所有signed terms必须保留；单cone/phase real part不保证正。完整group sum才是真实原H_epsilon，没有frequency、shape-expectation或短时截断误差需另补。epsilon本身是定义中的time cutoff，不能把H值当full G值。

## 2. Exact repeats、real/complex等价与B^6

sorted exact knots若子表两端相等，中间全相等，用f^(k)/k!；否则分母是准确非零rational/integer差，不可由rounded balls相减再阈值判零。confluent值

`f^(k)(b)/k! = −i^k/[k(epsilon−ib)^k]`。

minus在i^k外。all-zero top DD给1/(6epsilon^6)，all-one给1/[6(epsilon−i)^6]，都与直接time integral一致。部分重节点遵同一Hermite recurrence，不按perturbation方向保持“等号分支”作假导数。

DD recurrence coefficients对real knots均为real，因此Re DD[f]=DD[0.5log(epsilon²+b²)]。Hermite real data是上式real part，rational epsilon,b时为准确Gaussian-rational real数；b0的odd data0、even2m data(−1)^(m+1)/(2m epsilon^(2m))。Terra real path与complex scalar path是等价算法，不是把complex branch忽略后无根据取实部。

令common B=d_x lcm(q)>0、a_i=B b_i integers、beta=epsilon B。f_epsilon(a/B)=log(beta−ia)−logB，constant在positive-order DD消失，且每一阶node差带1/B，因此

`DD_b^6 f_epsilon = B^6 DD_a^6 log(beta−i·)`。

Hermite derivatives必须对a求、用beta，不在每个cached suborder再乘不一致的B。run.py每top只乘weight*B^6一次，正合公式。allzero test自动B^6/(6beta^6)=1/(6epsilon^6)。所有sorting/equality/差均在common integer knots上完成，精度不改变等节点判断。

cache只在固定physical denominator/epsilon/precision上下文内存sorted multiset；identical top multiset可exact sum alpha。source每phase cache新建，先按top key加Fraction weights，再做signed Ball DD。没有把仅有同一46-value unordered collection但simplex incidence不同的phase合并。

## 3. Phase orbit与fibres：zero、e1、ones、generic

B7作用phase y=Q^Tx，orbit每point fibre为Stab(x)。若z个zero coordinates、各正绝对值重数mλ，则

`c_x=2^z z! Πmλ!`，`|orbit|=2^(7−z)7!/[z!Πmλ!]`。

zero coordinate可任意signflip，nonzero matching sign已唯一；permutation multiplicities正是上述factor。非零x的orbit无negation fixed point，选first nonzero coordinate positive的half orbit，F(−y)=conjF(y)使H为2c_x/(2π)^7乘half real sum。

- e1完整orbit14、half7positive axes、fibre46080。
- ones完整orbit128signs、half64first-sign+、fibre5040。
- fixed generic input七个非零绝对值不同，fibre1、half322560。
- x0 orbit singleton且negation固定；绝不能再乘opposite2。

对generic本轮所有physical coordinates正，所以first sign+确实等于half representative规则。若未来允许负physical input，固定first sign+须相应调整；本轮fixed hex均正，无此错误。

令A_sum=Σ212 alpha。origin使用全group645120一次，给

`H_epsilon(0)=840 A_sum/(π^7 epsilon^6)`，`p_t(0)=5040 A_sum/(π^7 t^7)`。

谱球volume=128A_sum，heat basis512给sanity5040A_sum≤1/512。这里sanity不能代替source decomposition证明。Terra origin direct route与allzero DD route同一归一，opposite_factor记录1正确。

## 4. Near-distinct32-extra-term余项

在δ=max|b_i|<epsilon时，monomial DD恒等式给

`DD6 f=Σ_(m≥0) f^(6+m)(0)/(6+m)! ·h_m(b)`。

可由生成函数Π_i(1−b_i z)^−1或DD递推证明monomial degree6+m的DD是h_m，低于6则0。analytic Taylor在节点disk绝对收敛，包括confluent节点。|h_m|≤binom(m+6,6)δ^m，log系数模1/[(6+m)epsilon^(6+m)]，故m项≤epsilon^-6/6·binom(m+5,5)theta^m。

retain m0..M后，首遗漏m=M+1，tail首项为epsilon^-6/6·binom(M+6,5)theta^(M+1)。后续项ratio≤theta(M+7)/(M+2)，因此余项正是contract所给

`epsilon^-6/6 ·binom(M+6,5)theta^(M+1) /[1−theta(M+7)/(M+2)]`。

M32表示base加32extra共33项；不是只算32项或首遗漏32。Terra准确Fraction形成theta/denominator并验证正性，complete homogeneous recursion按每个node做升序update得到相应generating function。real/complex scalar与该finite expansion加双分量absolute tail overlap是合适校准，近节点仍保持distinct integer equality，没有阈值coalescing。

## 5. 四个真实H区间与generic partial的数据来源

冻结run.py读取全部135cell records和212simplex alpha，无top重复；每phase编译46 exactknots并完整sum212 DD请求（允许相同multisets exact聚合）。e1每epsilon实际7half phases，ones每epsilon64half phases，按2*stabilizer/(2π)^7一次归一。两epsilon各两physical inputs合计142完整symmetry phase evaluations，因而四个输出是完整true H_epsilon点值interval，非sample或orbit partial。

来源报告中的近似值为：e1在epsilon1/64、1/256约3.38730835e−5、1.34139898e−4；ones约1.85656973e−6、3.75439258e−4。正式证据是results.json准确dyadic endpoints和绑定的source，不是这些小数。Terra raw interval保持signed Ball运算，直接测试positive/width，不把negative或width failure剪掉。expected cutoff monotonicity与H_epsilon(x)≤H_epsilon(0)是额外diagnostics；后者也由positive Fourier multiplier的triangle inequality成立。

generic只取π lex前四permutations、每个64sign patterns，phase_id=64perm_rank+sign_rank，恰0..255。完整half orbit322560、fibre1；源码每chunk只乘opposite2，不乘全orbit或extrapolate。其partial sum约5.71776e−6是signed partial，不能叫H、不能与g目标比较或外推未算phase的取消量。planning width target只是后续运行规划门槛，不能保证full run width。

source receipt/results记256bits全部完成、promotions空列表；terminal task IDs只有四symmetry和四generic pilot chunks，counts另把checker的四fixed scalar complex DD和三origin常数工作分开计费。没有fallback/promotion实际证据需补成成功，也没有隐形重复phase。R012 full generic已由主另release，本审没有启动、读取其新结果或重复其任务。

## 6. 锐ψ floor的convex-max证明与出处

P={ψ≤1}是compact convex polytope。任一点是vertices convex combination，convex norm²≤相同weights的vertex norm²平均，所以maximum在至少一个vertex达到。已完整认证的R00321 spectral-vertex orbits按B7覆盖所有vertices，且norm invariant。主线程primary receipt逐项给exact rational squared norms、source catalog SHA256 `856e15aeaa785a72d4bbcd8d590141a8ad7b14fe179fae4671d5ef063f61e333`，max=1/2048。

因此|ξ|/ψ(ξ)≤1/sqrt2048，即 `ψ(ξ)≥32sqrt2 |ξ|`。方向(1,1,0,…,0)中32条body s2+给dot2、32条s2−给0，ψ64、normsqrt2，恰达该常数，所以是sharp floor。该证明不引用未知外部Khintchine结论，也不混用R00946arrangement rays；21spectral orbits才是此convexmax所需的完整vertex来源。

主primary已实际执行exact21-functional，receipt hash与catalog source绑定。本审读其全部rational结果及证明，不再次运行该functional或catalog枚举；理论作者receipt的0运算也不能替主primary执行背书。

对unit direction v，Fourier absolute bound和areaS6=16π3/15给

`||D_v^k H_epsilon||∞≤Γ(6+k)/[120π^4 c^(7+k)epsilon^(6+k)]`，c=32sqrt2。

非unit v乘|v|^k；tube |Imx|<c epsilon由同一absolute domination成立。k1相对c8改善(c/8)^8=2^20。坐标square的spherical average1/7给J2可选bound Γ8/(840π4 c9 epsilon8)。这些是解析常数，无本审数值重组。

## 7. arbitrary-direction稳定子限制与导数公式

multiknot analytic derivative由Cauchy DD表示求导得

`D_v^k f[b0..b6]=k! Σ_(|alpha|=k) Π c_j^alpha_j f[each b_j repeated1+alpha_j]`，

c_j=v·Qξ_j/q_j。multinomial系数与每个denominator derivative的alpha_j!正好抵消，留下k!。公式延伸confluent knots，不能将value中相等节点的分支当perturbed时永远相等。

同一个value phase Q^Tx的fibre中Q^Tv可不同。故e1/ones等symmetry-point arbitrary directional derivative不能简单乘value fibre；须sum/average tensor或使用保留的Q。只有v被physical stabilizer固定（如radial v=x）才可直接使用该fold。generic stabilizer trivial；opposite pairing对导数仍成立，因为整函数关系F(−Q^Tx)=conjF(Q^Tx)可随x求导，但必须保留正确Q/c。作者明确此caveat，没有从四H值推出任意方向导数。

extensions的cone moments公式也吻合：展开coordinate^m并做Gamma monomial/time积分，factor是alpha m! h_m(c)/[(6+m)epsilon^(6+m)]。full group even-coordinate orbit给720系数；m0退回origin840A。完整positive frequency measure有second tensor J2I，故梯度norm≤sqrt(J0J2)，无额外sqrt7；这些moments本轮未数值计算或授权。

## 8. 其他extensions、checker证据层级与停止

epsilon0 generic扩展在所有group/ray节点严格nonzero时允许每个cone time integrand在0附近有界，再由finite sum和H↑G得到real DD log|b|公式。单zero node使该充分termwise论证失效，不证明真实G∞，仍须原S3/cancellation理论。real目标可先取real DD；若要保留faithful complex cone值，须保留principal boundary arguments。该扩展未执行。

compact-time u=1/t变换产生u5 integrand；complete real sum=p_(1/u)/u²非负，而单cone可signed。rho|x|/c≤1/2的disk denominator下界及U·M Taylor余项合同成立，但是未授权alternative quadrature，不是本轮数值结果。

checker重组four symmetry/generic partial时使用overlap，因此是consistency诊断，不单独证明一个任意更窄stored aggregate interval外包。完整生产程序直接sum全部已注册phase balls并归一，其outward arithmetic才是stored true H intervals的严格来源；manifest/prehash、完整coverage和actual receipt必须一起保存。本审没有把checker PASS或作者自审替代这条数学/程序链，也没有另跑whole kernel。

本次只进行源码、数据接口和元数据hash/coverage身份核验，没有新DD/phase/kernel数值。尚无global face g、I或Cn新证书。原full-G finite-source/generalL1/all-common-r interface仍是下游合法比较，四个point H或generic partial不会闭合它。

到此停止：0MATLAB、0commit/push、0R012重复、0新phase或后续任务。主集成人保留最终证明与采用责任。
