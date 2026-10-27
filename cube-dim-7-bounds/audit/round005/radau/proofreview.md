# R005 Radau 完整积分预算独立 Sol 审计

2026-10-27。独立 clone 起点 `73646ed3d9818eb1d87a157aac09b8a709718ad7`，仅写 `audit/round005/radau/`。全读 Sol `merge4_proof.md`、registration 及 Terra run/assemble/check/certificate/receipts/manifest。Sol proof SHA256=`f456482119e204b06226875f4a9ef3a3ba49e2c690e45cbc18f2997ed45a287e`，Terra certificate SHA256=`2a1ab8ece05bbd95cbede8b4df607a5ae49ec4538e20e2ecf9e3642bffec6bf3`。

**结论：接受 Radau 全域解析上界及交付数值端点的标量接口。** 独立核 complete monotonicity、Hermite 全半线余项、μ/ν端点方向、R004完整q逐字复用、三组s各自的Cauchy及无限尾参数。没有重复1792个科学面板；给定交付面板系数区间，本审准确账本及768-bit标量重组证明

`I < 7730422817512/10^18 = 0.000007730422817512`，

严格优于R004 pair预算。没有全域g证书，因此没有新C7上界。面板数值生产由Terra负责，主集成保留最终采用及证明责任；本审不冒称独立第二套mesh验证。

## 1. 模型、四块恒等式与倾斜矩

八个不交Hadamard块及原weight-one body64共同时间Cauchy保持不变。T=detΣ_block^(1/7)>0的iid性与m(z)来自此前已审R002，m对Re z>0解析且负实矩有限。U=T1+T2，q(r)=EU^-r；q0=1严格来自probability normalization，不能把r=0塞入有Γ(r)极点的Mellin公式。

把八T分为两个独立四项和V1,V2，determinant Minkowski和两项AM–GM给 `M≤2^-1/2 [E V^-1/4]^2`。取s=1/4、a=s/2=1/8，U的独立a-tilted law密度为u^-a/q(a)。对两独立tilted U，Δ=logU1−logU2、τ=1/16，有

`q_four(s)=2^-s q(a)^2 E_a cosh(Δ/2)^-s`。

令W=cosh(τΔ)−1≥0，1/(2τ)=8，所以cosh(Δ/2)=T8(1+W)，g(W)=T8(1+W)^-1/4。这一变换没有漏cosh或2的尺度。

由于倾斜后U1,U2仍独立，`E_a e^(kτΔ)=q(a−kτ)q(a+kτ)/q(a)^2`。k1、k2分别给R1=q(1/16)q(3/16)/q(1/8)^2和R2=q0 q(1/4)/q(1/8)^2。由cosh²x=(cosh2x+1)/2，得到

`μ=R1−1`，`ν=R2/2−2R1+3/2`。

所需最小q指数为0，其余均正；所有矩有限。U连续非退化，正tilt保留非退化，故μ>0且ν≥μ²>0。有限区间仍需单独分离该正性，作者没有拿定性证明替代数值验收。

## 2. T8根与complete monotonicity

T8(cosθ)=cos8θ，有八个不同实根cos((2j−1)π/16)，leading coefficient2^7。因此

`T8(1+w)=2^7 ∏_(j=1)^8(w+c_j)`，`c_j=1−cos((2j−1)π/16)>0`。

本审另用准确多项式递推核T2、T4、T8与普通Chebyshev一致、T8(1)=1、leading128，没有用浮点根位置作证明。每个(w+c)^-1/4的k阶导为(-1)^k乘严格正rising-factorial及正幂。乘积的Leibniz各项同号且至少一项严格，故g在整个w≥0上complete monotone，特别g>0、g'<0、g''>0、g'''<0，g0=1。

这一步需所有c_j严格正，恰由根均小于1保证；不存在把T8负区间或复根错误用于非负半线的隐患。

## 3. 全半线Hermite majorant与Radau因子

d=ν/μ>0。唯一quadratic P匹配P0=g0、P(d)=g(d)、P'(d)=g'(d)。对任意w≥0，重复节点Hermite余项为

`g(w)−P(w)=g'''(ξ) w(w−d)^2/6≤0`，

ξ位于0,d,w张成的紧区间；在节点取等。因g'''在全半线同号，此结论对w>d也成立，绝不是只在[0,d]的局部插值估计。

W有有限二阶矩，P的期望合法。两点probability在0放1−λ，在d放λ，λ=μ²/ν∈(0,1]，恰匹配一、二阶矩。所以

`E_a g(W)≤1−δ(μ,ν)`，`δ=μ²/ν[1−g(ν/μ)]`。

μ>0、ν≥μ²且0<g(d)<1给0<δ<1。代回四块恒等式与两项AM–GM，因2^-1/2·(2^-1/4)^2=1/2，得到

`M≤[q(1/8)^4/2](1−δ)^2`。

相对R004确切pair上界表达式有严格乘法改善，I因子为(1−δ)^(7/3)。它仍是actual M/I的上预算，不是actual值的区间或下界。

## 4. μ、ν和δ的安全方向

k(d)=(1−g(d))/d。g递减凸且g0=1，所以1−g正、凹、在0为0，k正且非增。δ=μk(ν/μ)，故固定ν时随μ增、固定μ时随ν非增。直接导数同样成立：

`∂μδ=(μ/ν)[2(1−g)+dg']>0`，

`∂νδ=−(μ²/ν²)[(1−g)+dg']≤0`，

因为凸性给1−g≥−dg'。没有把ν上端误当增益上端。

对正q端点A=q1/8、B=q1/16、C=q3/16、D=q1/4，R1下端取B−C−/(A+)²，R2上端取D+/(A−)²。因ν对R1系数−2，ν上端必须用R1下端。作者full-ball实现比这些公式更保守，但方向正确；本审用Fraction直接计算极端有理式，逐项证其R1/R2/μ/ν均被交付区间包含。准确验 `μ_L>0`、`ν_U≥μ_L²`，还验更强ν_L>μ_U²。

δ使用交付μ.lower与ν.upper，是真实μ的下界、真实ν的上界。本审从准确dyadic输入，使用正递推T2−1、T4−1、T8−1与R=T8^(1/4)，按

`1−1/R=(T8−1)/[R(R³+R²+R+1)]`

稳定重算δ，并证明交付dyadic δ_lower严格低于本审768-bit表达式下端。这直接验证lower-bound关系，不能把“下端与真表达式必须overlap”当验收条件。

## 5. R004精确复用和三个新Mellin参数

R004证书原始字节与R005复制件完全相同，哈希与复制manifest及原manifest的certificate项一致；reused_q_1_8字段与原q对象逐项相同。本审没有重算其512面板。R004的q已包含Taylor误差和无限尾；q0=1不需数值积分。最终M使用R004 M_pair.upper，而非重新执行旧q积分；该dyadic upper≥(q_upper)^4/2，本审准确核验，使用它只更保守。

三个新s的参数分别为：

|s|a=s/2|h|ρ|面板数|a−ρ|a+ρ|
|---|---|---|---|---:|---|---|
|1/16|1/32|1/64|3/128|1024|1/128|7/128|
|3/16|3/32|1/32|3/64|512|3/64|9/64|
|1/4|1/8|1/16|3/32|256|1/32|7/32|

全部2hn=32、h/ρ=2/3、0<h<ρ<a。本审源码逐项核对，ledger有1792个唯一(group,index)键且全部completed。各s使用本身的πΓ(s)、M*=Γ(a−ρ)Γ(a+ρ)m(a−ρ)m(a+ρ)、Γ(a+16)和m(a)，不能沿用1/8那组常数。

继承R004的实正moment/Fourier证明，每个r>0半线系数1/(πΓr)正确；复strip |Imt|<a、h(v)=Γ(v)m(v)的log-convex对称端点乘积给全局Cauchy包络。degree128第一遗漏129，raw总误差始终 `32M*(2/3)^129/(1−2/3)`，因总长度32而非面板数n。N16 recurrence给normalized tail `m(a)^2Γ(a+16)^2/[πΓ(s)31·32^31]`，支付全部t≥32。

run.py准确中心(2j+1)h，ctx.prec512/cap129，F的conjugate coefficients相乘是正确holomorphic paired-series，偶项积分使用2h^(2k+1)/(2k+1)，负项/interval radii未被丢弃。Gamma-series库和Terra生成器承担面板数值生产；本审读源码及规格，不再次运行mesh。

## 6. 独立本机检查与真实失败记录

登记先于数学执行。minchecker是CLI可携带、stdout-only入口，参数--numeric-dir/--theory-dir/--r004-dir，无保存结果或修改源操作。默认路径适配合并后项目结构；本审运行时显式传作者独立目录。shell仅将stdout保存为本审results.json。没有调用作者checker，也没有启动新策略。

独立准确Fraction汇总三个ledger polynomial lower/upper，再以768-bit log-Gamma路径重评三组M*/Taylor/tail/q；这些区间均被交付dyadic区间包含。随后准确有理端点方向检查、稳定T8正公式δ复核、原始whole-boundary C^(7/6)H/896系数复核。完整结论是

`I768upper < Terra dyadic_upper < 7730422817512/10^18 < R004_budget_lower`。

M、I的相应lower endpoints只属于确定上预算表达式。小数μ≈.0017452683、ν≈1.032e−5、δ≈.02393288不作任何认证依据；results.json保存全部exact rational接口与源哈希。

本审首次执行遇到python-flint fmpq对象不能直接传Fraction的TypeError，改为其准确rational字符串转换后第二次exit0，无数学断言失败、无参数/计算范围改变。成功脚本0.07105037500150502秒，Python3.15.0、SymPy1.14.0、python-flint0.9.0。共享Arb基础库，不能称不同算术内核。

## 7. Terra checker首失败与修复范围

Terra receipt明确首次checker错误是要求保守dyadic lower endpoint与其所界定的表达式区间overlap。当前check.py改为 `delta.lower()>=dl`，并单独检查0<dl<1；这才是正确的lower-bound关系。本审以不同的稳定δ公式直接证明dl<独立表达式下端，确认修复后的数学验收正确。

prerun-sha256中的run.py、assemble.py、R004 certificate/manifest、Sol proof/registration哈希都与当前字节吻合；它们在checker修复前已锁定。当前certificate与assemble.stdout字节一致，面板账本与finalmanifest一致。没有旧失败checker字节可供本审直接diff，因此不冒称独立逐字证明历史checker只改一行；历史“无面板重跑、无数值更改”还依赖Terra的操作事实说明，而本审已独立核定final ledger/scalar证书有效。此历史边界不削弱final数学关系。

## 8. 完整I接口、未解决项与停止

原Gaussian表示仍同一个共同时间、64未单位化body方向过程。先全boundary shape积分与Minkowski，再用原law的14face×64orthant=896对称，得 `I≤M^(7/6)/(105·2^(25/6)π^(13/12))`。本审系数用原C,H重新组装，没有对单shape face强加对称，也不遗漏short time或singular neighbourhood；无需另加短时预算。

若今后有全域essential g0，下游已审全空间正部source、同一D对全部共同r>0、一般L1 monotone passage给 `C7≤7^(7/6)I0/g0^(7/6)`。R005只提供I0，不能填补g义务，不能声称新C7数值上界或排除原模型。

本批通过后停止，不开启任何新策略、更多panel、额外指数或样本。旧封存文件及其他工作者文件未改；无MATLAB、commit或push。
