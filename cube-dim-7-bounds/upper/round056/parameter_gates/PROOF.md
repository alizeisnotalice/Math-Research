# R056 · 新proposed tuple的准确常数与L≥2^20有限gates

## 状态与结论

本报告独立核验一个新的proposed参数包；**不声称新的C_n上界**。准确预算、coordinate moments和finite gates通过之后，仍须作者的新family完整operator证明及独审，才能把本包激活。旧R051的小d参数域不能自动套用到d=3/10。

准确Taylor lower给c_raw lower>1/40，因此安全展示系数可取1/40，无需向下改。完整scalar根式也已计入，见第8节。没有参数搜索、生产批次、MATLAB、旧grid或Git操作。

## 1. 参数与完整预算

L=logn≥L0=2^20，x=(93/100)L，t=n/x，
d=3/10，σ=7/100，χ=7/200，h=2^(−24)/(Ld)，
ζ=n^(−3)，κ=1+1/t−ζ，θ=exp(−χx)，k=nhθ，
qmin=31/40，qmax=37/40，
rlo=43/40+2^(−20)，rhi=49/40−2^(−20)。

εJ=εN=εRet=εS=εT=εD=1/1024，εA=1/64，P0=65275/65536，M=21/8。

定义

Ulo=P0(1−εN)/(1+εJ)，Uhi=(1+εN)/[P0(1−εJ)]，
b=(1−εRet)2Ulo/qmax−εS−εT。

准确有2Uhi/qmin<M、b>0。完整target mass为m_true=(1−εA)log(rhi/rlo)。令u=(rhi−rlo)/rlo∈(0,1)，用

log(1+u)>u−u²/2+u³/3−u^4/4

得到严格有理lower

mL=259465051078193863093554363/2018111068113594386771608802。

取ΓL=mL b/M、γL=(1−εD)ΓL。它们是可安全使用的预算lower，不能误称含log的m_true本身是有理数。准确有

γL=172531397335235056462877787662549380239/1643616601687301779182642153860479582208，
cL=γL−σ/(1−σ)=γL−7/93
 =4540103740365747796769139175593735286771/152856343956919065463985720309024601145344
 >1/40。

其准确正余量为

cL−1/40=3593475707213855800847480839340601290687/764281719784595327319928601545123005726720。

本轮不通过显示小数重建任何参数或结论。

## 2. 全q geometry、g与pin log

写reference R(q)=q+d−h、a0(q)=q(1−ζt)、g=R−a0=d−h+ζtq。Taylor给n≥L²/2、n≥L^8/8!。由x≥L/2，ζt/h≤2d/(2^(−24)n²)，基点严格≤1/2，随后递减。

L0处h<min(2^(−20),1/10000,d/8)，故所有L≥L0均有：R≥43/40−h≥10749/10000>1，0<g≤d<R/2，κ>1，t≥L/2，ρ=tq>4rhi。target inset使每个r∈[rlo,rhi]的完整matched interval

q∈[r−d−h,r−d+h]

包含于[qmin,qmax]，因此没有端点长度损失；其source log宽在[2h/qmax,2h/qmin]。

source pinned factor ell_*=1−g/(2R)，所以

−logell_*≤log(3583/3083)。

这里3583/3083正是R=10749/10000、g=3/10时2R/(2R−g)，不是小数拟合。

## 3. Baseline log与两count的四事件概率

baseline q_b=r−d+h，使R(q_b)=r，free-coordinate EJ=1恢复准确。baseline log目标为|ΣlogJ|≤1/256。

相同的非渐近coordinate moment界给

|EΣlogJ|≤[6d³L²/qmin²+dL/qmin]/n，
Var(ΣlogJ)≤12d³L²/(qmin² n)。

用n≥L^8/8!上包后，两者随L递减。基点分别≤1/512、2^(−26)，故Chebyshev坏率≤2^(−26)/(1/512)²=1/256=εA/4。

设q_±=r−d±h。每个source面数和Nsum(q)=Σ_i[1+1_{|w_i|<a0(q)}]，固定target面的一坐标必为1，其他坐标给

E_±=n+(n−1)a0(q_±)/r，Var Nsum(q_±)≤n/4。

两事件均取|Nsum(q_±)−E_±|≤εN n/4，各自坏率≤4/(εN² n)。n≥2^30时均≤1/256；L0处n≥L0²/2=2^39已保证。

第四事件是原target corner cap。rhi<5/4、t≥L/2保证ell_t(r)<3/4；n≥33给(3/4)^(n−1)≤(3/4)^32<1/256。四个坏事件并集≤εA，所以λ(A)≥m_true n≥mL n。

这些是完整uniform壳概率，不假设实际increment measure服从随机分布。

## 4. 两count校准与matching product的准确下包

Nsum与ell_*都随q递增。R_-=r−2h、R_+=r，a_±=(1−ζt)q_±。端点两事件给任意matched q

Nsum(q)/(2n ell_*(q))≥1−[2h+1/n+εN/4]，
Nsum(q)/(2n ell_*(q))≤1+[2hd+εN/4]。

第二式用

a_+/r−a_-/(r−2h)=2h(1−ζt)(d−h)/[r(r−2h)]≤2hd。

基点2h+1/n≤3εN/4，此后保持，故calibrated ratio在[1−εN,1+εN]。这里不把d=.3时真实lost faces错误地当成o(n)。

baseline到matching的q距离≤2h。coordinate导数界在g≤R/2、R,q≥1/2、ρ≥4R、T≤2d的正域使用。因为d=.3有K_w≤n≤4nd，原导数求和仍给log bridge≤1280Ldh=5/65536。于是

|ΣlogJ|≤1/256+5/65536=261/65536，
P0≤∏J≤1/P0。

此外pointwise |J−1|≤2d/ρ≤2dL/(qmin n)，基点≤εJ，尾部递减。与calibrated ratio合并，完整inverse row系数严格在Ulo n与Uhi n之间。

若新family采用完整target面耦合，L_i<1时恰丢一个面且L_i>1/2，故完整target shell系数≤n∏L_i；matched列≤2nhG/R。因此加soft weight后column≤2k/R≤M k，准确gate中也核2/(10749/10000)<M。该改善列公式应由新family原实例证明保留，不能只引用旧4k预算。

## 5. 精确inverse-coordinate一二阶moment

令W=−logell_R(Z)，Q_act=−logell_*。在一个captured source face上，free坐标的inverse-capture正规化因子J≥1−εJ。丢弃capture指标作upper，并换元w=−logell，有准确积分

E_full[W/ell]=R Q_act²/ρ，
E_full[W²/ell]=2R Q_act³/(3ρ)。

因为nR/ρ=x(R/q)≤(43/31)x，free-coordinate和的mean与variance分别至多

meanrate·x，varrate·x，
meanrate=[43/31]/(1−εJ)·Q²，
varrate=(2/3)[43/31]/(1−εJ)·Q³。

其中Q取pin log的准确Taylor3 upper：z=500/3083，

Q=z−z²/2+z³/3=13226208500/87910718361。

Fraction直接核Q<1/6、meanrate<4/125、varrate<1/300。pinned source坐标只加mean≤Q而不加variance。

重要scope：该variance是每个captured source face条件律下的界。不能未经证明称整个face mixture的variance也≤varrate·x；下节先在各面应用同一Chebyshev界，再混合即可。

## 6. 修复pin预算后的retention gate

正确分配是pin/x≤1/2000，而非1/400：Q<1/6、x≥1000已给Q/x<1/6000<1/2000。于是

χ−4/125−1/2000=1/400。

每个captured source face都有threshold χx与mean至少相差x/400，故

P_inv{G<θ}≤(x/300)/(x/400)²=1600/(3x)。

L0处x=(93/100)L0≥1000，且准确1600/(3x)<εRet；此后按L^(−1)递减。对各面混合仍保留同一失败率上包。

若错误地只取gap=1/2000，则概率会成为40000/(3x)，不能沿用1600/(3x)。等式χ−4/125−1/400=1/2000本身正确，但不能配成所需概率常数。本报告不使用那个错误搭配。

## 7. Source cap与τ的全尾

source q范围给min v/(2q)=3/74。n≥75时(1−1/n)3/74≥1/25，所以原source cap质量≤nd exp(−x/25)。由于θ=exp(−χx)、1/25−χ=1/200，

Cap/k≤(d/h)exp(−x/200)。

L0处d/h=2^24L0d²<2^44且x/200>2^12。Taylor8给基点upper 2^44·8!/2^96<εS，全尾按L^(−7)递减。

完整scalar所需τ≤n√L exp(−σL/2)给

τ/k≤(√L/h)exp[−(σ/2−χ(1−σ))L]。

准确σ/2−χ(1−σ)=49/20000>0。L0处prefactor<2^54、指数argument>2^11。Taylor16给upper 2^54·16!/2^176<εT，全尾按L^(−29/2)递减。

这些是proposed family需要的finite gates；它们不单独认证weighted能量引理、same-history或general-L1 passage。

## 8. 完整scalar余项与安全c

若new family确实供给预算γL n，则原cap scalar estimates在x=(93/100)L给

F(A0,B0−γL n)≤T−cL n+nR(L)，
R(L)=1/n+e^(−σL)+2√L e^(−σL/2)。

γL<1且B0>8n使根号预算正。L0处σL0>2^16，Taylor8及n≥L²/2给

R(L0)≤2^(−39)+8!/2^128+2^11·8!/2^120<cL−1/40。

全尾按L^(−2)、L^(−8)、L^(−15/2)递减，故完整scalar确实允许1/40，而不是只看负缺额γL。

## 9. 证据与停止身份

`verify.py`使用Fraction及最小SymPy积分检查，实际核全部注册gate与exact log上下包。receipt保留输入契约、来源及hash。作者新full family尚未接受时，本包始终是proposed-gates结果，不发布新Cn。

下一步仍须新family完整确认校准Ni、精确column耦合、inverse row probability及原finite-a/generalL1桥接。任何一项失败，均不能由本包常数通过覆盖。旧R045/R051 operator定理的小d域未被本报告默默扩展。
