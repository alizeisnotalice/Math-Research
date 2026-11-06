# 第22轮接续：对数级系数与独立比例机制
2026-11-06。用户未叫停，goal active。前一goal turn是progress：r021一般上界严格改善并已推送。

## 权威工作树
专用clone /Users/zhengzhihao/Desktop/T/Math-Research-strong-maximal，main=2218fdfe23d4d8a01904743f93b03cc2fb6465e3。
远端origin/main包含已确认，工作区干净。项目strong-maximal-endpoint/。
原来的/Users/zhengzhihao/Desktop/T/Math-Research脏工作区勿动。

## 最新严格成果
r019: 平台width .5logn+loglogn，使真实模数factor→1，qpl=qsoft−.5+o1，固定mu平台score~exp[nmu−dn+a z_n²+b0−.5]。
a=(3−mu)/6,b0=(5mu−6)/6；固定ratio小偏移cn=mu−sigma/2 sqrt(L/n)使softroot额外gL，g=(mu−1)/8。
实际下界logC≥logQ≥nmu−dn+(9−mu)L/24+oL。
r020: 一般upper nmu−dn+O(L²)（已推送4267233）。
r021: 一般upper已严格到nmu−dn+O(L)，已推送2218fdf，两独立agent全稿ACCEPT。
因此一般logC−nmu+dn夹在两个正常数倍logn之间，最优系数及一般乘法等价仍未知，绝不mark goal complete。
定义mu=logB=1+1/B,sigma=sqrtmu，barPhi_G(z_n)=sigma z_n/sqrt n，dn=sigma sqrt n z_n。
Q≤C≤calR桥在r011已严格证明；Q是固定比例vector+标量monotoneprofile见证，calR允许predictable feedback，不能混同。

## r021 proof关键点（权威在rounds/r021/report.md）
linear Phi−Hf在左中尾失去单调性，已明确失效；不能再把它作为候选。Exact nonlinear shiftedCDF是成功路线。
H=K[psi_L(z)+2log(t/t0)]，z=x/(sigma√t),t0=C0L,x=q−mu k−B0L。
psi=z² |z|≤2A√L，C3 taper到constant44A²L/7 |z|≥3A√L。
W=Phi_G((x−H)/(sigma√t)), global H≥0≤CKL, |Hx|≤1/4。
Exact LW /F_y=−K/t−Kx²/(mu t²)+2KxH/(mu t²)+2K²x²(x−H)/(mu²t³)≤−cK(1+z²)/t in wideinner，C0≫K²A²。
按NEW shifted densityF_y比值，relativeweightedspaceTaylor仍uniform，D=mu−m≥0与|v−mu|≤Cm√D完成平方给actual−Gaussian≥−C*(1+z²)/t F_y，C*在C0大条件下独立K。
ExactDuhamel P1W(t−1)−W(t)=−∫0¹Ps LW(t−s)ds。
|U|≤1固定正概率子事件取≥cK...正收益；其余仍在wideinner贡献有利，只有escape到外cap区用e^-c tL≤e^-cC0L²损失δ。
base B0L=4sigma√C0 L，outerbands tails，V=W−N^-4−kN^-5，backward arbitrarypredictableπ。
得到全状态CDFdom normal meanNmu+O L,variance mu(N+C0L),error2N^-4。
积分至T=(mu+1)N+uniformmgftail得excess d+O L；profilemiddleb=√N L，budgetfailure≤d/N+O(L^1.5/N)。thresholds=Nmu−dn+K0L，Gaussianfailure≥β+cL²/N，矛盾。一般upperO L闭合。

## 本机软件与状态
r021 MATLAB sm_shift_barrier_checks.m已实跑，exactPDE scans L10…1e6 positive，实际一步GaussLaguerre64/128+stableCDFdiff9states positive；差~4e−13。
Wolfram shift_barrier_checks.wl exactPDE残差0、C3matching6残差0、80位actualkernelquad与MATLAB fixedmu差≤1.67e−13。
CodeAnalyzer clean。无live shell或数学进程。
诊断参数K100,C0=1e10非认证存在性常数；严格性靠报告分析证明。

## 两个正在进行的代理任务（先收消息，勿重复探索）
/root/multilayer_theory刚被followup激活，关闭Q sharpL coefficient proof：期望
 logQ=nmu−dn+(9−mu)L/24+oL。
已有初步笔记/work/strong-maximal-r021-q-sharp-log-upper-audit.md（文件名历史r021，下一正式轮应r022）。
/root/global_upper_adversary已followup激活，审计generalupper sharpcoeff的Gaussianbase机制：t0=C0L和baseB0L本身O L，不能仅取K接近局部最小就宣称sharpcoeff；研究positive-supportGammaCDF/exactshorttime Bellman替代base。
upper_audit已完成r021，等正式新稿再唤醒审查。

## Q sharp-L route：已发送给theory的关键bootstrap
Independent triangularvector，meanM=nmu−D,varV，|V−mu n|≤Cm√nD，uniformmoments≥n variance。
Uniform Edgeworth absO1/n可直接r018同法：Σ log char standardvarV≥n，每步Exp damping，local n*t4/V²=O(t4/n)，outer exp−cn，不需新文献。
PointwisebudgetJ≤e^-s gives h≤Phi^{-1}(e^{t−s})。
低侧别用s−2L（会丢失sharpL）：E(1−h)+≥F(s)−∫e^-u f(s−u)du。
TiltedLLT g≤C/√n,gcenter~c/√n gives f(s−u)≤C f(s)e^-λu，so低积分O f(s)；f(s)~λF(s) atcriticalλ~ρ/sigma。
粗预算Epositive≤d_eff+O√n、independent MD first forcesD=o√nL,s<M,λ~ρ/sigma，并给F(s)≤Cρ。
再由tiltedtail integral E(s−S)+=O(F/λ)=O1，保留meanD。
高侧inversePhi二阶 h−1≤u/n+(1+o)u²/n²，中段+far控制得F(s)≤d_eff/n+O(L/n)，d_eff=M−s。O f(s)=O L/n已吸收。
独立Edgeworth转阈值：s≤M−droot(V)+aL+oL，a=(3−mu)/6，droot(V)满足barPhi(z)=√V z/n。
先用V stability：s≤nmu−dn−D+C√(DL)+aL+small，nearoptimal s forcesD=O L。
随后κ3sum=nκ3mu+O√nD，κ3mu=mu(3−mu)；varroot展开droot(V)=dn+(λ/2)(V−mu n)+oL。
最后用r015 exact global minLaplace取λ~z_n/(sigma√n)：
λD+λ²(V−mu n)/2 ≥−n mu(mu−1)λ³/8 +κ3deviation λ³/6+O(nλ4)。
divideλ且D=O L givesD+λδV/2≥−g z_n²+oL，g=(mu−1)/8。
得s upper nmu−dn+(a+g)L+oL，若所有bootstrap严证则匹配r019 Q lower。
目前theory初始笔记仍标3缺口，最新followup已要求关闭，不可提前当定理。

## 一般C sharpcoeff路线的真实障碍
r021 base Gaussian smoothing+ε=N^-p要求initialshift≥sigma√(2p t0 logN)(1+o)。uniformglobalstriprelativeTaylor需要t0≳L，所以base必付ΩL，不锐常数掩盖真实coeff。
形式Bellman diffusion曾给adaptiveCDF quantile feedbackgain g[z²+logn]+O1（仅heuristic，严禁当已证）；fixedvectorgain仅g z²。
可能Q Lcoef a+g=(9−mu)/24，calR Lcoef a+2g=(3+mu)/12，但实际C位置未知，不能把calR当C。
研究应保持真实strongrectangles与feedbackbridge差别。

## 提交纪律
每轮最终报告+必要code+最小summary+source，独立commit；本地未完成草稿留work外部。明确证明/数值/家族障碍。fetch/rebase/push后确认远端包含，不强推、不动其他研究目录。
