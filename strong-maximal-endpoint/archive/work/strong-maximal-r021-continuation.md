# 第21轮接续：一般上侧从log²n到log n？
2026-11-06。用户未叫停；goal active。前轮为实质progress，r019/r020均已完成推导、实际软件核验、独立审查及推送。

## 权威状态
专用clone /Users/zhengzhihao/Desktop/T/Math-Research-strong-maximal，main，当前4267233bfec8e723279a2e779c711678d6f9ab7f，远端包含已确认，工作区干净。
r019 3f13bb1：窄平台w=.5logn+loglogn，J=e^-q[1−e^-w+e^-w I_n]，K→1；qpl=qsoft−.5+o1，固定mu平台score~exp[nmu−dn+a z_n²+b0−.5]。
a=(3−mu)/6，b0=(5mu−6)/6；cshift=mu−sigma/2 sqrt(L/n)额外gL,g=(mu−1)/8。
因此一般下界logC≥nmu−dn+(9−mu)L/24+oL，系数正。
r020 4267233：严格一般upper logC≤log calR≤nmu−dn+O(L²)，完整证明在rounds/r020/report.md。
总结 logC=nmu−dn+O(L²)，一般logn系数和乘法渐近等价仍未解决。切勿mark goal complete。

## r020已证明的关键新工具（无需再猜测）
给任意固定大N，L=logN,t0=C0L，H0=4sigma√(t0L)，h_k=K L/(k+t0)，H累加。
W_k(q)=NormalCDF((q−mu k−H_k)/(sigma√(k+t0)))，V_k=W_k−N^-4−kN^-5。
统一一步Gaussian moment comparison靠D=mu−m≥0, |v−mu|≤Cm√D完成平方；Taylor余项相对density用E|Z−mu|^5 e^{eta|Z−mu|} uniform控制。
选择次序 Cstar→K→C0（继续增大C0不改变Cstar）。explicitGaussian density ratio给右移≥K/4 L/τ f，吸收误差；中央A带映到2A带。
外带tail+base闭合，任意predictable策略ALLq有F_pi(SN≤q)≥NormalCDF((q−muN−HN)/(sigma√(N+t0)))−2N^-4，HN=O L²。
积分只到T=(mu+1)N再uniformpositiveMGF处理远尾，得E(S−s)+≤d+O L²，避免旧Vmax负excess代理精度不够问题。
profile中段b=√N L²，far tail exp−cL4，给可行lowfailure≤d/N+O(L^2.5/N)。
测试s=Nmu−dn+K0L²，Gaussianbarrierfailure≥β+cL³/N，β=dn/N，矛盾得一般upper。
两独立审查upper_audit/global_upper_adversary均ACCEPT全文。所有软件真实执行；无live shell/MATLAB作业。

## 当前代理
/root/global_upper_adversary刚收到r021候选任务，可能仍running；请先收其消息再重复研究。只外部work文件，不改repo/git。
upper_audit已完成r020，可在有正式新证明后再唤醒审查。multilayer_theory已完成r019。

## r021候选：状态相关曲率补偿
r020把z²≤A²L统一上界造成ΣL/τ=O L²。
精确热算子identity（已本机Wolfram15实际运行/work/state_barrier_heat_identity.wl，符号残差0）：
f(t,x)=Normal density variance mu t，H(t,x)=K[x²/(mu t)+2log(t/t0)]，
 (∂t−mu/2∂xx)(H f)=K[1+x²/(mu t)]/t * f。
这正匹配momentcomparison一步损失C(1+z²)/t f。
候选barrier：x=q−mu k−B0L，t=k+t0，W≈Phi(x/(sigma√t))−H(t,x)f(t,x)，或更安全shift Phi((x−H(t,x))/(sigma√t))。
在|z|≤A√L中H=K[z²+2log(t/t0)]=O L，H+B0L总位移O L；若能离散化且余项被正heat residual吸收，可得一般upper O(logN)。
选C0≫K²，H/√t相对z小、density比常数；t≥C0L使时间/空间高阶误差有small参数。对非线性shift的heat残差需完整算，不能只用形式PDE。
外区把z²光滑截断在O(A²L)，在更大strip接上constantshift；所有patch点需要残差tail≤delta，保持base q<0 W0≤epsilon。
直接线性W=Phi−Hf可能右带为负或不保CDF单调，但subsolution本身可以signed，只需V0≤indicator和bounded；是否超出1/negative会影响farright和stochasticdom积分，需要仔细。
当前仅候选，尚无完整离散证明。代理正独立检查。

## 软件与归档
r020 code sm_gaussian_barrier_checks.m解析Exp卷积erfcx稳定计算63状态，每点优化c及端点，JSON只存每horizon最小摘要（3行）；MatlabAnalyzer clean。
Wolfram gaussian_barrier_checks.wl三符号残差0、三80位quad吻合45位；参数C0=1e4,K100只是数值诊断，非认证证明常数。
每轮只最终报告/必要代码/最小摘要，单独commit，fetch/rebase/push后确认远端包含。

## 下界下一线（未立新定理）
完整KKT profile比固定斜率平台更灵活：a(t)=e^t fS(t) monotone，optimal h=0 ifa<ν；h=1 ifν<a<nν；soft满足Phi'(h)=a/ν。
ν需使Eh=1且growthrate<1保证Jfinite。可直接Fenchel支持线证明最优，不必泛引无限维强对偶。
当前窄平台只是一类可行输入，不宣称Q全局最优。后续可本机计算KKT或平台+斜率α/n的参数搜索，但应同时保持真正generalupper进展。

## 最新审计修正（重要，优先于上方候选）
global_upper_adversary找到linear W=Phi−Hf实质blocker：t0=C0L时左中尾z~t^(1/6)已W_x<0，tail exp−cL^(1/3)远大于delta，不能patch。不要再使用linear近似作为屏障。
但nonlinear shifted W=F(t,y), y=x−H(t,x)并未被此否定。主代理随后发送exact公式，代理已重新激活验证：
L W=(∂t−mu/2∂xx)W
=[−(H_x−H_x²/2)y/t +mu H_xx/2−H_t] F_y(t,y)
=[−K/t−Kx²/(mu t²)+2KxH/(mu t²)+2K²x²(x−H)/(mu²t³)] F_y。
若按new density F_y衡量，错误项相对K(1+z²)/t由K|z|/√t≤KA/√C0，及K log(t/t0)/√t≤2K/(e√t0)控制，固定大C0可能足够；不需要linear化所要求的H|x|/t小。
W_x=(1−H_x)F_y，strip内H_x小所以positive。需验证exact公式及所有离散Taylor时间误差/patch，尚未定理。代理正在此具体任务上研究，先收其结果。
