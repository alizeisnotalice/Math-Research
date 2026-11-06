# 第23轮接续：独立机制已闭合后的相关输入搜索
2026-11-06。目标active，用户未叫停。前一goal turn为progress：r022严格定理+独立审查+本机验证+推送。

## 权威工作树与结果
专用clone /Users/zhengzhihao/Desktop/T/Math-Research-strong-maximal，main=1601f75108c7cf525e47d1a3eb8d00d0b3d02eab，远端包含已确认，工作区干净。
项目strong-maximal-endpoint/；其他/Desktop/T/Math-Research脏工作区勿动。
r021一般upper已证明logC≤nmu−dn+O(logn)。
r022已严格证明（双审ACCEPT）
 logQ=nmu−dn+(9−mu)logn/24+o(logn)
适用于任意确定性独立比例vector+任意非减标量profile，非仅等ratio或有限阶梯。
一般C只有相同lower与O(logn) upper；最优logn系数和乘法渐近等价仍未知，绝不mark goal complete。
mu=logB=1+1/B,sigma=sqrtmu,barPhi_G(z_n)=sigma z_n/sqrt n,dn=sigma√n z_n。

## r022关键证明（无需再列未闭合bootstrap）
正式在rounds/r022/report.md。粗预算+低band s−2L先排D≥2d，再BE得D=o√nL,y=M−s~sigma√nL。
负tilt λmean=s，uniform triangularLLT g0~1/√n,supg~1/√n,supg'~1/n，给f/F~λ与densityratiof(s−u)≤C f(s)e^-λu。
精确Gamma分解S=Gamma(n,1)+B≥0：R=E e^{S−s}1S<s=e^-s E[e^B(s−B)_+^n]/Γn+1 ≤(s/n)f(s)。
因此F≤Cρ，再negativeexcess≤Cf/λ²=O1，保留D的逆模数二阶预算F(s)≤(d−D)/n+O L/n。
统一独立三角阵Fourier Edgeworth absO1/n通过Σlogchar+每步Exp damping证明，无需同分布。
先implicitnormalroot得到s≤M−d(V)+aL+oL，再variance stabilityforcesD=O L。
d(V)=dn+λ_nδV/2+oL；global minHλ(ri)=μ−m+λ(v−μ)/2 =−μ(μ−1)λ²/8+Oλ³，求和给gL补偿。a+g=(9−μ)/24。
上界不能转给反馈或一般C，独立char product是实质scope。
MATLAB真实运行n2…256 two-groupexactGamma-binomial CDF与Edgeworth及lowR bound，Wolfram符号残差0、costcoeff与80位GammaGaussian反例核验，Analyzer clean。

## 当前代理及尚未完成的新probe
/root/directed_certificate 正在新相关网格probe任务，已确认agent running。它只写work外部，不改repo/git。先收其消息，不要重复启动同一MATLAB作业。
主代理没有live工具句柄；代理可能正在执行求解，需以其实际句柄/结果确认状态，不能仅凭agent running断言MATLAB运行。
其ownership：/Users/zhengzhihao/Desktop/T/work/sm_correlated_grid_probe.m、同前缀JSON/notes。
任务：本机MATLAB intlinprog，总solver预算30–60秒，5–10sec/iteration最多6次Dinkelbach。
输入4×4单位cells[0,4]^2，每cell选rational幅度[0,1,3/2,2,4,8,16]。
receiving非均匀axisbreaks[-16,-4,-1,0,1,2,3,4,5,8,20]，halo cell weights真实面积。
rectangles所有grid边界，prune最高幅度质量仍不能使avg≥1.001者。
binaryonehot amplitudex、certrectanglez、halow；constraints实际mass≥(1+1/1000)*area*z (无需BigM),w_i≤Σcover z,至少1halo。
Dinkelbach maximize haloArea−λ actualOrliczCost(n=2,amp(1+log^+amp))。
必须记录solverstatus/gap/tol/time，roundinteger后独立重算所有strictrectmargin/union面积/实际分母。
同一grid枚举singlefilledaxisblock×constantamp基准，判断相关input是否有finitefamily增益。此不是全局C_n改进，也不是超出asymptoticQ的证据。
特别n2已知smallthresholdcontrolC2≥4；候选若<4只能作为有限网格/形状probe，不能登记为新bestC2 lower。
若没有严格新机制/反例，不机械写“轮次完成”报告；保持probe继续或转更有区分力路线。

## 后续证书方法
rationalamp与rectarea使strictavg与coveredvolume可用Fraction精确核验。
Orliczdenom upper可用logA的atanh有理级数，u=(A−1)/(A+1)，logA=2Σ_{k=0}^m u^(2k+1)/(2k+1)+R，0≤R≤2u^(2m+3)/((2m+3)(1−u²))。
大A先按2幂range reduction可加速（r004已有工具）。certified lower=haloArea/denomUpper。solverincumbent不是global finitegrid最优证书；浮点gap也不能当精确对偶。

## 一般upper初始化机制（r022 §8 已登记）
常数GaussianbaseW0=Phi((q−b)/sigma√t0)，terminalε=N^-p要求b≥sigma√t0 Phi^-1(1−ε)~sigma√(2p t0 L)。
在r021全带weightedmoment条件t0≳L下，总effectiveinitshift必ΩL，非锐C0污染sharpcoef；不能只锐化localK。
Gamma2CDF在q≥0为1−e^-q(1+q)，q<0=0，严格≤one-stepinfCDF tanh(q/2)。
但direct Gaussian热拼接在q<0失败：actual E G2(q−Z)=0,Gaussian(mu,mu)卷积>0。
正确q−1值0.001895398038066558…，旧agent曾把变量均值写反所得.159659918已修正且未入repo结果。
正支撑Gamma reference / exactshorttimeBellman splice仍未证明。不要把scopeobstacle推成所有方法不可能。

## 技能与工作方式
本turn已应用并读取math-l04-mathematical-programming-workflow SKILL及method/provenance，已告知用户。它仅约束模型翻译、求解状态与独立验证；其既有账本仅2变量LP证据，不提供本次MIP全局保证或本机API保证。
位置 /Users/zhengzhihao/.codex/skills/math-l04-mathematical-programming-workflow/SKILL.md。
用户偏好额度高效、大量本机数学软件；不要为简单重复验证滥增代理。
每轮最终报告/必要code/最小结果+source白名单，onecommit，fetch/rebase/push，确认remote包含；未完成试验留work外部。
