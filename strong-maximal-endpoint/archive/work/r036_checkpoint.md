# 当前恢复入口：r036已完成并推送

2026-11-06；持续goal active，用户未叫停。上轮为实质progress：自然管内nonlinear MAX全动作引理已闭合，非仅条件接口；globalouter仍open。勿标goal complete。
专用clone /Users/zhengzhihao/Desktop/T/Math-Research-strong-maximal，main干净。r036=48f3b6559d112cc83f63dd5f9ae903e4ec0c904b，commit/fetch/rebase/push及最终fetch/ancestor核验全部完成；远端origin/main包含本轮提交。12必要files，项目202files~0.93MB,max18.5KB。Formal authoritative rounds/r036/report.md，3WL+3JSON+README，r034维数勘误3files、r035libraryflag2行。无其他项目改动，无Matlab启动，所有本轮Wolfram已退出。

**维数勘误：receiver_seed实际n64,b32，不是n256/b128！256是外层quadrature节点数。r034report错误句及JSON/code已本轮更正；原计算始终读取n64，数值无改。upper horizon N256是真的256，与lower n64区别。历史checkpoint旧误述以此更正为准。**

r036 strict progress（详细proof见正式报告）：

1. φ_N正logconcave，左非L1用e^-eps|r|截断后极限，common-a f_t logconcave t≥4 ⇒v=φ*f正、递减、logconcave。δ=-v'/v nondecreasing，v(r-z)/v(r)≤e^(δz)。自然管δ≤C sqrt(logt/t)，Mills统一a,N。Doobdensity q_c v(r-z)/(Kcv)≤2e^{-(1-δ)z}≤2e^-z/2，uniformallc（不是FOSD）。
2. exact auxiliary identity f_t C_t=t g_(Gamma(t+2)+aBin(t))(r)+a²t p(1-p) g_(Gamma(t)+aBin(t-2))(r-a)。由当前tilt全局derivativebounds，U=αJ+s_t f (0≤s≤1)满足可达左移weightedspatialderivative界；w=U/v有sup O(log^1.5t/sqrt t)，weightedLipschitz L=C(f/v)log^1.5/sqrt t，κ<1/4。左尾用固定小negative offcenter tilt（非精确saddle）：sup tiltedGamma mixture≤C/sqrt t，得所有r'≤tm的density Gaussianupper。Jensen实际payoffdefect≤Cf_t log4t/t^1.5 uniformlyallc，reference-semigroup可求和O f_N。
3. **LOCAL MAX closed**：α>1任意固定；存在B,Λ,t0，s_t=BΣ_{j>t0}log6j/j^1.5，U_t=αJ_t+s_t f_t，V_t=v_t exp(U_t/v_t)，W_t=V_t+Λt b²/N（a matched2a-m=q/N,b=m-a,ell=2a-m）。自然管|r-tm|≤A sqrt(tlogt),t0<t≤N，sup_cKcWprev≤Wt。关键F=ve^w-v的全左移derivatives：F'≤Cf(log1.5/sqrt t+log3/t), F''≤Cf log3/t, F3≤Cf log4/t1.5。精确矩消去 A_F=p[PF-F(r-a)]=p[-bF'+.5E(Z-a)^2 F''+R3]=O(f log3/t)，不是粗U'。strict meanmax给exterior动作gap−κd，Ftransport=o f；compact强凹定位perturbedroot，Youngerrorf log6/t²。Referencegap Vt-v_t exp(PU/v_t)≥γ_t exact；B先选大，再t0大使BΣtail≤1，避免常数循环；Λ floor抵消Cb²/N。criticalW_N(q)=v_N(q)+αg f_N logN+o(f_N logN)+O f_N；exphigherpowersL4/N²=o f_N。**仅局部，不能跨exit迭代冒称globalupper。**
4. fixed-positiveη nonlinear necessary HJ test：Jrate=I_mu-αD，D=p/8(y-mu)^2η(1-η)，H=J-yJ'+Λstar(J')。A3=mu(mu-1)/8，Bλ=mu(mu-1)(9-2mu)/48，Bgap=mu(mu-1)(99-43mu)/384，BD=A3(2-mu)，K4=mu(mu-1)(3+5mu)/384>0。α1 H=+K4η4；α>1 H=-(α-1)A3η3+...只通过smallfixedeta必要检验，不是fullBellman。α1.01在η.4仍positive，α1.05 samplednegative（非interval）。η→1若α<8/p则H~log(1-η)→-∞。
5. numericalfast36稳定logsumGamma mixtures，6spots vs50digits maxabs1.1e-14，densityrel<1e-12。clippedOleft nonlinearα1.1: βH2f仍positive；探索log4 summable-source β.01/.05,t16/64/128,h.025 positive≤2.3e-16；refine128 h.0125 coremax−3.8156e-8，terminalN=t256 h.025 globalulp2.22e-15 farleft/coremax−6.0598e-9。gridmono。仅PWL machine诊断，tinyprobunderflow0，不是cert。数值log4source不等theorylog6+t0构造。
6. lower n64 jointdeformedcurve+source实际体积已结算。finitepositiveconeunion acceptedavg>1+1e-5节点和mirror，V=Σ(e^-xi-e^-xnext)e^-prefixminy；r026 K_L无c因子！c在核归一化已含。Base V16.684838334873092,score71.27164607990836。outer最高minavgα.015 V15.63217,Jscaled.06223805807938，score71.2090468732；真正最高score α0=71.20932898008，仍低base .06231709983。不要按maxminavg选scorewinner。Middleα-.015,0,.015只12/11/11groups可接收，scores70.45271/70.41438/70.492665。c0 scalar numerical reference71.71499208794934（非fullQopt，非recC same-c）。新middleα.015 source上3true-dependent directions（c1/c2=recC±.1v/16镜像交换，common c(v)=recC+.1(v/16)^2）均score70.3760618<uniform70.492665，samecost.06275686942107。无winningcandidate，停止此familyfinefit。源weights exp reweight仅first-order massbalance，非exacttotal。

全部终审通过：baselineDoob/Jensen，upperlocalclosure，globaladversaryrate，directedgeometry。work/r036_local_closure.md开头旧notclosed已被后面acceptedlemma取代，不要再用旧头当block。Formalreport包含完整proof及范围。

## r037两项任务（未证明的想法，勿当结论）

UPPER：globalextension是现在唯一upper主缺口，不再重做local。正moderateη∈[A sqrt(logt/t),δ(α-1)]：H margin−εη³，有限时η/t误差比η³为1/(tη²)≤1/logt可能可吸收。w=tD可能大，Doob Jensen需用单侧oscillation w(r-z)≤w(r)+Oη² 而非wholewsmall；大的η外repair用精确Chernoffscreen，weighted-L1 outer source的exp-t积分可能O f_N。不能只min一个underpredicting referencepoly。
LEFT coarseclip：a=L-v≈F，qratioU/a≈ακ tλ³。coarse cuts q≥1；activeq≤1⇒λ≤C t^-1/3。失败costC=F-U≈F(1-q)，C'≈f[(1-q)+cλ]，C''≈λf[(1-q)+c'λ]，故gain≈fλ²[(1-q)+Oλ]≤(1+ε)reference source。近clip leakage -Kc[(I-O)+]有利，需严格Mills+离散边界层，而非假设全局I单调。
备选two-sided logistic uplift：capL=Oleft≥v,a=L-v；W=L*logistic(log(v/a)+U L/(v a))。一阶v+U，右v≪a变v exp(U/v)，左a≪v变L-a exp(-U/a)，自动成功侧饱和。criticalU/a~L^1.5/sqrtN，extraU²/a=o f_N。a=0需极限定义，capdrift与globalBellman均未证。可用positivecostreferenceP ψ_+（ψ在[0,rho]正凹⇒logconcave）给第二Doob结构；v+Pψ_+=L^a+P b，另加drift/Γrepair。这是candidateonly。
LOWER：当前jointcurve6cases+3trueconditionaldirections全部负验收，别再finefit。需新receiver集合/比例耦合或非锚定crossrect，所有试验同时算实际Phi_n预算、acceptedunionV、同维scalarbaseline。Fixedglobalc路线r028已排除为突破Q主线。

一般C上下仍差OlogN。即使globalW闭合可匹配R的glogN模型系数，真实C下界仍缺机制，goal不得complete。用户未叫停，下一轮继续。

以下旧轮历史仅供引用：
