# 第030轮：共同选择历史给出的严格帽区扣减

工作始于2026-10-27，2026-10-28封存。用户已恢复持续上界目标，goal实际active；下界继续暂停。接续f65a5bf9的真实工作区，没有重启旧三四维项目或旧计算句柄。本轮主负责证明整合、实际MATLAB交叉诊断、最终准确比较及统一提交。Sol负责理论/独审，Terra负责本机计算，Luna负责有界来源与整合检查，各自独立clone；不将子任务当完整轮次推送。

## 本轮主结论

固定Lebesgue、中心轴平行立方体、全部坐标共同半径、一般L1输入、所有正尺度。本轮严格得到

C7 < 12277313907857/500000000000 = 24.554627815714。

参数t=3391749/524288沿用旧点，新增sigma=2081927/6783498。与同t旧完整F区间下端严格分离约0.000332710251031；不是仅把显示上界最后几位取紧。新结果来自真实共同选择历史的预算扣减，不是单纯同族调参。完整B、P、D和新旧成本双端见audit/round030/numeric/independent_n7.json；原R028完整Arb下端也已准确比较。n4,5,6,8,16,64各有一个固定旧t的新机制认证，完整参数和严格分数见accepted_bounds.json及compute/n30b_certificates.json。n64改善约2.66e−22，由分数比较得到，不由显示小数推断。各点均有a=2^512的有限参数实现；新模型参数全域最优未证明。

## 完整机制：两区域共同历史cut

保留原有限算法：非增半径、固定非负权重wj=|Ij|/(2arj)^n，真实不交增量Ij=Qj减先前并集，资格sum前序wi 1_bQi(yj)≤tau，b=a(1+1/t)。被拒覆盖的首失资格前缀给

|并输入Q|≤[1+(1+1/t+1/a)^n/tau] sum|Ij|。

固定负载点，在正贡献中取最小半径且并列最后下标k，令其半径1。捕获项≤tau，自项≤a^(−n)。其余未捕获项只有1≤r<t。令alpha=a/(a+1), beta=1+1/t, c=(a beta−1)/(a+1)>1, z=y/(a+1)，并逐坐标反射使末中心u≥0。实际未捕获测度

 dmu=alpha^(−n)2^(−n) sum rj^(−n)1_{Ij/(a+1)} dz

总质量恰等于未捕获权重。由rj≥m(z)=max(1,||z||∞)、不交性和原帽包含，mu被lambda_corner=alpha^(−n)2^(−n)m^(−n)1_Hcorner dz支配。这里有比积分单调更强的集合事实：c>1使正向zi−ui<c m恒成立，故Hu={存在i,zi≤ui−c m}包含于H_{alpha·1}。无需假定末中心实际在角点。

取1/t<sigma<1，后序区A=(-sigma,sigma)^n∩{min_i zi<−1/t}，前序区P_R=Hcorner∩{m>R}，R>t sigma。A全部位于corner帽内，且m=1、||z||≤sigma；所有面并集只算一次，归一体积为

D_n(t,sigma)=sigma^n−((sigma+1/t)/2)^n>0。

若mu(A)=0，总帽质量立刻损失lambda(A)。否则A中存在一个真实未捕获项j和真实增量点w。未捕获给rj≤||u−w||/c≤(alpha+sigma)/c。固定R>t sigma后，此上界在足够大a时小于R。P_R中所有实际项i有ri≥m(z)>R≥rj，故i<j。中心与增量的位移误差至多(ri+rj)/(a+1)；只要

1+sigma/R+2/(a+1)<alpha beta，

每个这样的前序项均捕获同一个yj，实际资格条件给mu(P_R)≤tau。因此未捕获总质量至多

alpha^(−n)Ja(t)−min{lambda(A),lambda(P_R)−tau}。

这是同一真实历史的一个分支论证；不把若干不同最坏位置或独立分解无代价相加。A不需位于每个Hu；实际A无质量本身就提供包络缺口。

精确极限尾量为

P_n(t,sigma)=n log(1/sigma)−(n/2)∫_{t sigma}^t ell_t(rho)^(n−1) d rho/rho，
ell_t=1+1/(2t)−1/(2rho)。

若B_n>D_n且P_n−D_n>sqrt(A_n(B_n−D_n))，先固定R>t sigma足够近，再固定tau=sqrt(A_n(B_n−D_n))；所有足够大a均满足guard及lambda(P_R)−tau>alpha^(−n)D_n。由原覆盖乘积和负载加上tau及自项a^(−n)，再取纯标量a→∞，得到

Cn≤[sqrt(A_n(t))+sqrt(B_n(t)−D_n(t,sigma))]^2。

R↓t sigma也是同一算子的标量界连续极限，不要求选择族收敛。一般L1和全部尺度桥接完整保留：紧坏集每点取原算子的真实见证半径R_x，以R_x/a小立方体作有限子覆盖，积分加权负载；固定半径平均连续，最后用内正则性，无预设水平集有限测度，也无输入支撑或半径上下限。本轮路线不依赖旧Green/S3/小时间论证。

自足有限a证明、一般sigma条件及严格guard见geometry_theory/REPORT.md §3.1；独立审计见audit/round030/uniform/geometry_review.md。七维另以a=2^512、R=t sigma+2^−80和明确有理tau直接认证有限覆盖/负载乘积小于同一U，所有数值见独立finite_n7.json。

## 真正维数一致的显式改进及其界限

令L=logn,l=logL,T=nL−n l。对每个整数维数满足L≥18，严格证明

Cn < W_n − (5^n−4^n)/(3^n t_n^n)，
W_n=T+(3l+1)n/L+(9l^2+1)n/L^2，t_n=n/(L−3l)。

W_n本身比旧T+(6l+5)n/L在共同范围L≥32严格小，并扩大了显式适用范围。证明用Beta换元

I_n(n/x)=∫_0^((n−x)/2)(1−v/n)^(n−1)/(v+x/2)dv，

质量m=1−h^n及首矩r<n/(n+1)<1，Jensen给I_n>m^2/[(x/2)m+r]。取x=L−3l，h^n≤1/n、x+2≤L，得I_n≥2(1−2/n)/(x+2)。再合并A≤n/L^3、B≤nL、交叉项≤2n/L和显式对数余项，得到W_n。选择sigma=5/(3t_n),R=7/4时，新机制D=(5^n−4^n)/(3^n t_n^n)；B≥T−4n/L>8n、tau_new≤n/L、D<1使尾条件自动严格成立。全部阈值在L≥18有解析证明和Fraction补充检查，详parameter_theory/proof.md及audit/round030/uniform/proof.md。

必须说明障碍：本单一unitcore-cut族每个可行点E=B−D满足E>A（因P≤B且P−D>sqrt(AE)），故旧新成本之差

D[1+2sqrt(A)/(sqrt(B)+sqrt(E))] < 2D < 2。

定义新族可行下确界V_n（空域时+∞），有V_n≥U_n−2；充分大维有可行点。因此这条新机制不能降低nlogn的领先系数1。细化证明另得旧U_n=T+(3l+1)n/L+O(nl^2/L^2)，任意旧族极小参数n/t_n^*=L−3l+(7l−2)/L+O(l^2/L^2)，上述V_n也具有相同精细渐近。这里的模型下界不是Cn下界，也不排除外壳后序区、不同cover策略或joint覆盖/负载新目标。持续目标未完成，下一轮转向新族参数优化和更强共同历史约束。

## 准确比较、公开与内部

原七维已证上界为4910992105193/200000000000，原完整成本下端在upper/round028/cap_candidate/result.json。新U不仅小于该显示分数，而且小于原完整Arb下端；同t新旧完整有理区间也严格分离。内部R029的n4/5/6/8/16/64成本使用revision_v3/numeric/table.json完整端点比较，不从显示小数重建。

公开比较范围沿用本日已核原始文献，未新声称穷尽：Stein–Strömberg一般范数球O(nlogn)、Aldaz固定eps充分大维(2+3eps)nlogn、测度类2^n的Lebesgue特例，及二维专门结果的版本范围，均见revision_v3/literature。本轮保持全尺度方体模型，未用欧氏球或稀疏尺度结果替代。新结果是本项目内部严格推导与证书，不称正式公开发表、全球首次或真实Cn最优。

## 实际运行与失败诊断

1. P30A：六维局部驻点/余项诊断，六进程，每维41–55次，属于新诊断而非全域结论；同族七维微调仅约3.87e−14，未当作主机制。实际资源与PID在compute/results.json。
2. G30A：512-bit Arb有限根隔离，n4为159盒，n7为286盒，两者零残余并各恰两根。与独立端点/尾部解析证明合并，旧Fcap族在n4和n7有唯一global极小；此唯一性不属于新二参数族。n64达到4096盒，有4041未判定，停止且不盲增精，不称完整。独立根审计另行栈遍历与有理叶覆盖重跑通过，见audit/round030/roots。
3. L30A首版把ell常数和beta−1/rho=2ell−1写错，遗漏+3ell^5，原排除结论无效并已撤回。修正版正积分、SymPy展开与有理证书翻转为L(7/4)>tau；主MATLAB三个指定点与修正版吻合。失效原件仅在独立clone外保留，其hash及错误在withdrawn_diagnostics.json；不能引用为失败路线的严格排除。
4. N30B：固定旧t、七维并行各55次sigma浮点探索；每维明确有理点随后准确认证。原搜索主动退回边界1%，不是求得全局最优。全部七点的B、P、D、成本和有限a均通过；新n7另经不导入producer helper的独立Fraction/500项atanh尾项/isqrt认证。主整合后实际执行标准库verify_accepted.py再次通过七点及有限a，没有新搜索。每批实际PID、CPU、wall、RSS保留；G30A的45秒观察上限未作为持续目标总时限，后续只按有限批次数量终止，不设程序总运行时限。
5. 本机15核/24GiB；已实际使用6进程、7进程数学批次，下一轮登记12进程新公式优化。没有空转或填充内存。MATLAB仅主串行诊断，无clear或path重置。
6. soft核三类同核预算已证明不优于旧U_n，有限a非共同角点存在准确反例，scope见soft_theory/proof.md，不能将障碍扩为所有覆盖算法。clone部分缺失无关Git blob导致子代理局部提交失败，主只复制其三份必要成果；不修复或修改无关项目，也不将该操作失败当数学失败。

最小复核：在项目根执行 python3 audit/round030/numeric/independent_n7.py 和 python3 audit/round030/numeric/finite_n7.py；只需标准库。七点整表可执行 python3 upper/round030/verify_accepted.py，同样标准库（原helper未使用的mpmath导入已移除，算术函数AST全部不变，见code_adaptation.json）。根完整性重放需要python-flint；候选探索才需SciPy/mpmath。一般维数数学依赖解析证明，有限样本不能替代它。

所有来源、必要代码和完整端点将由本目录及audit目录MANIFEST锁定。本轮结束但持续goal保持active，下界仍暂停。R029论文仍为有效较弱旧稿，尚未吸收本轮新cut；本轮交付为研究报告、准确账本与最小证据，不伪称已重编论文。


## 最新调度

用户已于2026-10-28明确要求聚焦一般高维覆盖改进，停止七维特异性改进。本轮七维点和已完成的R031数值只是此前计算的封存，不追加新七维实验。持续目标转向一般多尺度覆盖、联合cover/load及高维外壳；原下界不启动。
