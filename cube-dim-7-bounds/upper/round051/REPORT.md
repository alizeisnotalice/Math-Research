# 第051轮：灵活soft预算族及更强固定tuple的原算子定理

完整generic族已独审ACCEPT；不是单纯scalar优化。固定a>3、s>2a、0<d≤1/16、d<1/(32a)、sigma=s d²<1/2，epsJ/epsA/epsD在(0,1)、epsS/epsT>0，且下述bracket>0。每个完整tuple都独立于n，其共同充分大N可依赖全部参数；不把n自适应搜索当定理。

取P0=3579/4096，qmin=1−3d/4，u_lo=P0/[2(1+epsJ)]，u_hi=1/[P0(1−d/2)(1−epsJ)]，m*=(1−epsA)d/[2(1+3d/4)]，M*=max(4,2u_hi/qmin)，b*=(1−3/a)u_lo−epsS−epsT>0，Gamma=m*b*/M*，gamma=(1−epsD)Gamma，sigma=s d²，c_raw=gamma−sigma/(1−sigma)。完整结论Cn≤F(A_n(t),B_n(t)−gamma n)，t=n/[(1−sigma)L]；若c_raw>0，每个固定0<c<c_raw给Cn<T−cn eventually。未声称c=c_raw或连续全域最优。

几何与真实mu不改，只取theta=e^-aLd²、h=2^-20/(Ld)、b=min(1,theta/G)。原exactJ、conditional pin/free log矩、a.e导数与AC全q桥仅需d≤1/16。典型logproduct上界517/4096给准确P0≤∏J≤1/P0；另J→1在所有matched w/q上逐点统一，给任意固定epsJ，而非由product界反推。好目标概率趋1给任意epsA的质量m*n。全部source pin signs求和保留，得到U∈[u_lo*n,u_hi*n]。

capture/G逆行概率下E(−logG)≤3Ld²，故soft保留≥1−3/a；匹配log宽[h,2h/qmin]给raw row界。sourcecap/k为O((d/h)e^{−(d/32−ad²)L})，严格d<1/(32a)使其≤epsS；tau/k的指数(s/2−a)d²>0使其≤epsT。最终maxdegree≤M*k，energy同一缺额给Gamma*n，留epsD后减gamma*n。ca²/ca有限a桥、b≤1同一真实eligibility、一般L1 witness及全正尺度均完整重证，未用实际mu的概率分布。

完整scalar正代价用−log(1−sigma)≤sigma/(1−sigma)，交叉根式与A为o(n)，因此gamma和sigma必须一起付款，不能只优化负gamma。epsJ/A/D缺少上界或b*非正都可能制造双负假阳性，本轮已将它们写成必要域条件。

R052提出的一个准确tuplea73/8,s293/16,d1023/2^20、epsJ/S/T/D1/1024、epsA1/64，已由另一独立Fraction核全部域并给c_raw=6468305130094850733787980147/370780099533215320770594052505600>2^-16。因此本轮full族接受后，该tuple确激活eventual Cn<T−2^-16n。这个实例只需准确可行性，不需要R052全网格最优或top64完整性。它在R052运行当时是conditional，原操作记录身份不重写。

本轮source尚用draft审前字样冻结，当前接受状态由最终review负责。作者只作P0算术小核，独审实际one-tuple Fraction核而未重跑7.37m枚举。主独立重建上述完整预算及same-history桥，保留证明责任。源码、审稿、单tuple准确receipt及hash均归档；没有新n7/下界/MATLAB/旧网格或公开发表声明。显式域由R053/R055另轮量化，不能在本轮预设。
