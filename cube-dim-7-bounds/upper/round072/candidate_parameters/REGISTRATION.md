# R072 · 固定L4096候选及有限支付预算收紧（非operator证书）

唯一基础几何：s=3/2,qmin=2/3,qmax=91/100,sigma=19/200,chi=9/200,L0=4096,delta=2^-19,rlo=1+delta,rhi=273/200-delta,h=(1/512)/L。eBaseline=2^-16,eJ=2^-20,eN=2^-19,eRad=2^-20,epsA=1/65536,epsD=2^-20,epsRet=1/3,mL=311/1000。

首核给定epsS=3/40,epsT=1/4096,gamma=19/100,c=9/100。只试固定有限支付表：epsS∈{75,74,73,72}/1000，epsT∈{1/4096,1/5000,1/6000}，共12行；几何、h、profile参数和正态输入不变，禁止以未获证的H2域假设改变几何。每行gamma取低于(1-epsD)Gamma的1/100000有理格点，若恰落格点则再下退一格；候选展示c=gamma-1/10，实际正价格上包为999/10000，留1/10000余量支付完整scalar。基础原gamma和c单独核验。

新算术：Fraction正项exp Taylor前缀S64作cap/tau全域上包；log(200/181)用7项odd alternating upper；BE候选的mean/variance/z、pi界和normal integral下包只作给定条件下的准确数值接口核验。H2 Var区间和BE独立非同分布定理/常数目前都尚未合格，必须在报告与每条候选receipt显式保持conditional状态，不可声称Cn改善。

全尾L/S64(aL)、L^(3/2)/S64(bL)以解析log导数或系数证书证明；其余guard用明确Taylor powers。只核本有限表与当前新基点，不跑旧证书/作者程序/生产grid/MATLAB/Git，不改他人文件。任一失败原样记录，不能删行或误当通过。单进程低资源。若本表不提供严格更好的可用conditional预算，则如实报告无改进。

运行前补充：target annulus的当前新tuple log用20项even alternating lower；保留率normal积分用pi<22/7的正积分恒等式和exp(-t²/2)≥1-t²/2。roots刚传R070候选H2接口：h≤2^-20,n≥2000,x/n≤1/1000时baseline H2在[.00528,.005295]x，双边[.00527,.0053]坏率≤5760000/n；shift≤xh/5，去任意sourcepin须付1/(125rho)。最终lemma/hash未落盘前只登记为待绑定接口。六坏事件各epsA/8，不能漏第二profile；Var范围只在每captured sourceface使用，不赋给mixture。
