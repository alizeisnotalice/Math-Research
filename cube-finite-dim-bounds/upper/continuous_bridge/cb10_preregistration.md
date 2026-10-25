# cb10：仅e1的真实源密度区间首批

2026-10-25，原独立clone，只新增cb10_*。非唯一工作者，不撤销其他修改，不改旧cb9锁源/总账/其他项目，不提交推送扩招。已读主第10轮登记、snapshot、固定源certificate和源码；按cb10_sources实际核hash（若仅ledger metadata改变，核不可变base blob及tasks完全相同，保留历史hash）。不重跑主49点/LP或其他警告点。

固定c=(4log2−1)/(24pi³)，e=e1，V={(1,±1,±1)}。仅认证

mu_c(e)=(4/pi)∫R (c−G(e+s(1,1,1)))+/s² ds。

027独审负责Jensen下界，本分支复用该源界，不重复其全域审计。用G≥5/(16pi³|x|²)认证s∈[−4/3,2/3]正部为0（含s0及s=−1/2 bodyline）；把两无穷tail转r=1/|s|：积分区间分别[0,3/2]上的(c−r²G(1+r,1,1))+和[0,3/4]上的(c−r²G(1−r,1,1))+。不在跨line附近嵌套数值求积。

内G接口先沿(a,1,1)从cb9正有理积分作准确partial fraction，给r>0的log/rational闭式；另用真实Cauchy条件密度的两个simple residue加一个double residue准确核该专用曲线，避免把尚未全独审20锥式当无条件数值黑箱。其后Arb区间直接包围log闭式，不使用scipy inner error estimate认证。

有限两级，不能无限加分区：
L1=128bit，r尾[0,2^-8]只以0≤positivepart≤c给误差，余下每dyadic block分32等份；
L2=192bit，r尾[0,2^-12]同样处理，余下每dyadic block分128等份。
两branch均从delta倍增到固定endpoint（最后短block），无adaptive新增partition。每panel用Arb区间对原始phi与导数包围：全正panel可用midpoint加C² remainder M2 h³/24；全负为0；switching不确定panel以range界或1-Lipschitz/一阶导remainder支付，不能默认为正/负。端点无限tail和全部switching误差另列，严格区间包含完整真实mu。若这些参数不足只交失败/误差预算，不扩第3级。

Python本机venv实际运行，flint Arb任意精度ball与SymPy。MATLAB只交cb10前缀便携脚本主串行；代理不调MATLAB，不clear/path/reset。最终标清true G代数身份、Arb内层log界、外层正部区间误差、无穷tail误差，浮点值与严格证据分开。首批停止于此一点严格区间或可验证接口/明确缺口，不作全源可行性或新U3/U4声称；交原独审。

同一两级运行前的解析尾加强：不改delta/分区，用专用log闭式三角界得到
∫_0^delta Qσ ≤ (8pi³)^−1{(3+delta)/(2−delta)²·delta³/3[log(2/delta)+1/3]+(2+delta)/(2−delta)·delta²/2}=Btail。
因(c−Q)+≥c−Q且≤c，tail积分包围[c delta−Btail,c delta]，不需预判tail内正部符号。此解析加强是同一无穷tail预算，不新增精度/网格；仍保留0≤positivepart≤c作基本校验。
