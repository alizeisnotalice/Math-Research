# 第30轮正向尾凸性余项：独立有限预算

核点诊断之外，根代理负责这个新增解析下界的有限算术验证。执行前冻结：alpha=47/40，R=33/20，beta=47/40,6/5,5/4，r=1/4,1,8/5，共9条。每条只计算尾补偿，不重算完整生成元/旧积分；Arb192bit包围和mpmath70位另表达交叉，后者不是证明。2项SymPy恒等式（f''与D0展开）。最多9条Arb+9条mpmath+2项符号，60秒hard，无LP/积分/采样加密，失败不扩预算。beta=alpha控制直接准确0，无灾难消去。

D0=r^(-beta)*[(p−1)−p*z+z^p]，z=(r/R)^alpha，p=beta/alpha。新增严格下界T=4D0/[alpha*(R−r)^alpha]。只保留两坐标向外t>=R−r跳的非负Bregman余项，其余舍弃，故全A>=chi*A_alpha+T，inside r<R；outside原chi界。另曲率下界Tquad将括号换成p(p−1)(1−z)^2/2，消去更稳定。

输出T/Tquad向外区间、T>=Tquad>0（beta>alpha）与mpmath直接D0表达。有限值不证明全域，通用解析证明见tail_bound_review.md；不据此宣称任何新C2界或LP成本收益。单一进程，现存输出拒绝重跑。
