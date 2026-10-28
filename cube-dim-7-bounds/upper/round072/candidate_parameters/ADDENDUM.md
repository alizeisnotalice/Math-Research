# R072 增量：同几何的17/25保留预算及R070真实来源绑定

本文件追加于原42共同gates/12行支付表运行之后。原REPORT.md、REGISTRATION.md、verify_candidate.py、verification.json及其MANIFEST均保留原字节，不回写当时真实的H2 pending身份。

## 一、R070来源现已落盘

实际全文读取的新LOCAL源为
`/Users/zhengzhihao/Desktop/T/cube-r030-roots-audit/cube-dim-7-bounds/upper/round070/second_profile/PROOF.md`，SHA256
`4d2f4fc2d0f3353ba82e32487dba68804d3021d4c0d477137c11f6c1a004e2b8`。

该源§3–7给出原报告登记的H2 baseline mean/bad率、强AC传输xh/5、去任意source pin的B_j≤1/(125rho)及真实中心化variance范围。本tuple的h=2^(-21)≤2^(-20)、n≥N≥2000、x/n≤1/1000、eJ≤1/1024及profile坏率≤epsA/8均已在原运行明确检查，不能只复用该LOCAL源较粗的epsA分配。此绑定不是把LOCAL升级为完整operator或BE证明；LOCAL仍由独立审查者负责复审。

## 二、一个额外登记的支付点

原条件CDF/BE计算已经给出逐face保留率严格大于

85309/125000=0.682472>17/25。

因此在相同、仍待核准的BE定理假设下，可以把epsRet预算改为8/25，同时保持所有几何、h、两个profile/Var输入和S64支付函数不变。这个新单点先登记后实际运行，没有重算前缀或重跑12表。

固定

epsRet=8/25，epsS=37/500，epsT=1/5000，
gamma=39/200，安全display c=19/200。

准确新预算为

Gamma=236838078735008969160929918019/
1213671334958326123751014400000，

(1−epsD)Gamma−gamma
 =7211793325439359917802882614597/
50905065349010470941533747019776000>0。

cap/S64<37/500、tau/S64<1/5000仍是原回执的真实上包。完整root预算保留价格<999/10000，并有

gamma−999/10000−19/200=1/10000>R8(L0)。

8项新Fraction检查实际PASS，PID1010，见retention_verification.json；该回执还实际锁定了上述R070源。所有尾部沿原报告的确定性单调界保持，新保留率条件也沿统一CDF下包保持。

所以本轮比基础c=.09更简洁且更好的**条件**候选是gamma=.195、c=.095，即39/200、19/200。它仅收紧同几何/同soft weight的分析支付预算，不需要外推新的几何/H2域。较大的gamma会对应更低的候选完整cost F(A,B−gamma*n)，但在BE/完整新机制未获证前不得称为Cn改善。

## 三、仍然停止在条件身份

本次绑定解决了“无真实H2源hash”的归档问题；它没有解决合格独立数组BE定理与C≤1常数来源，更不替代新tuple的完整actual-history/generalL1/allpositive-radius算子证明。无论本单点的精确余量多大，目前均不能计入Cn账本。原12行有6个cap失败照样保留，不被这个增量覆盖或隐藏。
