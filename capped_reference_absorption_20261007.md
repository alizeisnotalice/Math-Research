# 截帽弱参考项与主账自吸收

2026-10-07。本引理只降低所需reference合同的强度；尚未证明实际移动参考核的弱界，不改变R_angle主账。依N04技能分别登记误差和费用。

## 精确假设及证明

E是有限Lebesgue测度集合，λ,W>0。实际非负交通密度T满足0≤T≤cλ，并且T≤B+e，其中B,e≥0。假设
sup_{s>0} s|{B>s}|≤AW，  ∫_E e≤DW。
B必须是同一原来源构造的一个整体非负函数，不能仅有每个参数/来源条件下的弱界。

逐点T≤min(cλ,B)+e。Tonelli层蛋糕给
∫_E min(cλ,B)≤∫_0^{cλ}min(|E|,AW/s)ds。
令X=λ|E|/W。若cX≤A，右侧=cXW；否则右侧=AW[1+log(cX/A)]。统一有
∫_E T≤W{D+A[1+log_+(cX/A)]}。
A=0时B=0几乎处处，直接用∫T≤DW，不使用除零式。

如果此前已建立主账X≤P+C(∫_E T/W)，C>0，则
X≤U+V log_+(cX/A), U=P+CD+CA, V=CA。
利用log_+(ab)≤log_+a+log_+b和log_+t≤t，得到
log_+(cX/A)≤X/(2V)+log_+(2cC)。
所以
X≤2P+2CD+2CA[1+log_+(2cC)]。

这不是用未知弱界反向付款。只有B的独立整体弱界已证，才能使用上述结论。若C,c绝对有界，不引入额外log n、logW或输入峰高费用。E有限可由现有O(n log n)一般基线保证，不需要目标根号界。

## 对当前实际账的条件用途

当前λ|E|≤γ[4B_paid+C0𝓕_nW+C0 R_angle]，C0=8192/49，γ≤147/143。
若未来合法分解R_angle=∫_E T，0≤T≤cλ，并有同源B及误差e满足上列条件，则代入P=γ(4B_paid/W+C0𝓕_n)、C=γC0，得到
λ|E|≤2γ[4B_paid+C0𝓕_nW+C0 DW]+2γC0 A[1+log_+(2cγC0)]W。
因此目标只要求新D,A都≤√n n^{o(1)}。若分解仅对R_angle子支成立，未处理互补支必须保留，不能套全账宣布完成。

目前候选B来自按真实continuation访问mask冻结原条件门的h_L*H参考项。它没有自动继承原fullfuture；移动L、c_t及first来源混合的整体弱预算仍未证明。固定H弱型不能经条件弱范数求和替代该合同。

## 新证书

注册后执行capped_reference_absorption_guard_20261007.py：三个有限分布规模8、32、128，各80项，总240项通过。单尖峰、Pareto阶梯、两高度、嵌套高度混合；Fraction精确离散层蛋糕、atanh展开48项加有理尾得到log区间。三个C包含1、C0、(147/143)C0，四个cap独立变化。
测试不涉及n渐近，不作阶数拟合；它核验上述通用代数常数和积分方向，并不是原cube/下界模型或actual reference的空间数值证明。
