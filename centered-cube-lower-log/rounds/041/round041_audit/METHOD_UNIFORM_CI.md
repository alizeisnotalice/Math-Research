# 固定target下对合法pointwise方法同时有效的KL区间

结论：父线提出的强化**成立**。它不把事后分析变成独立确认，但在固定target、固定真实cap及逐点夹逼前提下，方法选择本身无需额外union。原1008预算仍可保守沿用，不修改40/41冻结结果。

先已运行的256个有限样本枚举验证CI单调支配；另实际运行 `method_uniform_ci_checks.py/json`：1152个邻点单调检查、512个全样本自适应夹逼检查，以及data-tuned cap/target的有限反例，均通过。一般结论由下述自含证明给出。

## 1. 有界真实贡献的KL两端点

固定一个target。真实贡献Y_i iid，0≤Y_i≤C，**C>0是对总体分布有效的固定界**；µ=E Y。令x为归一化样本均值barY/C，β=log(1/α)/N。定义

`L_C(barY)=C inf{p∈[0,1]: kl(x||p)≤β}`，

`U_C(barY)=C sup{p∈[0,1]: kl(x||p)≤β}`。

其中kl(x||p)=x log(x/p)+(1−x)log((1−x)/(1−p))，端点按连续极限。x0时L0、U=C(1−e^-β)；x1时L=Ce^-β、U=C。

对Z=Y/C∈[0,1]，凸性给 `E e^{tZ}≤1−m+me^t`，m=E Z，任意实t。独立乘积、指数Markov及对t优化，得到两侧Chernoff界exp[−N kl(a||m)]。反演即

`P(L_C(barY)>µ)≤α`，`P(U_C(barY)<µ)≤α`。

这只要求Y_i iid有界，不要求Bernoulli；没有normal近似。

## 2. 端点随样本均值单调

内部根kl(x||p)=β处，

`∂_x kl=log[x(1−p)/(p(1−x))]`，

`∂_p kl=(p−x)/[p(1−p)]`。

lower根p<x时两导数分别正、负；upper根p>x时分别负、正。因此两者的隐函数导数dp/dx=−(∂_x kl)/(∂_p kl)>0。端点连续延拓后，L_C与U_C在全部[0,C]上均不减。

## 3. 任意合法方法的pathwise支配

同一批数据上，对每种方法m允许产生数据依赖的数字l_i^(m)、u_i^(m)，甚至每个数字依赖整个样本。**不要求这些bound数字iid。**只要求对每一种可能数据实现，

`0≤l_i^(m)≤Y_i≤u_i^(m)≤C`。

若原bound在[0,C]外，可在已知真实Y∈[0,C]的前提下合法裁剪。由均值及端点单调性，每个数据实现同时有

\[
L_C(\overline l^{(m)})\le L_C(\overline Y),\qquad
U_C(\overline u^{(m)})\ge U_C(\overline Y). \tag{1}
\]

于是同一个真实Y的两个浓缩事件，**同时罩住所有合法pointwise方法**，包括观察数据后选择方法、对合法upper取min、对合法lower取max、不同样本点选择不同合法方法。无需对方法数再union；也不需将选定bound函数的样本当iid。此处覆盖的是固定µ=E Y，不是自适应选定bound函数自身的期望。

`l≤u`的观测断言远远不够：必须证明它们夹住未知真实Y。通常应由所有counts/隐藏位置上有效的解析包络或正确的认证程序保证；仅在已见样本上检查数值彼此次序不能建立该前提。

## 4. 多target与当前预算

对预先固定的target族t，分别有固定C_t和各自iid真实贡献流；不同target可共享样本、相互相关。按每个target-direction分配α_t,L、α_t,U，union费用只需Σ_t(α_t,L+α_t,U)。在这些真实target事件上，(1)同时涵盖每个target的全部合法方法。

本轮252个固定target=4R×3model×21h，故504个target-direction事件已足够形成方法统一的数学陈述。继续使用既定α=.05/1008当然保守；无需为本lemma改变已注册1008预算或重算其端点。方法预览、MCD提出时间、replay与旧样本使用仍必须如实披露；该结果不是新的独立确认，也不提供“条件于研究者为何报告这个项目”的选择条件coverage。

改变source/model/h到原固定族以外，是换target，不能归入同target方法改进。有限反例：一次观察U~Uniform{1,…,100}后，挑target Y=1{U等于所见类别}；每个固定indicator均值.01，但被挑target的N1样本恒1，α=.05的lower为.05，失败概率1。未支付目标选择的union会失效。

data-tuned cap也不自动有效。N1、Y~Bernoulli(.5)，若选C=max(观测Y,10^-6)，在Y0事件中upper仅.95×10^-6<.5，失败概率至少.5，远大于α=.05。总体有效的固定source-only cap不能替换为观测最大值。

## 5. 对当前41分析的准确口径

固定真实贡献可取 `Y=λ D(N)^−1 P_hidden(max_{L∈J} A_LF>λ |counts=N)`。原Markov/source帽必须也上界这个真实Y；原/新conditional bounds逐点夹住同一tail时，旧lower、variance/3 lower、旧upper、refined-Cantelli、McD以及它们的合法min/max都属于同target方法。

完整FP、source均值/BE包、实际提议RN、underflow及尾函数计算仍未获完整向外认证；因此目前保存数值不能据本lemma升级为严格有限真概率或Cn证书。**合法lower下溢保留0可保守；upper的正数下溢为0并不自动保守**，要认证必须向上处理并传播误差。高精度重算CI只能验证统计反演算术，不能修复上游pointwise夹逼前提。

该lemma解决同一固定target上合法包络方法的事后选择费用，不解决换target、未证代理、随机cap或完整计算认证。
