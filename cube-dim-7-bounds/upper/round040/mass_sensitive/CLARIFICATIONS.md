# R040补记：f-weighting仍有树形index-overcount

本附记补充原`PROOF.md`第4.1节接口的范围；不改其合法性，不把原接口悄悄替换成一个已改善的uniform证书。原数学源hash为`4d6a22fcfb343cbdb851cea6d017c799e8b04ece5d9b9f0813cec76afba530b1`，两处文字scope修订另见`editorial_adaptation.json`。

## 1. 结论

原f-weighted共同star、unit-index-cover目标S_f= min Σ_J θ_J M_rep，在合法二叉层状cube历史上准确满足

S_f/||f||1=1+D/2，

其中D为树深度。把essential maximum换成共同f的平均，并不能单独修复对tiny但正increment的unit-index过度收费。

即使每个共同star列改收其真实volume cost

c_J=Σ_{i∈J}v_i/V_rep，

但仍只允许共同star链列，同一树上仍有S_{f,vol}/||f||1≥a^(−n)(1+D/2)。树的containment-root列则准确修复到a^(−n)Σ_{k=0}^D(2·8^(−n))^k。以下都是证书/有限支付范围结论，没有新C_n上界或下界。

## 2. 合法cube树与单一真实f

使用star_lp已证明的具体树：节点σ深k，star半径R_σ=8^(−k)，root中心0，子中心y_{σ±}=y_σ±(R_σ/2)e_1，S_σ=Q(y_σ,R_σ)，小Q_σ=Q(y_σ,R_σ/a)。每个star是共同半径cube。

取a>8/3，1≤β≤2，D a^(−n)≤τ。所有小Q两两不交，实际increment v_σ=V_σ a^(−n)>0。原b=aβ资格在深k节点只计k个祖先，值k a^(−n)≤τ，所以largest-radius算法全部选中。任何允许预算h_i≥V_i使p_i≤w_i，仍保留该合法全selected历史。

leaf stars两两不交。把一个总质量F均匀分给2^D个leaf，每个leaf内部用L1密度承载F/2^D。则每个深k star恰包含其后代leaf支撑，所以

M_σ=F2^(−k)，V_σ=V_root·8^(−nk)。

只要F>λV_root，便有

M_σ/(λV_σ)=[F/(λV_root)](8^n/2)^k>1。

因此一个共同非负L1输入严格见证全部stars。没有独立指定彼此不一致的M_i。endpoint h_i=M_i/λ也合法并使资格更容易。

## 3. f-weighted LP的精确值

不同分支stars不交，故每个共同star cluster是一个祖先链，代表为最高节点。深k代表的cost是M_rep=F2^(−k)。

下界：对深k<D节点赋非负dual质量s_i=2^(−k−1)，对leaf赋s_i=2^(−D)。若列代表深k，沿任一后代链都有

Σ_{i∈J}s_i≤Σ_{l=k}^{D−1}2^(−l−1)+2^(−D)=2^(−k)=M_rep/F。

所以任何可行unit-index覆盖满足

Σ_J θ_J M_rep/F≥Σ_i s_i=1+D/2。

该不等式覆盖所有fractional共同star列，不依赖选取某个pivot算法。

上界：令q_root=1，其余节点q_v=1/2。对每个深k节点v及每个后代leaf ℓ，使用完整v-to-ℓ链列，权重q_v2^(−(D−k))。深j节点覆盖量为

2^(−j)+(1/2)Σ_{k=1}^j2^(−(j−k))=1。

代表v的总权重为q_v，故f-weighted cost恰为

Σ_v q_v M_v/F=1+D/2。

这证明精确等式。原R040 arrangement-cell列保留代表且可扩大任意链，不改变代表cost，因此其有限列接口具有同值。

## 4. 真实volume cost与containment的区别

共同star列J包含代表j。由于本树v_j/V_j=a^(−n)，

c_J≥v_j/V_j=a^(−n)。

因此保持原unit-index共同star列时，无论如何选fractional覆盖，

S_{f,vol}/F= min Σ_J θ_J c_J M_j/F≥a^(−n)(1+D/2)。

所以仅乘c_J虽缩小了尺度，仍不完全压掉链列对分支的重复覆盖。此处不声称该下界排除所有参数下的任何uniform改善；其精确含义是：该LP在同一树上仍可比真实支付多一个深度因子。

另一方面，树的全部实际increments都包含在root star中，containment-root列J=全部节点合法。其准确cost为

c_root=U/V_root=a^(−n)A_{n,D}，A_{n,D}=Σ_{k=0}^D(2·8^(−n))^k。

M_root=F，所以一列给S_{f,contain}/F≤c_root。反向也精确：本树每个anchor j都有

M_j/V_j=(F/V_root)(8^n/2)^depth(j)≥F/V_root。

因此任何按真实Σv_i/V_j收费的列覆盖，其总f-cost至少(F/V_root)U。故有等式

S_{f,contain}/F=a^(−n)A_{n,D}。

实际selected-union比率为λU/F=[λV_root/F]a^(−n)A_{n,D}，在F↓λV_root的严格witness极限逼近该证书值。Tiny increments的实际尺度被保留，而不是每个index都被视为一个完整volume任务。

## 5. Scope与来源身份

树历史的合法深度条件D a^(−n)≤τ明确保留；不声称在固定a、原weights及τ下可任意增加D。允许a随D增大时，volume-free S_f/F可无界，但这只是证书overcount，不是C_n增长下界。

本附记的树构造与共同star LP来源：
`/Users/zhengzhihao/Desktop/T/cube-r040-star-lp/cube-dim-7-bounds/upper/round040/star_lp/PROOF.md`

SHA-256：`6aece373926c7dfc65e7373b7245d117ba3fe4ae8f17fcd09492d1acc55a2a95`。

这里给出的f-weighted精确值和volume对照是解析补记；没有运行新的toy或生产LP，也没有改写原`verification.json`的历史run身份。原报告第4.1节作为合法有限接口保持成立，但不能把“目标对远离分量加性”误写成“uniform bound已经改善”。新的positive uniform theorem仍需要额外几何、真实volume或containment控制。
