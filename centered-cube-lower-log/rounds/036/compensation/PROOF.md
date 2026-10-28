# 低频投影补偿：合法来源与剩余对角误差

## 1. 定义与来源合法性

固定K个原子a_j∈T、对称宽delta箱beta_delta及p∈(0,1)。令

    mu_delta=K^-1 sum_j beta_delta(·−a_j),
    P_delta'=K(mu_delta−1),
    z_m=sum_j exp(−2pi i m a_j) sinc(pi m delta),
    phat_m=z_m/(2pi i m),        1≤m≤M<K.

令Pi_M P_delta为含正负m≤M和零频的真实Fourier投影，Q=P_delta−Pi_M P_delta。Q均值0，Q'=P_delta'−2 Re sum_(m≤M)z_m exp(2pi i mx)。desired低频平均项及其密度导数取

    B_L=sum_(m≤M)c_m sin(pi m L)cos(2pi m x),
    D_M=sum_(m≤M)pi m c_m cos(2pi m x),

其中c_m为原来源的实系数。修正后的真实密度定义为

    rho_proj=1+(p/K)(Q'+D_M)
      =1−p+p mu_delta+(p/K) Re sum_(m≤M)(pi m c_m−2z_m) exp(2pi i mx).    (1)

这是对来源密度的改动，不是只修改Gram。其积分准确为1。若

    1−p − (p/K) sum_(m≤M)|pi m c_m−2z_m| > 0,                    (2)

则(1)对全部x严格非负，因mu_delta≥0。box峰即使相交也不影响该证明。条件(2)是充分条件；失败不等于真实密度必负。

真实cube平均准确为Y_L=1+p( H_Q,L+B_L)/(KL)，其中H_Q,L=Q(x+L/2)−Q(x−L/2)。投影由同一个来源统一完成，不随query L选择。

## 2. 严格正交与完整Gram

H_Q,L的全部非零频率在|m|>M，B_J在1≤|m|≤M，因此对所有L,J严格有

    <H_Q,L,B_J>=<B_L,H_Q,J>=0.                                  (3)

这消去HB，不等于方差已经平坦。记Q的自相关为C(t)，C0=∫Q²>0，则

    G_proj(L,L)=2(C0−C(L))+BB_M(L),
    C0=R_delta(0)−(∫P)^2−2 sum_(m≤M)|phat_m|².                 (4)

在原窗口I上，完整谱恒等式给BB_infty(L)=(gamma L²−1)/3、gamma=27/29。定义非负尾

    T(L)=BB_infty(L)−BB_M(L),   0≤T(L)≤T_M/3,

则相对原目标gamma L²/3的准确误差分解是

    G_proj−gamma L²/3=(2C0−1/3)−2C(L)−T(L).                    (5)

投影去掉原HB的同时，还去掉了正的低频高高能量8 sum|phat_m|² sin²(pi mL)。所以投影后误差可增可减，不能只将原误差表中的cross删除。

## 3. 实测C0归一化必须改真实来源

先移去primitive空间均值，C0才是不随加常数改变的物理能量。令alpha=6C0、s=alpha^-1/2。

方案A：真实来源改成rho_H=1+(p/K)(s Q'+D_M)。其峰质量变为ps、背景变为1−ps，必须另验

    1−ps−(p/K)sum|pi m c_m−2s z_m|>0.

此时Gram为s²HH_Q+BB_M，对原目标的误差准确为

    G_H−gamma L²/3=−2C(L)/alpha−T(L).                           (6)

不能继续把该来源的峰权重描述成原p。

方案B：保留峰质量p与背景1−p，真实来源改成rho_B=1+(p/K)(Q'+sqrt(alpha)D_M)。另验

    1−p−(p/K)sum|sqrt(alpha)pi m c_m−2z_m|>0.

此时desired B同时乘sqrt(alpha)，Gram为HH_Q+alpha BB_M，整个目标也变成alpha gamma L²/3，误差为

    G_B−alpha gamma L²/3=−2C(L)−alpha T(L).                     (7)

对应原参考v=n p²gamma/(3K²)的实际方差参考变成alpha v。若将来使用鞍点阈值，必须据此重新定标。方案B和方案A的相对Gram误差相同；只比较绝对误差时，方案B较小的一部分来自目标整体乘alpha，不能当作额外平坦收益。

两方案均保留质量1，因所有增减项是零积分导数/三角函数；来源正性由各自的条件保证。

## 4. 现有预算没有自动突破

对随机等格jitter，已有来源事件只给远lag的统一绝对误差约sqrt(log K/K)。投影后的C(L)为

    C(L)=R_delta(L)−(∫P)^2−2 sum_(m≤M)|phat_m|² cos(2pi mL).

低频移除不提供这个剩余远lag函数的新小帽。现有事件加三角不等式至多给原绝对误差、平滑项、均值平方和移除能量之和；不能因cross归零就把其中远lag误差归零。归一化只消去(5)中的共同常数，仍有(6)/(7)的C(L)。例如

    ||G_B−alpha gamma L²/3||_infty
        ≥2||C||_infty−alpha ||T||_infty,                       (8)

所以若剩余远lag在所选来源上仍为K^-1/2量级而谱尾较小，当前v/sqrt K费用仍在。这是本修正的失败条件，不是对所有可选来源或其他补偿的不可行性证明。

来源正性与幅度也有单独费用。由随机相位的已有二阶界，可在正概率事件上把sum|z_m|控制到O(M²/sqrt K)，故补偿密度的幅度预算约p M²/K^(3/2)。sum|phat_m|约O(M/sqrt K)，并给|H_Q,L|≤2+4sum|phat_m|。在旧beta<1/5窗口这些来源幅度费用都可小，但此处没有证明更宽窗口完整概率链。归一化还要求C0远离0；本轮具体样本对此另作有向包络检查。

## 5. 本轮具体来源的正性核算范围

NumPy先生成三份U并保存其精确二进制浮点值。positivity_certificate.py将这些U视为精确dyadic有理数，按a_j=(j+U_j)/K重新定义解析来源；不把NumPy三角函数输出当精确Fourier系数。

脚本用Machin恒等式和交错有理arctan余项包住pi；将相位化到[−pi,pi]，用49阶Taylor多项式及保守1e−35余项包住sin/cos。Decimal基本运算采用70位向下/向上舍入，sqrt结果向两侧各扩一ulp。R0、均值和box平滑损失用Fraction准确计算，得到C0与两种缩放的区间，再对(2)及两种变体的三角正背景下端作外包。

保存的36个下端全部严格正，最小大于0.0168，故这些解析定义的具体修正来源都有正密度及单位质量。该包络只支付来源正性，不支付NumPy对角曲线/误差表的浮点误差、连续L补点、尾概率或任何Cn下界。
