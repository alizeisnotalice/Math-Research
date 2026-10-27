# R035 energy_theory 独立审计

结论：**ACCEPT，共同充分大n版本。** 未认证具体数值起始门槛。没有七维搜索、生产参数计算、MATLAB或Git/网络操作。

受审文件：
`/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round035/energy_theory/PROOF.md`

受审SHA-256：`f7d98eb21b89d83264053fda2c594387387facfdcbdbec7cf83a41e5dd55f0bd`。

## 1. 已独立重建的结论

令L=logn、l=logL、d=√(l/L)，原同参数A0,B0,Wn如受审稿定义。对所有充分大n，严格成立

C_n≤F(A0,B0−nd/65536)<Wn−nd/65536，

以及

C_n≤nlogn−nloglogn−(n/131072)√(loglogn/logn)。

该扣减比R033/R034的n/√L多√l增长因子。它不与本目录`PROOF.md`中的单次容量障碍冲突：新energy机制的row degree约k≪nd，不要求每row有至少δ的前序质量。

## 2. 后序质量接口

受审稿使用full post annulus∩D∞，不保留R033仅负最大面的限制。固定每张post面，一坐标确定贡献d，其余独立。精确均值m_w=d+(n−1)d²/(2r)，方差≤nd³/(3r)。容忍t/100给坏率O(L²d³/n)=o(1)。K_w均值≤3nd，Markov坏率≤1/4。因此每面至少半数good点对全部充分大n一致成立。

原post annulus质量≥nd/4；原角点捕获的壳系数为(n/2)ell_t(r)^(n−1)，并且ell_t(r)≤3/4对充分大n均匀成立。扣去O(nd(.75)^(n−1))后，λ(A)≥nd/16成立。这里没有把全cube条件律替换固定face条件律。

## 3. Source分布及全局删除

独立核验详见`COLUMN_REVIEW.md`。准确有

Var S_z≤(2/3)Ld³，
K_z=1+Bin(n−1,2d/ρ)，P(K_z>16Ld)≤4/(121Ld)。

source annulus质量≤nd，故坏统计质量≤(20000/3)n l²/L+(4/121)n/L=o(k)，其中k=n/[100L^(3/4)√l]。原corner捕获source质量O(nd exp(−Ld/32))也为o(k)。因此可在每row直接扣除同一全局坏质量，而不需宣称good事件在捕获条件下独立。

## 4. 全域uniform product核验

反对角匹配|r−q−d|<h、κ=β−ζ令阈值shift b=e+ζρ，|b|≤h+ζt=o(d)。K_w,K_z明确记录r−2d或ρ−2d以上的全部受影响坐标。

Row归一化缺失和与nd²/(4rρ)的误差uniform limsup≤13/200：t/100项贡献1/200，12nd h项贡献3/50；其余项趋零。基准

nd²/(4rρ)=(l/4)x/(Lrq)=l/4+o(1)

因l d+l²/L→0而在整条匹配带uniform成立。最大row损失O(d/t)，平方和o(1)。故全乘积及删除任意一个坐标后的乘积均在e^(±2)L^(−1/4)之间。

Column同理，归一化缺失和误差uniform limsup≤17/200；最大单坐标损失O(d)，平方和O(dl)=o(1)。因此同一e^(±2)包围对全部good source和匹配r成立。受审稿没有遗漏固定壳面坐标，也没有把删一坐标的变化当成固定比例损失。

## 5. Degree与radial边界

源壳每坐标至少一面被捕获，至多两面。故row完整壳系数介于(n/2)e^(−2)L^(−1/4)与ne²L^(−1/4)。反对角带在I端点最短宽度仍为h；相应log长度介于h和4h。删去source坏统计与corner部分的o(k)后，所有good post row均有

e^(−2)k/4≤R(w)≤4e²k。

反向column的post log长度≤2h，捕获面系数≤ne²L^(−1/4)，故C(z)≤2e²k。restrict到A或D∞只减小上界。没有遗漏负face/fullface的factor2。

## 6. Energy支付与实际历史

A与P的空间径向支撑不交。对任何μ≤λ，Δ=λ(H)−μ(H)，联合Borel核K有

Eλλ−Eμμ≤Rmax(λ−μ)(A)+Cmax(λ−μ)(P)
 ≤max(Rmax,Cmax)Δ。

一次实际eligibility对μ-a.e.post给rowμ≤τ，故Eμμ≤τλ(A)。由row下界得到

Δ≥(rmin−τ)λ(A)/max(Rmax,Cmax)。

τ/k≤100√l L^(−1/4)→0，足以取τ≤rmin/2。这给Δ≥e^(−4)nd/512>nd/65536（e<3即可）。这一步正确保留column capacity，未把不同post的局部cut任意相加。

有限a中，在A/P上λ_a=α^(−n)λ准确成立。两侧degree同乘c_a、energy乘c_a²、mass乘c_a；除以c_a M后总缺额不小于极限测度所得常数。无需TV或历史选择收敛。

对a≥8n^4，post实际半径<3、source实际半径>4，所以源项确实先选。κ+2/(a+1)<αβ使增量点捕获传到中心捕获；每个source增量子集的μ质量不超过该项原weight，故同一post实际eligibility支付全部rowμ。该核联合Borel，Fubini合法。

## 7. 预算正性、一般L1与渐近结论

δ=nd/65536=o(n)，旧B0>8n给根号正性。A0≤n/L³、B0≤nL给τ≤n/L。有限a负载和覆盖因此对全部原有限实例合法；固定n,t,τ后只对数值上界令a→∞，一般L1/all radii桥接保持原真实立方体见证。

F(A0,B0)−F(A0,B0−δ)>δ，接Wn得到严格比较。Wn−T=O(nl/L)=o(nd)，最终负修正系数1/131072正确。

受审稿仅声称一个共同充分大n阈值；所有o(1)、o(k)均已有明确uniform估计与正固定余量，所以无需先认证logL≥256或512的具体门槛即可成立。若日后要声称具体门槛，须另外把各o项统一数值化，本次ACCEPT不覆盖那一加强。

## 8. 计算核验范围

本审计独立运行了SymPy/Fraction：column一维一二阶矩、sourcebad系数及一个保守l=256幂界；相关精确分布核在本目录`verify.py`中另有代数核验。全维数energy结论由上面的解析重建给出，未把数值均值或有限网格当成证明。
