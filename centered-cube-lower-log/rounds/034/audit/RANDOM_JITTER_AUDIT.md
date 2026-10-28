# 等格独立jitter核：独立审查

结论：所提自含source二点核证明闭合，可解除旧field/Sidon存在性充分门槛。固定w下平滑费用乘v后是O(1/w)，不是趋0。新cubic核与现行原01中的cubic局部形状一致，匹配原Brownian强度58/(9L²)；完整概率链仍需逐项迁移。本审查没有证明新的weak或C_n下界阶。

只读round034_random_jitter/probe.py/results.json，未重跑FFT或大MC。另用本机Fraction实际核对了5个有理平移下的bridge积分，以及v100/1000/10000时的局部斜率；数值与严格代数分别列明如下。

## 1. 期望核

对x=(j+s)/K，0≤s<1，P_U(x)=1{U_j≤s}−s。因此对每个固定x均值为0，绝对值≤1，不同格的随机量独立；每条样本路径的空间均值一般不为0，这不影响以下恒等式或来源总质量。

当0≤t≤1/K、d=Kt时，只有x与x+t落在同一格的部分有非零期望乘积，其协方差为s(1−d−s)，s∈[0,1−d]。积分为

    ∫_0^(1−d) s(1−d−s) ds = (1−d)^3/6.

K≥2时结合周期和偶性，ER_U(t)=Phi_K(t)=(1−K dist_T(t,0))_+³/6。K1需要另计wrap部分，不在本轮范围；K≥10安全。实际Fraction检查d=0,1/7,1/2,6/7,1均精确相等。

## 2. 一致集中与常数

改变一个U_j只使P在该格内相差一个有符号区间指标，L1差≤1/K。将R_U−R_U'分为两个乘积差，使用两条路径均有|P|≤1，得每坐标有界差分c_j=2/K。Doob/Hoeffding给

    P(|R_U(t)−ER_U(t)|>eta) ≤ 2 exp(−K eta²/2).

在K²个固定等距周期网格上取eta=4 sqrt(log K/K)，并集失败≤2K^-6。TV(P)≤2K和BV平移不等式给R_U的Lipschitz常数≤2K；Phi_K的常数≤K/2。到网格点距离至多1/K²即足以得到补偿≤5/(2K)<3/K（取最近点可更小）。故以概率≥1−2K^-6有

    sup_t |R_U(t)−Phi_K(t)| ≤ epsilon_R
       :=4 sqrt(log K/K)+3/K.

网格补偿及并集常数方向正确。对较小K该界虽很松，但仍为合法界。

## 3. 全低频交叉项与good realization

对整数1≤m<K，Z_m=sum_j exp(−2pi i m(j+U_j)/K)。root-of-unity相消给EZ_m=0。独立性及到U=1/2处相位值的平方距离界给Var Z_m≤pi²m²/K（甚至可改善常数）。周期分布导数P'=sum delta_a−K给Phat_m=Z_m/(2pi i m)，故E|Phat_m|≤1/(2sqrt K)。

对|c_m|≤A/m²，定义单一随机帽S=4A sum_(m=1)^(K−1)|Phat_m|/m²。它同时支配所有M<K、所有这些系数及全部L,J的两个H/B交叉项绝对值，其期望≤A pi²/(3sqrt K)。取S≤2A pi²/(3sqrt K)的Markov事件，失败≤1/2；与上一事件交集概率至少1/2−2K^-6>0。probe使用4倍期望阈值更宽，也合法。这里必须先选一个good realization作固定来源，而不是把U平均后换成另一来源。

平滑的Fourier乘子为sinc，模≤1，故同一事件保留交叉帽。无需声明随机来源的所有实现都好，也无需相位/尺度交叉事件独立。

## 4. 平滑及全I的HH核

令beta_delta为对称宽delta均匀箱。P_delta=P*beta_delta，可选此周期primitive；任何加常数在H差中抵消。其自相关R_delta=R_U*beta_delta*beta_delta。对两个独立箱变量的和V，有E|V|=delta/3，故

    ||R_delta−Phi_K||_infty ≤ epsilon_R+K delta/6.

准确恒等式为

    ∫H_L,delta H_J,delta
      =2R_delta((L−J)/2)−2R_delta((L+J)/2).

L,J∈I=[10/9,11/9]时，第二个自相关参数模1距0至少1/9，故K≥10时未平滑Phi项严格为0。无需假定平滑后的该项也严格为0，统一误差已支付它。于是

    |∫H_L,delta H_J,delta−2Phi_K((L−J)/2)|
        ≤4 epsilon_R+2K delta/3.

加上全低频交叉项，真实G核相对「2Phi+BB」的全I误差至多

    4 epsilon_R + 2A pi²/(3sqrt K) + 2K delta/3.

因此v² log K/K→0足够使前两项乘v趋0；delta=1/(wKv)的最后一项乘v为2/(3w)。固定w时应写o(1)+O(1/w)，随后如需消去再取w→infty。v约n^beta、K约n^((1−beta)/2)时充分条件为beta<1/5，明显宽于旧q>H的beta<1/9。其余正背景、M<K、log幅度和概率近似资格仍须分别核对。

## 5. 与现行原01局部核及常数的核对

理想HH核为(1−tau)_+³/3，tau=K|L−J|/2。逐行复核现行01_source_kernel_detailed.tex第4/5步可见，原标签模型也是C(tau)约(1−tau)^3/6、HH约(1−tau)^3/3。在tau=r/v下，

    v[1/3−(1−r/v)^3/3] = r−r²/v+r³/(3v²) → r,

实际Fraction显示r1、v100/1000/10000时约.990033/.999000/.999900，趋于1。对角G为gamma L²/3=9L²/29，故v(1−kappa)趋于29r/(9L²)，对应sqrt(lambda)B_r−lambda r/2的lambda=58/(9L²)，与原02一致。对角HH仍为1/3，原BB的方差平坦机制也可继续使用。固定w的局部平滑形状可从新Phi卷积推得，但完整局部极限、全距离pair bounds、参数与误差一致性仍须分别迁移。

审查修正记录：初版曾误将原标签局部核记作triangular，进而提出3倍强度差的错误警告；收到主线指出后已逐行核对原01并撤销。该错误不影响上面独立推导的cubic期望、集中、Fourier交叉和平滑预算，但不能将初版的形状比较用于主线证明。

probe给K15/63/255/1023的FFT网格误差约.05135/.03604/.01635/.00973，均远低于很松的理论帽。它检查单一固定seed的离散FFT而非good事件的严格区间证书；完整m<K交叉good事件也不能只由所测最多16频率代替。其BV求积预算4K/Q保守可用，但明确未支付FFT及其他浮点舍入。因此脚本只作与公式相容的探路证据；上述存在性结论依赖解析论证。

## 6. 最终整合稿复核

本轮收尾已全文只读复核round034_random_jitter/PROOF.md。K≥2的期望恒等式、Doob差分常数、K²网格并集、BV补偿、m<K完整低频Markov事件、正概率交集、箱平滑与全I核误差均与本审查一致。对K≥10的全I式使用|L−J|而非圆周距离是合法的，因为I内差位移绝对值≤1/18；反射位移的Phi为0，但平滑后的反射误差仍已支付。

最终稿已正确写明现行原01也是cubic核，并核对lambda=58/(9L²)，不存在此前初审误报的3倍差。固定w局部卷积轮廓、O(1/w)平滑费用、beta<1/5仅为来源二点预算、未完成完整旧概率链及未证明新下界阶的量词均保留。未发现整合新增的实质数学错误。probe实现仍未获全good事件或浮点区间证书，该限制不影响解析存在性证明，但不能据probe表声称这些具体seed均已认证合格。
