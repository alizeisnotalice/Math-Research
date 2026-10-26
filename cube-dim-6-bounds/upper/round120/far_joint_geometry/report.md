# 120：真实远区的先时间积分与椭球角预算

本轮已实际完成两个固定配置，随后一次保存证书核验PASS，现封存停止。没有重跑037 pair envelope、旧时间格、116表、MC或native。

## 合法全空间与真实far接口

原body32/common-time Cauchy模型、Lebesgue、q=6/5保持。远区严格使用059的 F={x=(1,y),0≤y≤1,dist(x,S2)≥1/64}；S2包括全部496实际unordered body-pair planes，与已付middle/CEH分配不重叠。所有pair交叠保留，没有只选最近plane，也没有把31 unsigned距离等式用于整个有符号空间。

先用106正确正Gaussian混合完成t∈(0,∞)积分，再做det加权Jensen。对同一个真实far的整个cube边界，QΣ(x)^−3为临界−6次齐次；精确椭球角积分与384对称因素得到

    actual I_F ≤ J_F,
    J_F^5 = (3π³ p1(0)/2^53) · P_F^5,
    P_F = Pr[X1/||X1||∞ 落在完整有符号far边界]。

这是先时间、后空间的上预算，未对每个t取Lq再Minkowski。完整证明在proof.md，包含Gaussian协方差2t²Σ、π³角归一化、各个极限、真实S2及一般L1/全部正共同尺度接口。扩展非负不等式允许角奇异点处为∞，但完整积分由概率归一化支付，没有删除奇异层。

## 两项固定配置的严格结果

配置1固定R=256。六个真实Cauchy(32)边缘给maxnorm尾上包；全部496个真实rank4投影用108峰和四球体积给tube事件上包。union只用于估计排除事件，未假设pair无交叠或坐标独立。准确得到

    P_F > 1−3/(2π)−1519/(1152π²) >3665/10368。

完整保存π端点在同一R上给更锐P_F下表达式约0.3889354502792523。没有扫描R。

配置2以ψ(ξ)≤32|ξ|所得 p1(0)≥15/(2^33π³)作解析系数见证，严格认证

    J_F^5 > (45/2^86)(3665/10368)^5。

因此 J_F>5.02696160042959e−6；它超过111原始完整较宽余款R_U≈4.81959607749471e−6约1.043025498倍。使用同R、更锐π表达式的下界约5.53100571107247e−6。所有比较使用FULL111有理端点及精确五次方残差，不使用这些显示小数。

该严格障碍即使改用exact p1(0)与exact P_F仍成立。**被排除的是本轮指定的det加权Jensen预算族，不是actual I_F、原模型或C6。** 它丢失的是不同Gaussian形状在G^q内部的空间补偿；单一Jensen权重后的临界空间积分只剩P_F。它也不是对其他保留形状关联的估计的反证。

108全6D真实峰同时给一个有效、完整、但不竞争的actual far上包：

    I_F ≤ (81/(2^80π³))^(1/5) <1.84895844055132e−5。

该数上包的完整有理端点保存于results.json与certificate.json。它支付全部时间和全部actual F的有限粗上界，但不能支付111现有竞争额度R_L≈4.53710531716404e−6。没有将该完整far粗上包登记为竞争成功或新C6结果。

## 实际执行、核验与封存

主完整数学/source静审release SHA6a3e8dc793ebb4d0fc2614e55fa85e7233122329ea8f36fa016049dc0e4f5aa3之后，才实际唯一运行observe_once.py。tool chunkae3d62、exit0，科学PID71808、observer71800；两个配置全部完成。完整科学子进程从spawn前至最终收据、解释器退出及reap后，wall0.24142566599766724s，CPU0.22662s。observer自己的收据写入/退出不包括在该child范围。无运行时间上限，最多一科学CPU，不使用117的native资源。

唯一saved-only checker工具chunk8964f2、exit0、PASS_SAVED_ONLY。它不导入生产生成器，0新科学配置，核验精确系数、fifth-power残差、根端点及结论范围。复核入口：

    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B check_saved.py

本机数学venv实际使用Fraction、Sympy与Arb384；17项来源pins结束时全部未变。原098/104/107/111/116均未改，没有新方向、代理、下界算子研究、MATLAB、stage或push。现仅完成这一指定族诊断。全域g0、竞争far额度、完整S仍未闭合；下一轮CEH回收需主另行调度。
