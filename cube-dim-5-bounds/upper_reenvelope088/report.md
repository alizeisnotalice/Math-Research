# upper_reenvelope088：56139节点零H新均值重包络

**全080+082既有数值层级的56139个节点已全部用085新均值界重包络，零个新H、零个新anchor、零个DD/stencil。** geometry、全部parent/child边及55819个physical frontier leaves保持不变；下端逐字保持，未重包节点始终支付旧self上端。

完整正面F=[0,1]^4上的新严格证书可简写为

\[
\boxed{\int_F H_{1/512}^{5/4}\le 66783889/1000000000000,\qquad
 I=\int_F G^{5/4}\le 70301737/1000000000000.}
\]

内部细上端约UH=6.678388883850708e−5、I=7.030173690972618e−5，准确有理端点保存在summary与各block receipt。**新I上包严格优于原完整Icap71555899/10^12（约改善1.752703%），但仍未通过正式C5<32门。**

g=1/10000本阶段使用主接受的全域证书，不再标为conditional；g105未用。新cost上端约52.56280830284697，准确测试5^5 B_U^5 < 32^4 g^5为false，所以正式C5上界仍32。这是上包未过门，不能推断实际C5≥32。也未达到I≤1/25000。

## 1. 冻结输入与执行范围

本轮仅新增upper_reenvelope088，不改079/080/082/085或主clone。085 manifest994becf91e7d5d5c48817dc1dbfcdf0e551a12b84a3aac0326f470a325fe4102，report23f8d31bc7a56f20c347f0979e6c05806f135fbd99a41e1fd25fb79ffacdc7e5，moment967b1e8a93d161a946476d06111b1066d8c620b02a89c829f48a60cf21ceba94，源manifest每项实际核hash。moment085.py是相同SHA的最小实现快照。

082 manifest1f63457e9452ccd7654fc59c636d5bb6980efbb70cefc9e8234f550f0115ae97、080 manifest28eeeaa38ae47399413ab76f16630ffb43318d77dde3bb53bd2ee1dffc8ac55b、079 manifest28e42068c0619a0769ed2130966c29d4f99860222a436d8298cd622c24cab3db全部保持字节锁定。复制的082初始state SHA为1ab46299468867980df34281ec20403881b3f0ea6acc70589c7d0e0b6545a034，只在088 scratch中更新复制品。

56139节点包含52360个080 L5 roots和3779个082新节点。所有原anchor、真实fullbox、73时间端点、eps=1/512、q=5/4、body16未单位化方向及逆Fourier因子(2π)^−5相同。nodes.json.gz冻结每个node的旧selfL/selfU、level/indices/weight、完整旧raw文件路径及SHA。物理leaves按旧weighted gap降序先重包，历史320个internal node最后处理，完整顺序在生产前锁定。源节点没有被选择性删除。

主round012/globalg_acceptance.json SHA为0b5251bbad3fb507e1965d3ee7eb3fbab6526f7019880aed34a907a9d5c6beae，实际记录g=1/10000、6590个新独立分量严格包含、52360 roots与97960 leaves闭合、083/084v2修复被主接受；其原字节快照为globalg_acceptance012.json。该接受与新积分证书使用同一固定模型，g105不属于本输入。

## 2. 真正的零H均值消费

085已证明所有面内u1=0的方向都有

\[
 D_u^2p_t(x)\le\frac{25}{24t^2}\|u\|_\infty^2p_t(x).
\]

真实box和anchor首坐标均1，因此位移u1=0。085 Fraction接口从真实offset计算M2=E||u||∞²、μ、旧η=E(uᵀQu)及体积。新的Taylor均值二阶系数是(25/48)M2；其对所有盒不劣于旧(21/8)η。对称四尾盒新/旧系数比25/56仅是曲率项比较，**没有把它当作真实I误差或全局积分改善比例**。

每个旧raw文件先核完整SHA，再只读取其73个H球、几何与旧范围/均值。所有相邻H差及H_1尾段严格正，并保持同anchor同band。程序没有H评估函数，不调用任何热核生产器、DD或H对数和。ratio R没有单独保存在原证据中时，从同一已保存h∞、BQ²重新向外计算075 min(expLip,oldasinh)，仍为384-bit Arb；这是合法ratio计算，不是重算H。

对a为band左端、RU为合法全box-to-anchor线段比上端，采用

  C=min(RU,1+||μ||∞/a+(25/48)M2 RU/a²,old072 Taylor factor)。

Σm_upper C给新均值上端，与保存的旧mean及range U取min，向上量化2^−160。保存的range L、U完全不变，再在该区间用z^(5/4)的递增凸弦支付真实mean。每个node存selfU_new=min(oldselfU,weight·new chord upper)，selfL不改。56135个self上端严格改善，4个保留原上端，全部结果仍保存。

## 3. 数值层级、完整空间账本与g配对

不创建children、不动physical frontier。历史self约束和完整children约束仍界同一parent full orbit：effectiveU=min(selfU,Σchild effectiveU)，effectiveL=max(selfL,Σchild effectiveL)。每个有限block只更新已完成node；其余节点继续用旧selfU，故全部14个prefix都是完整F证书，UH单调不增。

所有selfL保持原值，因此全域L_H恰保持

  38079990466813616787770097274397379916416525 / 1461501637330902918203684832716283019655932542976 ≈2.605538679816875e−5。

这仍是截时H积分下证，不是实际I的双边数值评估。weighted physical体积准确为1，56139节点与55819 leaves及全部320展开parent的边/守恒不变。

使用完整076范数上端E（同源同eps），B=UH^(4/5)+E，I≤B_U^(5/4)；没有064 tube加项，也没有把E当作I本身相加。与合法全域g=1/10000配对，cost≤(5 B_U/g)^(5/4)，严格32门用准确有理数5^5 B_U^5 <32^4g^5直接比较，不用显示小数作门。

|完整阶段|UH上端显示|I上端显示|g100配对cost上端显示|32门|
|---|---:|---:|---:|---|
|082输入|7.188416690829e−5|7.545313351352e−5|56.41437561938|未过|
|088第8块|6.785488858825e−5|7.138372001288e−5|53.37177936015|未过|
|088全56139|6.678388883851e−5|7.030173690973e−5|52.56280830285|未过|

原独立完整Icap是71555899/10^12，新严格简端点70301737/10^12更小。精确门记录在每块receipt及summary；正式C5仍取旧32，不把较大的cost上端宣称为C5进展。

剩余physical最高weighted gap约3.17480987652e-08，node为5:31,31,17,17；top256仅作为remaining_priority诊断保存，没有据此再启动空间refine。

## 4. 有限程序、实际资源与完整独立读回

14个独立程序，前13块各4096node，最后2891。每程序在派worker前冻结parent state/receipt SHA、确切prefix及代码/输入SHA。12worker逐完成node原子保存，有匹配source/fingerprint的成功输出直接复用；pool每程序close/join并退出。新程序没有wall alarm/timeout，没有无限controller。失败时不推进完整prefix，未处理node保持旧bound；本轮没有故障或续算，不声称测试过恢复。

生产receipt wall合计34.523896秒，单块1.749784–2.853989秒。receipt从程序入口输入加载前计到pool join，最后cert/state落盘在其后，未把此口径混成整链calendar wall。首主PID74173（workers74185–74196），末主PID74799（workers74808–74819）；实际PID和时间逐块存receipt。全部worker现在均不存在。

**全量独立读回**主PID74961、workers74971–74982，实测34.672370秒，全部join。M2改用“条件于唯一最大绝对值坐标”的准确积分，把该坐标x²乘其余三坐标CDF再积分；与085 survival-CDF四因子积分是不同路径，56139个moment逐一准确相同。另独立从实际offset重建Q角点、协方差η、Chebyshev oldR/半exp平方、band mean及另一凸弦表达式。全部4,098,147个旧band质量严格正，全部量化上端通过。

读回另用递归层级（生产为level逆序迭代）逐个重建14个全域prefix，核所有self lower不动、未重包节点付旧self、UH单调不增、physical frontier/体积1及g32准确门。零新H。只读可携读回又用约2.000372秒检查56139个量化凸弦和全部14prefix层级，不读取本机旧raw缓存。

## 5. 精简证据与复现

global registration SHA4e36eba85bb16b9e9c6cdc373c89cfd2c133f5d0bc69a588f2b78f72b239ecbd；engine SHA0f6e799fc0e1ca14167a7c6cab1b4360d8348c4c72623eaf819a049e455540c8；完整独立verify SHA0473b03610e71343842863321adfd3934ffe093887f86a6327321d45aac0361e。所有source/块/cert及摘要hash以manifest为准。

`python -B verify_portable.py`默认只读，需python-flint及相邻冻结079/080/082证据包，不读取旧1.5GiB raw。作者本次--record结果在portable_verification.json。

`python -B replay_one.py --rank R --output NEW_DIRECTORY`可在新目录纯重mean回放一个已登记node；它需要冻结旧raw endpoint文件，仍不调用H。若迁移机器，应映射source_evidence_path到同SHA旧raw，而不能从不同anchor或eps的H替代。本轮没有为了测试入口重跑已成功node。

本包保存全部14个有限登记/receipt/紧凑新mean和self端点、085最小moment源码快照、节点source索引、数值internal hierarchy、g接受小快照、独立及可携读回、样本/诊断、报告与manifest。既有raw/stencil不复制进仓库。所有旧079/080/082/085文件最终逐manifest hash复核不变。没有MATLAB、提交、暂存或推送；未消费g105，未增加空间点。下一阶段由主决定。
