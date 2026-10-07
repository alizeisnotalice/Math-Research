# Critical source-tail MC：整数内部类 pooling 的只读理论审计

2026-10-07。审计者 actual_profile_budget；父任务指定 6.1-sol high。本稿只读新 [run.py](critical_source_tail_mc_20261007/run.py)、[registration](critical_source_tail_mc_20261007/registration.json)、原 [arrival module](critical_source_tail_arrival_20261007.py)、已有 [独审](critical_source_tail_mc_20261007_independent_review.md) 与 [报告](critical_source_tail_mc_20261007/report.md)，没有导入/执行 sampler、MC、oracle 或后处理，没有更改原文件。

**结论：按物理来源位置定义的整数内部类，其全部 S(y) 确实相同；给定所有抽中来源，把该类的全部16N_C个 fresh receiver scores合成二阶 U-statistic，再以 N_C/N_s 加权，与原非该类的双-replica产品组合，理想 fixed-N 模型下无偏。** 本稿还证明理想-law方差不增。完整有限来源边界仍由非该类原样估计；不存在丢弃边界后仍报全源量的问题。新的 pooled I₂ 统计结构已改变，不能沿用旧 outer CI；浮点法则/rounding边界也没有被 pooling 消除。

## 1. 原完整来源及物理 class

固定一轮 n，A=16n，k_ℓ=1,2,4,8，q_ℓ=2Ak_ℓ+1，w=(1/2,1/6,1/6,1/6)。原完整质量 W=1 的 measure 是

\[
 \mu=\sum_{\ell=0}^3w_\ell q_\ell^{-n}
 \sum_{m\in\{-Ak_\ell,\ldots,Ak_\ell\}^n}\delta_{m/k_\ell}.
\]

重合物理位置的四层贡献真实相加。定义

\[
 \mathcal C=\{y\in\mathbb Z^n:\ |y_i|\le A-2\ \text{对全部 }i\}.
 \tag{1}
\]

这取决于 physical y，不是 source_component=0：粗层边缘点不属于 C，细层抽到整数内部位置则属于 C。同一整数点抽到多个 labels，都有同一完整响应与 proposal，仍只是完整 μ 的正常抽样；不以 label 新领来源质量。

若保存的 source_k=k、integer_codes=c_i，正确成员判别为全部 c_i可被 k整除，且 |c_i|≤(A−2)k。`source_grid_membership[0]` 可核整数条件，但仍须另外核 margin。解析 class 质量为

\[
 p_{\mathcal C}=\mu(\mathcal C)
       =(2A-3)^n\sum_{\ell=0}^3\frac{w_\ell}{q_\ell^n}. \tag{2}
\]

每个 C物理点的质量相同；整数点数是[2(A−2)+1]^n，且每个点存在于全部四层。A=16n 时粗层主项趋于 (1/2)e^{-1/8}，细层项趋零，所以 p_C约0.44 的来源是严格质量公式，非估计或“重复源点”假设。本审计不重新计算三个样本的 class counts。

## 2. receiver/query halo 与严格整数平移

原 score 只在 source被 selected Q_R捕获时非零。R∈[1,2] 保证这个 receiver的 offset u=x−y满足 ||u||∞≤1，恰由原 core/cone prior覆盖。任何 candidate R'∈[1,2] 的 query都满足

\[
 z\in Q(y+u,R')\ \Longrightarrow\ |z_i-y_i|
          \le|u_i|+R'/2\le2. \tag{3}
\]

因此 y∈C时，所有 receiver及其 **全部 candidate queries** 都在原[-A,A]^n来源盒内。margin等于2的边界情况也合法：query闭面可以恰好到±A，该面原来源被包括；query中不需要任何盒外 grid点。没有只检查 selected query却忽略 runnerup/fullfuture尺度。

对每层，m/k_ℓ→m/k_ℓ−y 是整数 shift的格点双射，因为 k_ℓ y_i是整数。在 (3) 的局部 halo内没有 finite-source clipping，故对所有 u∈[-1,1]^n、全部 R'∈[1,2]，

\[
 h_{R'}*\mu(y+u)=h_{R'}*\mu(u). \tag{4}
\]

归一化始终是原 q_ℓ^(-n)，不是平移后另造的小盒概率。完整响应、全部 arrival边长、相等候选的最小 R tie规范、登记 τ、strict E=M>τ、以及原 source capture，逐 offset全部相同。因此

\[
 S(y)=\tau\int_E\frac{h_{R(x)}(x-y)}{M(x)}dx
      =S(0)=:s_0\qquad(y\in\mathcal C). \tag{5}
\]

式 (5) 没有假设全局 finite-input E 是平移周期集。积分中只有 y的捕获 halo有贡献；在这个 halo上用 (4) 就足够。盒外的其它 receiver、非C source、来源边缘及全部原质量，仍由原定义保留。原物理尺度窗口没有扩为[0,∞]，原 threshold也没有改为 infinite-grid threshold。

## 3. proposal 与 scores 的共同 conditional law

原 core/cone prior在 offset坐标本就不依 y。倾斜使用的 finite k=1 reference是 R*=3/2、high condition=捕获两个粗格点。free offset区间半宽最多1，加 reference query半宽3/4，halo至多7/4<2，故同样无 finite clipping。

整数平移后 `high_sets` 的 exact intervals减去 y即与 y=0完全相同；p_i、global product-mixture normalizer、free high statuses、固定 signed face及其诊断 status、RN 都有相同法则。T只随本轮 n改变，不随 class成员改变。给定相同 offset，oracle用了exact source+dyadic offset构造，故 source_component/source_k标签不会产生另一物理点或另一响应。于是 C内不仅 conditional means共同为s₀，理想 scores 的整个 conditional distribution共同。

run.py 的两个 rrng在每轮source循环 **之前** 初始化一次，各自在后续所有 receiver calls中向前消耗 fresh draws；未逐 source重设同一seed，也未把一个 saved offset block复制给多个 source。两个 streams与source stream分开。所用 independence仍是已登记 ideal experiment的条件概率假设，不由“seed不同”或输出不同证明实际伪随机数的数学独立性。

在此理想法则下，条件于全部 Y₁,…,Y_Ns，所有 receiver calls独立；每个 C source有两replicas各8 draws，全部16N_C个 scores彼此独立且同law。各调用在core/shell与base/tilted branch下用随机数的个数可以不同；这没有重用已消费的随机数，在ideal iid随机流中仍给 fresh suffix。若后处理发现重复整块 offsets是人为复制或重seed，则不能使用这个 conditional independence合同；本审计没有执行保存records的复查。

若source sampling本来重复抽到同一物理y，仍须保留该次正常样本及其fresh receiver draws。iid来源采样允许 coincident y；不应以去重改成另一 source law。

## 4. 随机 N_C 下的精确无偏校准

令 n_s=N_s为固定原source样本数，c=N_C为物理class成员个数。条件于全部 Y，若c≥1则N=16c≥16，不会出现N=1；记pooled scores Z₁,…,Z_N。原 Z∈[0,4C_n]，其 conditional mean都是s₀，C_n=1+nlog2。定义

\[
 U_{\mathcal C}=\frac{(\sum_{r=1}^NZ_r)^2-\sum_{r=1}^NZ_r^2}
                        {N(N-1)}
       =\frac1{N(N-1)}\sum_{r\ne t}Z_rZ_t. \tag{6}
\]

独立性给 E[Z_rZ_t|Y]=s₀²（r≠t），所以 E[U_C|Y]=s₀²。这是同一物理profile的交叉估计；未把任意不同profile来源相乘。实现 (6) 的大和减平方若产生舍入消去，可等价用正的pair累积2Σ_r Z_rΣ_{t<r}Z_t；不能因逐样本矩约束而裁剪估计、再宣称原无偏性。

若c=0，class项直接定义为0，不计算 (6)。对于非C source s，保留原两独立replica均值 A_s,B_s。新的全源 estimator为

\[
 \widehat J_2^{\rm pool}
       =\frac c{n_s}U_{\mathcal C}
         +\frac1{n_s}\sum_{s:Y_s\notin\mathcal C}A_sB_s. \tag{7}
\]

于是，包含c=0的全部source configurations，

\[
 \mathbb E[\widehat J_2^{\rm pool}\mid Y_1,\ldots,Y_{n_s}]
          =\frac1{n_s}\sum_{s=1}^{n_s}S(Y_s)^2.
\]

原source iid law为完整μ/W，故

\[
 \mathbb E\widehat J_2^{\rm pool}
          =\frac1W\int S(y)^2d\mu(y)=I_2/W. \tag{8}
\]

随机c不是bias来源：class项无条件期望为E[c/n_s]s₀²=p_Cs₀²。不能把c/n_s擅自换为解析p_C、仍在c=0把项设0；那会得到p_CPr(c>0)s₀²，少了p_CPr(c=0)s₀²。极小概率不等于数学上可删。原I₁ estimator不改，仍为所有source两replica平均的原均值；新的I₂/I₁比值仍一般有偏。

## 5. 理想方差确实不增，而不只是保持无偏

对固定c≥1，原C block的平均是

\[
 H_c=\frac1c\sum_{s\in\mathcal C}A_sB_s
       =\frac1{64c}\sum_{\text{同source跨replica的64 unordered pairs}}Z_rZ_t.
\]

把全部N iid scores随机置换，再分成c组、每组两replicas各8个。给定scores的multiset，每个unordered pair被纳入此H_c的概率相同；总有64c对。因此 H_c在这个随机分组下的conditional expectation恰为U_C。iid score law保证原分组在给定multiset下可交换，故U_C是H_c的 Rao–Blackwell化，conditional variance不增。

也可直接自证。令单score均值s₀、方差v。二阶 complete U-statistic给

\[
 \operatorname{Var}(U_C\mid Y)
        =\frac{4s_0^2v}{16c}
              +\frac{2v^2}{16c(16c-1)},
\]

而独立8-draw replica均值的产品给

\[
 \operatorname{Var}(H_c\mid Y)
       =\frac{s_0^2v}{4c}+\frac{v^2}{64c}.
\]

第一项相同，第二项的pooled/原比为8/(16c−1)≤1。推导来自相同index、只共享一个index、全部不同index的pair展开，无需新的外部U-statistic定理。

非C block条件于Y与C scores独立且未改；c/n_s条件于Y固定，所以 (7) 的conditional variance≤原全源I₂ estimator的conditional variance。二者conditional mean完全相同，total variance公式再给无条件variance不增。c=0时两个estimator相同。这个结论限定于已证明的C同law及fresh independent scores；只有common means却没有交换性时，(8)依然可以成立，但上述 Rao–Blackwell/方差公式不能原样宣称。

## 6. Post-hoc、固定预算、CI及最终边界

本class是input几何与source位置的固定规则，不按receiver scores选择成员。尽管在主run后才发现并应用结构，(7)作为对全部保存source/draws定义的固定统计函数仍满足 (8)。应按同一规则报告所有轮次，不能因改善原矩残差才选择某轮或某个score子集。

证明使用已完成的固定N_s、每source固定16 draws；若中断、缺失draws、按score截断或重复stream，则须重审条件采样和权重，不能把被保留样本继续视为本fixed-N模型。解析成员规则、完整源边界和coincident labels都不能在后处理中更换。

原registered结果保持原样；pooled值只能标为 **post-hoc结构降方差**。C项跨多个source共享一个U-statistic，已经不是原outer每source独立产品均值的形式，不能把旧I₂ Hoeffding区间半宽直接套在新值上。这里未推新的CI，也未认证实际浮点proposal law。原I₁估计没有改，但不能据旧联合CI宣传新的pooled矩对已有相同联合覆盖率。

所有无偏/方差结论对应原登记的ideal continuous exact模型。源格点和realized dyadic oracle的严格整数平移不补偿浮点offset/RN/score与连续法则的未认证误差；pooling不把既有理想统计区间升级为浮点实现的certified CI。本稿也不复查保存的N_C、Z sums或新估计值，root的只读后处理须另保存formula、counts、原profiles哈希与新结果。

审计结论为解析接口通过，无新的MC执行。它改善这个冻结四grid输入族的二阶估计，不提供任意L¹ source的尾合同、actual FIRST/CPGP/LCA/history认证或一般平方根阶数证据。
