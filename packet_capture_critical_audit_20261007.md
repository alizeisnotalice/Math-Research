# Critical truewinner 的固定原 packet capture 审计

2026-10-07；指定6.1-sol high；只写新prefix。已读 [general_source_packet_square](general_source_packet_square_20261007.md)，这里只检验其固定packet资格/覆盖，不重新证明平方引理。registration在执行前保存。

**结果：固定四个全局tensor grid层packet在三种eta下全部fullness失败。** 24个原保存receiver、72个参数行、288个packet incidences均精确核算；所有q_dom为0。另有统一几何证明表明，这不是24点外推：对这些特定全局packets，窗口[1,2]内任意receiver都不可能满足eta>=1/100的fullness。

## 1. 原合同完全冻结

只读 critical_source_tail_arrival_20261007_results.json 与n4/n16/n64 exact gzip profiles，先核保存SHA-256及输入hash。原source参数、tau、complete continuous winner R(x)、full M、原E={M>tau}全部保留。没有import或执行旧oracle，没有新的receiver/source采样，没有重找候选或调整lambda。

原保存字段充分：有理receiver x、真实winner R、full M、n、A、k、q、正weights。预固定来源为原四个完整n维grid component measures

\[
 \mu_i=w_i\,q_i^{-n}\sum_{z\in\{-A,-A+1/k_i,\ldots,A\}^n}\delta_z,
 \quad(k_i)=(1,2,4,8),\quad A=16n,
 \quad q_i=32nk_i+1,
 \quad(w_i)=(1/2,1/6,1/6,1/6).
\]

所以M_i=w_i、sum mu_i=mu、W=1。coincident物理atom上的source fractions为a_i(z)=mu_i({z})/mu({z})，sum a_i=1；没有重新分配或归一化。packet标签固定为这四个全局来源层，没有按receiver细分、挑新支撑或拆同位标签制造合作。

postprocessing只用原闭Q_x计算

\[
 c_{ij}=\#\{m/k_i\in[x_j-R/2,x_j+R/2]:-Ak_i\le m\le Ak_i\},
 \quad m_i=w_iq_i^{-n}\prod_jc_{ij},
 \quad m=\sum_i m_i=M R^n.
\]

最后等式在24点全部以Fraction准确核对；所有同距离/closed边界由exact ceil/floor计数，不使用容差。它核capture数据一致性，不重新认证原winner。

## 2. 资格与全部失败量

三轮n4/16/64对应theta=1/2,1/4,1/8，theta²=1/n；eta=1/2,1/10,1/100。dominance使用m_i>=theta*m，fullness使用m_i>=eta*M_i，两条等号均进入paid branch，失败为严格小于。记录

\[
 d_i=m_i/m,\quad f_i=m_i/M_i,\qquad
 q_{dom}(x)=1_E(x)\sum_i d_i1_{\{d_i\ge\theta,\ f_i\ge\eta\}}.
\]

完整原m作分母，没有改成selected packet质量。另保存忽略E时的raw qualified capture fraction，以区分资格失败与原超阈gate。每行捕获份额准确分为四项：paid、仅dominance失败、仅fullness失败、两项皆失败，总和为1；再保存dominance失败总份额和fullness失败总份额、各packet pass/fail及个数。没有把两个重叠失败量相加后当成残余份额。

结果按每一种eta相同：

|n|原receiver点|原E内点|全部点dominance通过incidences|E内dominance通过incidences|fullness失败incidences|
|---|---:|---:|---:|---:|---:|
|4|8|4|7|4|32/32|
|16|8|1|8|1|32/32|
|64|8|1|11|1|32/32|

跨三eta共有72个row records、288个packet tests；E内共18个参数行、72个packet tests，fullness仍72/72失败。q_dom与raw qualified share准确均为0，fullness失败捕获份额准确为1，所有captured incidence留在residual。

原E内reference k1层的dominance份额：n4约.63047至.77237；n16约.9950612；n64约.999999995445。与此同时，全部保存点中最大的单packet fullness仅为

|n|最大观测m_i/M_i|
|---|---:|
|4|5.86823e-8|
|16|4.80669e-39|
|64|1.06220e-193|

这些小数只是已保存精确Fraction的显示，不参与任何判定。可见真实rarelevel receiver可能几乎全由一个全局层供给，却只截获该层极小的总体质量；dominance高与fullness高完全不同。表中点/incident计数不能解释为Lebesgue receiver覆盖概率，也不是I_dom/I1的样本估计。

## 3. 该固定allocation的uniform fullness障碍

任何R<=2的闭坐标interval在spacing1/k_i的grid中至多捕获2k_i+1个点，有限source边界只会减少个数。因此，对每个receiver，不依赖其winner或tau，

\[
 \frac{m_i(x)}{M_i}\le
 \left(\frac{2k_i+1}{32nk_i+1}\right)^n. \tag{1}
\]

对n>=4,k>=1，该比值随k减小达到最大，且(2k+1)/(32nk+1)<=3/129=1/43。故三轮的任何packet都满足

\[
 m_i/M_i\le43^{-4}=1/3418801<1/100.
\]

于是三种eta的E_i全为空；该allocation的S_dom也准确为0，其平方引理没有截出正交通。这个uniform结论来自实际full n维source与cube count几何，不是把24个有理点当Lebesgue覆盖。若原I1>0，该allocation不可能给general_source_packet_square所需的正kappa平均覆盖；并不需要假定这24点是随机样本。

它**不**反驳一般packet平方引理或其它合法预固定空间packet的覆盖可能。四个全局分辨率层横跨整个[-16n,16n]^n，窗口cube只访问其中极小比例，本来就不具备local fullness。后续若改用固定空间cells/反链/共享source权，必须另登记真实mu_i、M_i及不重计source，重新实算完整分母；不能把本批层标签任意细分后仍声称沿用旧资格，也不能逐receiver临时重组packet。

## 4. 完成与范围

run正常exit0，耗时约.0254秒；只有精确Fraction后处理，无新seed/CI，无oracle重跑。input与profile hashes、24个saved full response一致性、全部72份四类失败partition均复核，原文件不修改。

本次回答的是：**现有critical样本及其完整source的这四个原global-layer packets不给fullness/dominance joint覆盖。** 未计算Lebesgue积分、I_dom/I1、S_res平方或一般kappa，也未验证general lemma的证明、宣称一般sqrt(n)预算或actual geom反例。registration/script/results/receipt以同prefix保存。
