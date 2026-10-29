# 高维计数试验的来源、误差与尾包络

本文件说明实数算术下的组合链。程序的普通浮点四则、三角函数、BLAS、erfc及随机采样实现尚未向外认证；经验guard不是证书。严格条件矩引理和已核验BE定理，与这次有限数值pilot的证据资格分开。

## 1. 具体输入与真实质量

K255、B4、p=1/2、v8、M3、δ=1/(4Kv)。旧有限域标签y_j∈{0,1,2,3}由代码中全部域状态核验。plus格为j<127，a=127/255。第b表在plus使用(y_j+b) mod4，在minus使用(y_j−b) mod4，位置为 \((j+(y+1/2)/4)/K\)，箱峰全部留格。

令 \(\mu_b\) 为K个箱峰的等质量平均，采用同一个原低谱
\[
 D_M(x)=\sum_{m=1}^M\pi m c_m\cos(2\pi mx),\qquad
 \rho_b=1-p+p\mu_b+(p/K)D_M .
\]
系数在source_tables中完整指定。由实际公式有
\(c_m^2\le(23328/87)/(\pi^4m^4)\)，而
\(\sqrt{23328/87}<20\)、\(\pi>3\)，所以
\(\pi m c_m\le20/(3m)\)。故所有表积分1，且
\(\rho_b\ge1-p-(p/K)(20/3)H_M>0\)；这是来源正性，不是Gram平坦结论。

各表在第j格的质量相同：
\[
 q_j=\frac1K+\frac{p}{2K}\sum_m c_m
 [\sin(2\pi m(j+1)/K)-\sin(2\pi mj/K)].
\]
所以每表的真实plus概率均为
\[
 \eta=\frac{127}{255}+\frac{p}{2K}\sum_m c_m\sin(2\pi m\,127/255).
\]
系数三角帽足以保证η∈[.49,.51]，程序值约.4980478009。不能用1/2代替η。

固定整数切点为注册中的3338016、3339311、3340606。R4用四个连续整数段，R2只用中间切点，分别选表0/2。若S(y)为源坐标落在plus的个数，定义
\[
 F_R(y)=\sum_k\mathbf1_{\{S(y)\in I_{R,k}\}}
            \prod_{i=1}^n\rho_{\nu_R(k)}(y_i).
\]
每个产品律下S~Bin(n,η)，因此 \(\int_{\mathbb T^n}F_R=\sum_k\pi_{R,k}=1\)。
真实先验混合为 \(\bar f_R=\sum_k\pi_{R,k}\rho_{\nu_R(k)}^{\otimes n}\)，其中π是这个真实Binomial段概率。切点只是预先固定的整数，不声称为精确分位点；真实π由确定误差包携带。没有向gated来源再乘π，也不需要未知归一化Z。

## 2. query bins与条件独立

query每坐标原律为均匀圆周。先给定有序bin标签，各隐藏位置独立均匀。F_R、各产品及先验混合都对坐标置换不变，所以给定histogram counts后，可以任取一个规范有序bin分配计算条件矩；不声称坐标在histogram条件下本身独立。

类核 \(z_{q,b,\pm}(x)\) 是同一长度L窗口中 \(\rho_b\mathbf1_{C_\pm}\) 的质量除L。给定bins，取其真实均值m。正multilinear恒等式精确给
\[
 E(A_LF_R\mid N)=\sum_k A_{\nu_R(k)}
 P_{\mathbf p_{\nu_R(k)}}\{S\in I_{R,k}\},
\]
其中 \(A_b=\prod_q(m_{q,b,+}+m_{q,b,-})^{N_q}\)，
\(p_{q,b}=m_{q,b,+}/(m_{q,b,+}+m_{q,b,-})\)出现N_q次。
Bernoulli只表示多项式系数，不是把原非乘积来源误设为独立。条件方差的逐monomial和分层Minkowski证明见conditional/PROOF.md。

## 3. 核均值误差到PB参数

假设某表/尺度的两个真实类均值分别落在 \(\tilde m_\pm\pm e\)，且下端正。置s=Σ\(\tilde m\)、s_min为全bin下确界。比例对分子递增、对另一类递减，因此
\[
 \left|p-\frac{\tilde m_+}{s}\right|\le e/s.
\]
极端恰为 \((\tilde m_+\pm e,\tilde m_-\mp e)\)，其分母仍为s，不需另加虚构因子。故n个坐标的PB均值与方差误差都≤ne/s_min，后者用p(1−p)的1-Lipschitz性。总核产品落在 \(\tilde A_be^{\pm t_b}\)，其中
\[
 t_b=n[-\log(1-2e/s_{\min})].
\]
程序以这两个参数区间的角点包住normal CDF的自变量。真实η的prior也保留参数误差，使用实际Bin(n,η)，未改成Bin(n,.5)。

## 4. BE、均值与分层方差

任意L≥1的类核概率满足
\(\eta/(2-\eta)\le p_i\le2\eta/(1+\eta)\)，取类均值后仍成立。因此本源可安全使用p_i∈[.3,.7]。独立非同分布BE常数0.5606给CDF误差≤1/√n：代入三阶绝对矩预算，只需准确核对
\[
 (2803/5000)^2(29/50)^2/(21/100)<1.
\]
对整数段的两个半格CDF端点分别加减该帽，空段/全支持精确处理。normal参数矩形角点给其CDF上下包；两个有限端点最多付2/√n。这里使用的是已读原始BE定理，不是未支付误差的Gaussian proxy。

每表gated项真实均值于是有区间 \(l_b\le W_b\le u_b\)。若逐坐标两类核相对方差帽最大为t_(q,b)，则
\[
 r_b=\prod_q(1+t_{q,b})^{N_q}-1,\qquad
 \operatorname{Var}(A_LF_R\mid N)
 \le\left(\sum_k u_{\nu_R(k)}\sqrt{r_{\nu_R(k)}}\right)^2=:\widehat V .
\]
没有R²或monomial数量因子。固定产品取W_b=A_b；prior使用π的真实区间乘A_b，再作同一分层组合。区间的下/上权重分别正相加，不把π近似重新归一化而破坏包含关系。

## 5. 条件尾与尺度

记总均值区间[l,u]。l>λ时，Cantelli给下包
\((l-\lambda)^2/[\widehat V+(l-\lambda)^2]\)，否则取0。u<λ时，用
\(\widehat V/[\widehat V+(\lambda-u)^2]\)作上包，否则取1；再与Markov上包min(1,u/λ)取小。

有限尺度真max的下包为各尺度下包最大值；上包为各尺度上包之和截到1。下选择由counts决定，同一个L在形成全部表/类核后才选择，不假设尺度独立。这里J233为旧宏观自然网格，不是连续窗口的上界。

实际模型可另外允许L=1：单位质量周期源的A₁恒为1。程序的weak下界锚1来自λ向1的左极限；λ≥1时，该尺度的严格超水平事件为零，所以已计算的尾仍是原233点的尾。锚不能冒充裸I-window的结论。

## 6. counts提议、贡献帽与统计

以名义总核s_(q,b,j)>0定义Z_(b,j)=Q⁻¹Σ_qs。counts提议相对原Multinomial(n,1/Q)的准确理想RN为
\[
 D(N)=\tfrac12+\frac1{2RJ}\sum_{b,j}
       \prod_q[s_{q,b,j}/Z_{b,j}]^{N_q},\qquad W=1/D\le2.
\]
实现抽uniform分支或先选一个b,j再抽其product tilt counts。没有抽输入F_R本身，也无需把gated源当product抽样。

若C_Z≥maxZⁿ、C_E≥maxe^{t_b}，各模型每尺度的均值上包≤C_EΣ_b\(\tilde A_{b,j}\)，而
\(D\ge\sum_{b,j}\tilde A_{b,j}/(2RJ C_Z)\)。
下包为正的active尺度有λ<真实均值；上包由union/Markov控制。因此两方向的贡献均满足
\[
 0\le\lambda WB(N),\lambda WU(N)
       \le\min(2\lambda,2RJ C_ZC_E).
\]
不把仅对下方向成立的active条件套到上方向。

128个独立counts样本使用冻结阈值；8模型、21阈值、双方向共336项单侧Bernoulli-KL界，总α=.05。Chernoff预算来自[0,c]变量的凸性mgf支配，非观测最大值。区间包络是确定误差，不另付统计α。以上coverage以理想实数核/提议为前提；程序的完整FP未认证，故不能据本次数值CI宣称严格的新C_n下界。

## 7. 真正L¹输入的截断

取整数T≥2，令 \(F_{R,T}=F_R\mathbf1_{[-T,T]^n}/(2T)^n\)。周期质量1保证其非负L¹范数准确为1。中心在[-T+1,T−1]^n时，全部L≤11/9及L=1的cube包含于截断箱，平均值等于周期平均除以(2T)ⁿ。内域含(2T−2)ⁿ个完整周期；原周期尾系数乘 \((1-1/T)^n\) 后仍是Rⁿ真输入的下界。对固定n取T→∞恢复周期下界，无待证弱型常数的循环使用。这里只传递下界，不把周期有限网格上包变成全尺度上界。
