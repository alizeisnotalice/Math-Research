# U09：体对角Cauchy finite-frame的完整比较迁移

2026-10-25。n=3或4，V={(1,s2,...,sn):sj=±1}，方向未作单位化，m=2^(n−1)。这是U08已审坐标Cauchy链的真实迁移，逐项付frame假设；不乘二维3.879，不把frame当cube与diamond的线性等价。Fourier乘子ψ(ξ)=Σ_v|v·ξ|，L=−ψ(D)。所有积分背景为Lebesgue。

## 定理

(a) 对f≥0、有界紧支撑、κ>0，存在u≥0，u∈L1∩L∞∩Dom_2(L)，及sigma∈L1∩L∞，满足

f=sigma−Lu，0≤sigma≤κ，sigma=κ a.e.于D={u>0}，∫sigma=∫f，κ|D|≤||f||1。

(b) 若K≥0、K(−x)=K(x)、K∈L1，且LK≥0于R^n\{0}分布，则M=∫K给

λ |{sup_(r>0)K_r*|f|>λ}|≤M||f||1，K_r(x)=r^-nK(x/r)，f∈L1。

所有r使用共同一个D；一般L1不加H1输入假设。若K≥2^-n1_cube，则cube中心共同半径极大算子常数≤M。

## 1. 真实半群、密度与能量域

取相互独立标准一维Cauchy C_t^v，令X_t=Σ_v v C_t^v。其特征函数为exp(−tψ)。因此P_t是对称、保正、保常数/质量、L2收缩的真正全空间卷积半群；跳率为每方向ds/(pi s²)，不是谱分数Dirichlet Laplacian。反向等价类方向仅计一次。

选n条独立方向作列矩阵B。对应n项和的密度为

p_t^B(x)=|detB|^-1 Π_j t/[pi(t²+(B^-1x)_j²)]。

剩余独立项为概率测度ν_t，真实p_t=p_t^B*ν_t，故p_t>0、连续、质量1，并且

||p_t||∞≤C t^-n，C=1/(4pi³)（n3），C=1/(16pi⁴)（n4）。

所用det为4及16，输入主frame证明已给具体basis。共同时间缩放使p_t(x)=t^-n p_1(x/t)。每个p_t∈C0（连续L1密度的上述卷积，或ψ~|ξ|的Fourier逆），且对正时间s，s↦p_s在L1连续：概率密度的缩放在L1连续。

Gram恒等Σ_v vv^T=mI给√m|ξ|₂≤ψ≤m|ξ|₂。故form域H^(1/2)，L2 operator域H1；没有用局部H²或C0嵌入。对h∈C_c∞，

Lh(x)=(1/(2pi))Σ_v∫_R [h(x+sv)+h(x−sv)−2h(x)]/s² ds。

二差分使Lh∈L1∩C0∩L∞，并给大cutoff||Lχ_R||∞=R^-1||Lχ||∞。小cutoff g_ε=η(x/ε)φ，若φ≥0、φ0=0，Taylor给||g_ε||∞=O(ε²)、||D²g_ε||∞=O(1)；按|s|=ε/|v|分割每方向积分，得||Lg_ε||∞≤Cφ εΣ_v|v|。所有方向线均已计入。

## 2. 分解存在、uniform L1与真实生成元域

记M_f=||f||1、F=||f||∞，f=0另取零解。0<t≤1，V_t=tΣ_(j≥0)P_(jt)f。由

||P_sf||∞≤min(F,C M_f s^-n)

及n>1，V_t有界、V_t−tf∈C0，P_tV_t=V_t−tf。递增迭代u_(t,0)=0、u_(t,k+1)=(P_tu_(t,k)+t(f−κ))_+被V_t控制，极限u_t为最小非负固定点且u_t≤(V_t−κt)_+。f紧支撑及V_t−tf∈C0使每个固定t的u_t有界且essential compact support，特别L1∩L2；支撑不要求对t一致。

定义L_t=(P_t−I)/t、sigma_t=f+L_tu_t。固定点给0≤sigma_t≤κ及sigma_t=κ于u_t>0；P_t保质量给∫sigma_t=M_f。准确配对为

E_t(u_t,u_t)+κ||u_t||1=∫f u_t，E_t=<u_t,(I−P_t)u_t>/t≥0。

右Riemann和上界为

||V_t||∞≤tF+[n/(n−1)](C M_f)^(1/n)F^(1−1/n)≤C_f。

所以||u_t||1≤M_f C_f/κ、||u_t||2²≤M_f C_f²/κ、||sigma_t||2²≤κM_f、||L_tu_t||2≤||f||2+√(κM_f)。V_t本身不要求L1，真正uniform L1来自能量。

ψ_t=(1−e^-tψ)/t。对Q>0及t≤1/Q，ψ≥Q处ψ_t≥(1−e^-1)Q，故u_hat_t的L2频率尾≤const Q^-2。低频H1有界、局部Rellich加该尾，给t_j↓0的strong L2_loc/a.e.子列；同时u_t、sigma_t弱L2。Fatou给u≥0、L1∩L∞。对Schwartz测试，L_tφ→Lφ于L2，得Lu=sigma−f；其Fourier右端L2给ψu_hat∈L2，即u∈Dom_2(L)=H1。

sigma仍在[0,κ]。在u>0的a.e.点u_t最终正，于是sigma_t最终κ，帽+DCT和弱极限给sigma=κ于D；没有用连续收敛。先以χ_R↑1测试弱极限得∫sigma≤M_f；再由u L1及||Lχ_R||∞=O(R^-1)，测试f=sigma−Lu回收∫sigma=M_f，排除质量逃逸。故κ|D|≤M_f。

V=∫_0∞P_sf ds为有界C0，短时Fδ/长时s^-n/中时L1核连续分块给V_t→V于L∞，故u≤V、essential-sup尾趋0；不能由此代替能量L1。方向截断|s|>ε总率2m/(pi ε)，符号ψ_ε≤ψ且趋ψ，故L_εu→Lu于L2。在u=0处L_εu为非负跳入积分，得Lu≥0 a.e.于D^c，一般不为0。

## 3. 给定K的无限源、同一D和一般L1

令Bφ=K*Lφ。第一节小cutoff配K L1，使原点外LK正性升级为φ≥0、φ0=0时<K,Lφ>≥0。对Σz_i=0的p(x)=Σz_i exp(iξ_i·x)，用χ(x/R)|p|²；二阶导及sup一致有界，逐方向二差分DCT再K L1给

−Σ_(i,j)z_i conjugate(z_j)rho(ξ_i−ξ_j)≥0，rho=ψ k，k(ξ)=∫e^-iξxK(x)dx。

因此rho连续实偶条件负定，rho0=0，且0≤rho≤Mψ。明确标准依赖是实偶Lévy–Khintchine表示（U08已核Gao§5.1）；线性增长排除Gaussian项，rho0排除killing，实偶排除漂移。得到rho=∫(1−cos(ξ·z))ν(dz)，ν允许总质量无限。

u∈Dom(L)且rho≤Mψ使u∈Dom(B)，Bu=K*Lu于L2。截ν到|z|>δ为有限ν_δ，B_δu=ν_δ*u−ν_δ(R^n)u，其符号0≤rho_δ≤rho趋rho，故B_δu→Bu于L2。在u=0处B_δu≥0；有限测度dx-Fubini合法，L2子列给Bu≥0 a.e.于D^c，未作无限常数相减。

每个K_r有rho_r(ξ)=r^-1rho(rξ)≤Mψ，因此同一u、sigma、D给K_r*f=K_r*sigma−B_ru≤κM于D^c。取κ=λ/M，先正有理r去掉可数零集；有界f下||K_r*f−K_s*f||∞≤||f||∞||K_r−K_s||1，L1缩放连续性扩至全部r。一般f用f_j=min(|f|,j)1_cube_j↑|f|：非负卷积和尺度sup交换单调极限，bad sets及质量单调收敛。此还给全尺度sup可测性，不假定一般L1输入本身有H1分解。

## 4. 明确范围

以上为有限spanning Cauchy frame的实际迁移，固定V_n已验证所需basis/Gram/密度界。恒定Lebesgue实际用于对称Fourier与保质量、dx-Fubini、cutoff自伴配对、K*sigma≤κ∫K及κ|D|≤∫sigma；不适用非恒定背景或独立矩形多尺度。本文件不提供低成本K。U09正部势核另见u09_positive_part.md；两者合并后仍需主g/J成本准确低于8/16才可能更新上界。
