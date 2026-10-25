# 全连续参数rho的Young证书族障碍

本文件只排除明确的辅助证书族，不给真实G的上界或C5的下界。模型固定五维alpha1、16条未单位化body方向，beta1/2。令M=A/T，A=ΣS_i v_i v_iᵀ，S独立且LT=e^(−sqrt lambda)，R由T^(−1/2)倾斜；Ztilt=1/(8sqrt pi)，c1=1/(32pi³)。原点以外G=c1 ER[d^(3/2)B_x^(−2)]，B_x=xᵀadjM x。

对每个x∈[-1,1]^5，令W_x=Σ_(j<k)M_jk xj xk、u_x=(W_x+2)/12。五维body多仿射顶点给0≤u_x≤1。考虑

Z=a0(1−rho u_x)^(−3/2)，rho<1，a0≥0。

允许a0、rho随空间x选择，但不允许依赖尚未积分的随机形状M。rho可为负；本证明覆盖全部rho<1。设lambda_i(x)=1−rho[(v_i·x)²−|x|²+4]/24>0，L=Σlambda_iP_i。

Young给G≥c1[3a0²U−2a0³V]，U=ER sqrt(d)L^(−3)、V=ER B_x L^(−9/2)。最大化a0后的该族预算为B_rho(x)=c1 U³/V²。以下给该下函数的上限，不可读作G≤B_rho。

## 1. 正lambda的完整闭式与Cauchy上包络

令H=Σlambda_i S_i，D=Σlambda_i^(−1/2)v_i v_iᵀ，C=Σsqrt(lambda_i)。尺度完全抵消为

U=Ztilt^(−1) E[sqrt(det A)H^(−3)]，
V=Ztilt^(−1) E[b_A H^(−9/2)]。

多线性Cauchy–Binet及Laplace导数给

E[detA e^(−tH)]=detD e^(−C sqrt t)/(32 t^(5/2))，
E[b_A e^(−tH)]=xᵀadjD x e^(−C sqrt t)/(16 t²)。

正Mellin积分和Cauchy给

U²≤Ztilt^(−2) E[detA H^(−5)] E[H^(−1)]
   =8pi detD/C^7，
V=(128/35)xᵀadjD x/C^5。

所有lambda严格正，负矩与带指数的半正矩均有限；积分均正，无T/P独立假设。该证明与030固定两rho结果共用基本LT模型，但此处重新推导后用于不同空间点及整个连续参数域。

## 2. x0=(1,1,1,0,0)的整个参数域

16方向的u_i(x0)有4个5/12、12个1/12。设lambda_s=1−5rho/12、lambda_l=1−rho/12。预算对全部lambda同乘正常数不变，因此可归一lambda_l=1，并令

t=sqrt((12−rho)/(12−5rho))，0<t<sqrt(11/7)<9/7。

实际16方向加权矩阵D的三维head对角为4t+12、非对角为4t−4，另两个tail对角为4t+12、其它交项0。因此

detD=16384(3t+1)(t+3)²，
x0ᵀadjD x0=12288(t+3)²，C=4/t+12。

代入上一节，严格有

B_rho(x0)≤[1225 sqrt2/(32768 pi^(3/2))]·[(3t+1)sqrt t/(144(t+3))]。

函数(3t+1)/(t+3)递增，t<9/7给其值<17/15、sqrt t<8/7。利用pi>3与sqrt(2/3)<5/6，得到完全有理的统一上端

B_rho(x0)<2975/31850496<1/10000，全部rho<1。

这里真实G(x0)=∞，并不矛盾：我们控制的是Young所生成下函数的预算，不是Green。

## 3. 正面上正测度patch：不能借奇点零测绕过

取eps=1/1000及||x−x0||2≤eps、x∈cube。每个W_i(x)=((v_i·x)²−|x|²)/2的梯度范数≤4sqrt5<9。对rho≥0、rho<0分别利用u_i(x0)∈{1/12,5/12}，均有

|lambda_i(x)/lambda_i(x0)−1|≤|W_i(x)−W_i(x0)|<9eps=:eta。

于是L_x≥(1−eta)L_0，U_x≤(1−eta)^(−3)U_0；D_x≥(1+eta)^(−1/2)D_0，C_x≤(1+eta)^(1/2)C_0。

对SPD矩阵，xᵀadjD x=|x|² det(D restricted to x-perp)，故D_x≥c D_0推出该二次型≥c⁴ xᵀadjD_0 x。D_0的head对称性又使x0始终为其（及adjD_0的）特征向量，因而

xᵀadjD_0 x≥(1−eps)² x0ᵀadjD_0 x0。

合得V_x≥(1+eta)^(−9/2)(1−eps)² V_0。因此全rho统一地

B_rho(x)≤[(1+eta)/(1−eta)]^9 (1−eps)^(−4) B_rho(x0)
        <(6/5)(1/10000)=3/25000。

最后一处纯Fraction比较为(1009/991)^9(1000/999)^4<6/5。

正面上的明确patch x1=1、x2,x3∈[3999/4000,1]、x4,x5∈[0,1/4000] 满足距离≤1/2000<eps，有正四维Lebesgue measure。因此即使空间自适应选择rho和a0，该族也不能在a.e.完整面上认证g≥15198/10^8（当前alpha1 I预算所需门槛更高）。证明没有对G作S2线外连续延拓，只用有限代数预算及统一lambda比较。

不排除不同lambda模式、随机形状自适应参数、Z求和/条件机制、不同I预算或真实G的新界。尤其该族失败不能反推C5≥任何数。

## 实际复核

profile_barrier.py用SymPy从16个方向独立构D、det与adj，正Mellin积分与预算平方作准确有理恒等检查，随后Fraction核所有参数与patch余量。实际PID40228，约0.108秒。首个结构相等断言及一次误用系统Python的修复均失败、没有输出证书，已单列failure回执；最终改用数学venv及simplify的准确零差检查后通过，没有改变数学输入。源码、输出及来源hash在JSON；还需独立审查完整patch链。
