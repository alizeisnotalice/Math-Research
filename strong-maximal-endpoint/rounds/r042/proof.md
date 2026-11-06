# 第42轮证明：无限终端参考下的左中间带全动作比较

日期：2026-11-06。**范围限于 a=mu、N=∞、无额外密度源和floor。** 本文不把此结论直接提升为匹配有限N的全空间定理。

## 1. 定理与剪切方向

令 mu=1+exp(−mu)、p=(1+exp(mu))⁻¹、A₃=p mu²/8=mu(mu−1)/8。沿用 Z_c=Exp(1)+c Bern(p_c)、K_ch(r)=Eh(r−Z_c)。记

v_t=K_mu^t min(1,exp(−r))，A_t=1−v_t，f_t=f_(t,mu)，

J_t=K_mu J_(t−1)+(p/8)f_t(r)(r−t mu)²/t²，J₀=0，

I_t=v_t exp(alpha J_t/v_t)，F_t=1−I_t，G_t=(F_t)_+，W_t=min(1,I_t)=1−G_t。

**定理。** 固定 alpha>1、C₀<∞。存在 A₀、t₀，使 t>t₀、

A₀ sqrt(log t/t)≤lambda≤C₀t^(−1/3)，r=t Lambda_mu'(−lambda)

时，对全部 c≥0 有 K_cW_(t−1)(r)≤W_t(r)。这里 A₀可依赖alpha；常数不声明对alpha趋1一致。

若F_t≤0，W_t=1而W_(t−1)≤1，结论自动成立。若F_t>0，只须证明

K_cF_(t−1)(r)≥F_t(r)。                                           (1)

因为K_cG_(t−1)≥K_cF_(t−1)，而G_t=F_t。**剪切有利**，无需把G当作三次可微函数。

下文均在活动点F_t>0，令q=alpha J_t/A_t、s=1−q。由exp(w)≥1+w知0≤q<1。记nu_lambda=Lambda_mu''(−lambda)=mu+O(lambda)、d=q/(nu_lambda t lambda²)。

## 2. 负倾斜 Mills 及真实代价导数缓冲

精确地

A_t(r)=∫₀^∞(1−e^(−u))f_t(r−u)du，
h_A(lambda)=∫₀^∞(1−e^(−u))e^(−lambda u)du=1/[lambda(1+lambda)]。

用负倾斜中心密度g_t，f_t(r−u)/f_t(r)=e^(−lambda u)g_t(−u)/g_t(0)。第41轮中心Fourier证明原样适用：倾斜参数在mu附近的负小区间，增量仍有统一指数矩、正方差和Exp平滑因子；终端正权(1−e^(−u))e^(−lambda u)/h_A的第j矩为O(lambda^(−j))。因此A/f=h_A[1+O((t lambda²)⁻¹)]，前三阶倾斜/物理导数可由同一积分微分。A'=∫₀∞e^(−u)f(r−u)du、A''=f−A'还给

A'/f=1/(1+lambda)+O(1/t)，A''/(lambda f)=1+o(1)。          (2)

Poisson精确式（第36轮）为

J_t/f_t=(p/8){t(y−mu)²+(C_t/t)(H_t−1)}，y=r/t，
f_t C_t/t=g_(t+2,t)+mu²p(1−p)g_(t,t−2)(r−mu)。

g_(j,k)为Gamma(j,1)+mu Bin(k,p)密度。辅助Fourier乘子及其前三阶导数有统一界；由此C_t/t=nu_lambda+O(1/t)，相应参数导数也有界。于是

q=alpha A₃t lambda³[1+O(lambda+(t lambda²)⁻¹)]+O(lambda log t)。 (3)

所有o(1)以下均随t→∞、在定理范围一致。将区域分两段，不能在整个范围直接声称前态q变化都是o(q lambda)：

- q≤1/2时，s≥1/2。由(2)–(3)微分可得F'/f=s+o(1)、F''/(lambda f)=s+o(1)、|F'''|≤C lambda²f；调和项产生的误差O(lambda log t)→0。
- q≥1/2时，(3)强制lambda≈t^(−1/3)。此时log t/(t lambda²)→0，且对q=alpha(J/f)/(A/f)做商式微分得到

  q'/q=−3/(nu_lambda t lambda)[1+o(1)]，
  q''=O(q/(t²lambda²))，q'''=O(q/(t³lambda³))。            (4)

这些式中的prime是r导数。利用(A−alpha J)=(1−q)A，保留q的导数，得到尖锐插值

F'/f=s[1+o(1)]+3d[1+o(1)]，
F''/(lambda f)=s[1+o(1)]+6d[1+o(1)]，
|F'''|/(lambda²f)≤C(s+d)。                                 (5)

在q≤1/2段也可用同样形式理解(5)：d=O(lambda+log t/(t lambda))=o(1)，所有次项由固定s缓冲吸收；在q≥1/2段，d=O(lambda)，(5)的余项相对s+d趋零。故(5)统一足以推出F'>0、F''>0、F''/F'=O(lambda)、|F'''|/F'=O(lambda²)，以及

F''²/(lambda² f F')≤s+o(1)。                               (6)

这里(6)不是从粗界F''/F'≤2lambda直接推出。接触段精确用

(s+6d)²/(s+3d)=s+9d+9d²/(s+3d)≤s+12d，d=O(lambda)。       (7)

粗比值界单独使用会损失因子4，无法处理alpha接近1。

**指数化余项。** 令E=v(exp(alpha J/v)−1)−alpha J，则F=A−alpha J−E。全部左侧v≥c>0，且J≤C log t/sqrt(t)；后者由左Gaussian密度上界及辅助密度sup界得到。活动带f_t(r)≤C t^(−1/2)exp(−c t lambda²)，而J/f≤C(t lambda²+log t)。取A₀充分大，直接微分E的解析式给E及所需前三阶导数相对(5)的正项为任意所需负幂误差，故(5)–(6)确实对真实F成立。并有

|E_t(r)|+K_mu|E_(t−1)|(r)=o(f_t(r)lambda²)。                 (8)

## 3. 归一化左移尾与前态一致性

令x=r−mu。对z≥0，负倾斜Fourier的全局中心sup界给

f_(t−1)(r−z)/f_t(r)≤C e^(−lambda z)，
A_(t−1)(r−z)/f_t(r)≤C e^(−lambda z)/lambda。

Poisson公式和同样辅助密度界给

J_(t−1)(r−z)/f_t(r)≤C[t+z²/t+log t]e^(−lambda z)。          (9)

r−z≤0时三者都为0。全左侧v≥c、J/v有界，因此|F_(t−1)|≤A_(t−1)+C J_(t−1)。由核密度q_c(z)≤2e^(−z)，在z>K log t外的归一化绝对积分可压到任意指定负幂；选择K后再增大t₀。这是相对f_t的尾控制，未拿绝对概率小量除以极小f_t。

在z≤K log t内，前态负倾斜与lambda相差O(log t/t)。接触段的q前态变化为O(log t/(t lambda))=o(d)；非接触段用s≥1/2吸收，而不声称同样的o(d)。所以(5)在x附近整段|u−x|≤K log t成立，且f_(t−1)(x)/f_t(r)=1+o(1)。F'在此段的比值为1+o(1)，因为lambda log t→0。由此所有后续Taylor余项均可用当前F'(x)归一化。

## 4. 全动作局部化与锐增益系数

以下F表示F_(t−1)。置L(c)=K_cF(r)。对有限c>0，直接对显式核微分得

L'(c)=p_c[L(c)−F(r−c)]，
L''(c)=p_c{(2p_c−1)[L(c)−F(r−c)]+F'(r−c)}。              (10)

mu处EZ=mu、Var Z=mu；m(c)=EZ_c在mu有唯一严格最大值，m''(mu)=−p。由(5)、(9)分割后Taylor，

L'(mu)=p mu F''(x)/2+O(lambda² F'(x))。

在固定小邻域|c−mu|≤delta，

L''(c)=p_c F'(x)[1+O(delta)+o(1)]>0。                       (11)

因为L(c)−F(r−c)=(c−m(c))F'(x)+O(lambda F'(x))，而c−m(c)=O(delta)。因此该邻域唯一极小点c_*满足c_*−mu=O(lambda)，并有

L(mu)−L(c_*)=p mu²F''(x)²/[8F'(x)] [1+o(1)]
             ≤A₃ f_t(r)lambda²[(1−q)+o(1)]。               (12)

(12)精确保留了(7)，没有粗化为四倍系数。

还需排除全部远动作。由于所有c≥0的Z_c前三阶矩统一有界，(9)尾分割及同一Taylor给

L(c)−L(mu)=[mu−m(c)]F'(x)
 +{E(mu−Z_c)²−mu}F''(x)/2+O(lambda²F'(x))。               (13)

误差对所有c一致，未对无限大的c本身作Taylor。邻域外mu−m(c)有固定正下界：c→0或∞时m(c)→1<mu。于是(13)的首项支配其余O(lambda F')项。全局最小动作确实在上述邻域，(12)覆盖所有c。

## 5. 精确源结算

因K_mu A_(t−1)=A_t及Poisson递推，

K_mu F_(t−1)−F_t
 =alpha(p/8)f_t(r)(y−mu)²+E_t−K_mu E_(t−1)
 =alpha A₃ f_t(r)lambda²[1+o(1)]。                         (14)

结合(12)，全部c满足

K_cF_(t−1)−F_t≥A₃ f_t lambda²[alpha−1+q−o(1)]>0。

先固定alpha>1，再增大t₀，即得(1)与定理。

## 6. 接触常数与未覆盖内容

在q趋1的接触尺度t lambda³→(alpha A₃)⁻¹，(5)给

F'/f=[3alpha(mu−1)/8]lambda+o(lambda)，
F''/(lambda f)=[3alpha(mu−1)/4]lambda+o(lambda)。

它们是实际指数化代价的缓冲，不只是线性替代式。上述证明闭合了无限终端参考的左侧活动中间带；但以下事项仍须另证：匹配a_N偏移、实际有限N的affine/call项及floor的统一误差；该带外早期左侧强制选择粗屏障；加入正式有界密度源和固定短时初始化后的全空间拼接。因此本轮不宣称一般C_n锐对数系数或乘法渐近等价已经得到。
