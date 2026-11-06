# 第43轮证明：匹配有限终端的早期左侧比较及带外选择

日期：2026-11-06。延续r042，以下保留实际有限N终端、匹配a_N、有界密度源及加性floor。全空间组合另见assembly.md。

## 1. 统一参数与真实仿射缺口

令mu=1+e^(−mu)、A₃=mu(mu−1)/8。匹配参数a=a_N趋于mu，满足b=m(a)−a≥0、b²≈log N/N；delta=mu−m(a)≥0、ell=2a−m(a)=m(a)−2b。取固定Lambda使kappa=Lambda b²−delta≥0。参考v_t=K_a^t phi_N、f_t=f_(t,a)，Poisson势J_t沿用r036。

固定alpha>1，取0≤s_t≤S<∞、Delta s_t≥0、Delta s_t=O(log^6t/t^(3/2))；正式r036源及常数源都满足。定义

U_t=alpha J_t+s_t f_t，V_t=v_t exp(U_t/v_t)，I_t=V_t+Lambda t b²/N，
L_t^a=1−r/N+t m(a)/N，L_t^mu=1−r/N+t mu/N，
A_t=L_t^a−v_t，E_t=V_t−v_t−U_t，
F_t=L_t^mu−I_t=A_t−U_t−E_t−kappa t/N。                       (1)

不能把A_t误写为1−v_t。对t≤sqrt(N)、r≤tm(a)，充分大N下r<rho_N；在全部可达前态r'≤r上有

0≤A_t(r')≤A_t^infty(r')=∫₀∞(1−e^(−u))f_t(r'−u)du。          (2)

因为实际终端缺口是psi_N(u)=1−u/N−e^(−u)，在0<u<rho_N非负，r'−u<0时密度为零；r'≤0时A_t=0。此早期区域B_t=0，故O_t^−=L_t^mu。

## 2. 早期有限N活动带定理

**定理1。** 对固定alpha>1、C₀<∞及上述源，存在A₀、t₀，使N≥t²、t>t₀、

A₀ sqrt(log t/t)≤lambda≤C₀t^(−1/3)，r=t Lambda_a'(−lambda)

时，全部c≥0满足

K_c min(L_(t−1)^mu,I_(t−1))(r)≤min(L_t^mu,I_t)(r)。          (3)

若F_t≤0，右侧为L_t^mu，(3)由K_cL_(t−1)^mu≤L_t^mu直接成立。若F_t>0，令q=U_t/A_t、s=1−q。活动性给

0≤q<1，kappa t/N≤A_t≤C f_t/lambda。                         (4)

只须证明K_cF_(t−1)≥F_t，因为正部只提高核平均，且

K_cL_(t−1)^mu=L_t^mu−[mu−m(c)]/N≤L_t^mu。

### 2.1 有限终端Mills及导数

负倾斜终端总质量精确为

h_(A,N)(lambda)=∫₀∞(1−u/N−e^(−u))e^(−lambda u)du
 =1/[lambda(1+lambda)]−1/(N lambda²)。                      (5)

相对无限终端缺陷为(1+lambda)/(N lambda)。N≥t²下统一有

b/lambda=O(t^(−1/2))，1/(N lambda)=o(lambda)，t/N=o(1)。      (6)

第一式由log N/N在N≥t²上的单调性得到。若N lambda≥2(1+lambda)，归一化后的有符号终端权总变差和第j绝对矩分别为O(1)、O(lambda^(−j))：分子直接由(1+u/N+e^(−u))e^(−lambda u)的Gamma积分界控制，分母至少h_(A,infty)/2。因此r041/r042的中心Fourier证明仍可逐次微分；不需要假装这个无限u权处处非负。其在实际支持u≤r<rho_N上是正的。

具体令R_t(x)=g_t(x)/g_t(0)，r=t Lambda_a'(−lambda)。有

A_t/f_t=h_(A,N)(lambda)[1+O((t lambda²)⁻¹)]，
|∂_lambda^j∂_a^k(A_t/(f_t h_(A,N))−1)|≤C/(t lambda^(2+j))，j+k≤3。

证明是对全轴|R_t(x)−1|≤C(|x|+x²)/t及归一化参数导数界，积分上述有符号权的绝对矩；时间导数界同r041，N固定。因(5)与无限终端相差O((N lambda)⁻¹)，其逆及前三阶导数也有相应相对误差。精确物理导数为

A_t'=R_exp−P(S_t<r)/N，A_t''=(1−1/N)f_t−R_exp，
R_exp=∫₀∞e^(−u)f_t(r−u)du。

由此A_t'/f_t=1/(1+lambda)+O(1/t+1/(N lambda))，A_t''/(lambda f_t)=1+o(1)。所有误差在(6)范围一致。

### 2.2 真正相位与源的缓冲

精确Poisson式为

J_t/f_t=(p_a/8){t(y−ell)²+(C_t/t)(H_t−1)}，
y−ell=−nu_a lambda+2b+O(lambda²)，nu_a=Lambda_a''(0)。

辅助密度乘子及其前三阶导数在a的固定紧邻域统一有界，故r042的微分论证加入O(b/lambda+(N lambda)⁻¹+t/N)扰动即可。源s_t f对q的贡献为s_t/(A_t/f_t)=O(S lambda)。在q≥1/2接触段，lambda≈t^(−1/3)，此项及调和项均为低阶；q≤1/2段由s≥1/2吸收。记d=q/[Lambda_a''(−lambda)t lambda²]，则在接触段

F_t'/f_t=s[1+o(1)]+3d[1+o(1)]，
F_t''/(lambda f_t)=s[1+o(1)]+6d[1+o(1)]，
|F_t'''|/(lambda²f_t)≤C(s+d)。                               (7)

非接触段前两式分别是s+o(1)，第三式有固定界。这里floor的空间导数为零。s_t在相邻时间的变化为o(lambda²)，故前态系数变化不破坏(7)；不能据此扩张到任意有界单调大跳源。

指数化余项保留：全部左侧v_t≥c>0，U_t≤C_S log t/sqrt(t)；当前f_t≤Ct^(−1/2)e^(−c t lambda²)。取A₀大，和r042同样的解析微分给E及前三阶导数的任意所需负幂归一化余项，且

|E_t|+K_a|E_(t−1)|=o(f_t lambda²)。                          (8)

于是(7)对真实F成立，而非仅对A−U。尖锐插值仍为

F''²/(lambda² f F')≤(s+6d)²/(s+3d)+o(1)≤s+12d+o(1)=s+o(1)。 (9)

### 2.3 全左移尾及全动作

r042的负倾斜全局密度比，在匹配a紧邻域统一成立。对z≥0有

f_(t−1)(r−z)/f_t(r)≤C e^(−lambda z)，
A_(t−1)(r−z)/f_t(r)≤C e^(−lambda z)/lambda，
U_(t−1)(r−z)/f_t(r)≤C_S[t+z²/t+log t]e^(−lambda z)。

r−z≤0时A,U均为0；F只剩常数−kappa(t−1)/N。活动性(4)控制该常数/f_t≤C/lambda。故物理核q_c≤2e^(−z)在z>K log t外的归一化绝对尾，对所有c一致压到任意负幂；没有遗漏有限floor的远左尾。

在z≤K log t内，接触段前态q越界量是o(d)，非接触段用固定s缓冲，和r042同一分段即可。相位导数因(7)满足F''/F'=O(lambda)、F'''/F'=O(lambda²)。将动作目标居中于mu，精确核导数及远动作Taylor完全同r042，得

K_muF_(t−1)−inf_cK_cF_(t−1)≤A₃ f_t lambda²[(1−q)+o(1)]。

将参考mu换回a只增加

O((a−mu)²F'+|a−mu|lambda F')=O(b²F'+b lambda F')=o(f_t lambda²)，

其中使用m'(mu)=0、Var'(mu)=p_mu mu及(6)。无穷远动作仍用Z_c的统一前三阶矩，而不对c本身作Taylor。于是

K_aF_(t−1)−inf_cK_cF_(t−1)≤A₃ f_t lambda²[(1−q)+o(1)]。       (10)

### 2.4 精确时间结算

K_a A_(t−1)=A_t，故

K_aF_(t−1)−F_t
 =alpha(p_a/8)f_t(y−ell)²+(s_t−s_(t−1))f_t
   +E_t−K_aE_(t−1)+kappa/N
 ≥alpha A₃ f_t lambda²[1+o(1)]+kappa/N。                    (11)

floor贡献是+kappa/N，不能因空间导数为零而删去。由(10)–(11)，

K_cF_(t−1)−F_t≥A₃f_t lambda²[alpha−1+q−o(1)]+kappa/N>0。

对每个c成立。正部剪切及mu−m(c)≥0即得(3)。证毕。

## 3. 早期左带外分支选择

**定理2。** 可取固定C_out、t₀，使t₀≤t≤sqrt(N)、r≤tm(a)−C_out t^(2/3)时I_t(r)≥O_t^−(r)。此处源只需非负。

取r₀=tm(a)−C_out t^(2/3)>0。r027已用有理数证书证明共同a的f_(t,a)对全部整数t≥4对数凹。阈值的负鞍点lambda₀满足lambda₀=(C_out/nu_a)t^(−1/3)[1+o(1)]。中心Fourier奇次抵消给g_t'(0)/g_t(0)=O(1/t)，对a紧邻域一致，所以

sigma₀=(log f_t)'(r₀)=lambda₀+O(1/t)≥c C_out t^(−1/3)。

对数凹性把score下界向整个左半轴传递。对0<r≤r₀，

A_t/f_t≤A_t^infty/f_t≤∫₀∞(1−e^(−u))e^(−sigma₀u)du
 =1/[sigma₀(1+sigma₀)]。

又|r−t ell|≥C_out t^(2/3)−2tb≥(C_out/2)t^(2/3)，且

J_t/f_t≥p_a(r−t ell)²/(8t)。

于是alpha J_t/A_t≥c alpha C_out³。先取C_out大，再取t₀大，可令alpha J_t≥A_t。early r<rho_N使B_t=0，而O_t^−−v_t=A_t+t delta/N。故

I_t≥v_t+alpha J_t+Lambda tb²/N≥O_t^−。

r≤0时A=J=0，floor直接支付t delta/N。于是整个左带外min(O_t,I_t)=O_t，且可直接使用全局粗屏障比较。此结论不要求Ominus小于Oplus。

本轮定理1上常数C₀可任意先固定。选择它严格大于C_out/inf nu_a，便与定理2有空间重叠；下端与r036自然管同样留固定重叠。全局组合不能省略这一常数选择。
