# U08：真实坐标Cauchy离散时间障碍的连续极限

2026-10-25。仅n=3,4、Lebesgue。以下为分析证明，不把有限矩阵诊断当R^n极限。标准基础是Plancherel、L2弱紧性/局部Rellich及非负积分收敛；没有调用Delta局部H²/接触零集结论。独立审计状态在批报告更新。

## 1. 目标、符号与空间

令psi(xi)=Σ_j|xi_j|，L为L2 Fourier乘子−psi，

P_t=e^(tL)，p_t(x)=Π_j t/[pi(t²+x_j²)]。

P_t是对称、保正、保质量、L2收缩的真实全空间卷积半群。其坐标跳率为ds/(pi s²)，总率无限且集中在坐标轴；不是Dirichlet Laplacian谱分数幂。L2域是H1，因为|xi|≤psi≤sqrt(n)|xi|；能量E(u,u)=∫psi|u_hat|²的form域是H^(1/2)。L2 Fourier变换取unitary归一。

**定理。** f≥0有界紧支撑，kappa>0。存在u≥0，u∈L1∩L∞∩Dom(L)，以及sigma∈L1∩L∞，满足

f=sigma−Lu，0≤sigma≤kappa，sigma=kappa a.e.于D={u>0}，∫sigma=∫f。

故kappa|D|≤||f||1。不声称u紧支撑、sigma紧支撑或Lu在D^c为0。一般L1算子比较另经单调逼近，不先假定一般L1的该分解带H1正则。

## 2. 固定t的最小正障碍解真实存在

若f=0取u=sigma=0。其余记M=||f||1、F=||f||∞。0<t≤1，定义

L_t=(P_t−I)/t，V_t=tΣ_(j=0)^∞P_(jt)f。

由||P_s f||∞≤min(F,M/(pi^n s^n))及n>1，j≥1的级数对固定t一致收敛，V_t有界。V_t−tf是C0函数：每项P_(jt)f∈C0，且一致收敛。其精确半群恒等式为P_tV_t=V_t−tf，非负Tonelli给逐点/a.e.合法性，不要求V_t∈L1。

令u_(t,0)=0，u_(t,k+1)=(P_tu_(t,k)+t(f−kappa))_+。该映射保序，且T_tV_t=(V_t−kappa t)_+≤V_t，所以u_(t,k)单调增加、有上界V_t。非负积分收敛给有限最小固定点u_t，

u_t=(P_tu_t+t(f−kappa))_+，0≤u_t≤(V_t−kappa t)_+。

f紧支撑且V_t−tf∈C0，因此V_t在无穷远的essential sup趋0。于是固定t下u_t有essential compact support且有界，特别u_t∈L1∩L2。这里只对每个固定t成立，支撑半径允许随t恶化；绝不将它作为连续极限的统一紧支撑假设。

## 3. 帽、饱和和固定t质量均是准确恒等式

定义sigma_t=f+L_tu_t。u_t>0时固定点等式给sigma_t=kappa；u_t=0时

sigma_t=f+P_tu_t/t≥0，且固定点的零分支给sigma_t≤kappa。

所以0≤sigma_t≤kappa全空间成立，饱和只在D_t={u_t>0}。因为u_t∈L1且P_t保质量，sigma_t∈L1并且∫sigma_t=M。它在D_t外一般不等于f。

E_t(u,u)=<u,(I−P_t)u>/t≥0。与u_t配对的每项均有限，得到

E_t(u_t,u_t)+kappa||u_t||1=∫f u_t。

这条能量式是后续L1界的来源，不把V_t的点态有限误当V_t可积。

## 4. 统一L∞、L1、L2及生成元界，付清尾项

g(s)=min(F,M/(pi^n s^n))为递减函数，右Riemann和满足tΣ_(j≥1)g(jt)≤∫_0∞g(s)ds。因此

||V_t||∞≤tF+[n/((n−1)pi)]M^(1/n)F^(1−1/n)≤C_f，

C_f=F+[n/((n−1)pi)]M^(1/n)F^(1−1/n)。

故||u_t||∞≤C_f，能量式给||u_t||1≤MC_f/kappa，并且||u_t||2²≤MC_f²/kappa。帽和质量给||sigma_t||2²≤kappa M，所以

||L_tu_t||2≤||f||2+sqrt(kappa M)=C_2。

这些界对t∈(0,1]一致。没有将1/t总跳率视为一致有限：L_t本身的compound-Poisson跳率确为1/t，极限发散，只有作用于这些u_t的L2量被一致控制。

还可明确尾行为。V=∫_0∞P_s f ds是有界C0函数：短时间块的L∞范数≤F delta，长时间尾由s^-n可积控制，中间时间积分属于C0。对s远离0的L∞连续性用||p_s−p_r||1，而非假设P_tf→f在L∞成立，可证V_t→V于L∞。因此极限u≤V a.e.，有essential-sup尾趋0；V自身仍可能不在L1。真正u∈L1由上述能量界和下一节Fatou得到，不能由u≤V单独推出。

## 5. 强局部紧性和真实Dom(L)，没有局部H²捷径

psi_t=(1−exp(−t psi))/t。生成元界等价于||psi_t u_hat_t||2≤C_2。固定频率阈值Q>0，当t≤1/Q且psi≥Q时，psi_t≥(1−exp(−1))Q，因此

∫_(psi≥Q)|u_hat_t|²≤C_2²/[(1−exp(−1))²Q²]。

这是沿t↓0的渐近频率尾界，不声称每个固定t有统一H1范数。低频投影的H1范数有界，在任意有界球上用Rellich，结合该尾界得到u_t的strong L2_loc预紧性。取t_j↓0子列，同时u_t弱L2、sigma_t弱L2，并令u_t→u locally L2及a.e.。Fatou给u≥0、u∈L1∩L∞及上述L1界。

对任意Schwartz测试phi，psi_t≤psi且逐点趋psi，故L_tphi→Lphi于L2。于是从L_tu_t=sigma_t−f得到Lu=sigma−f的分布/弱L2恒等式。等价地在Fourier中psi u_hat=f_hat−sigma_hat∈L2，因此u∈Dom(L)=H1。也可用逆Fourier的光滑紧频率测试作为core；未对未经定义的Lu作形式代数。

弱L2极限的sigma仍满足0≤sigma≤kappa（对有界集合上的非负测试核）。对u>0的a.e.点，局部a.e.收敛使u_t最终正，故sigma_t=kappa最终成立；bounded sigma_t结合D上的紧支撑测试和DCT，给sigma=kappa a.e.于D。这个步骤不依赖H1在n3/4嵌入C0，也不把weak convergence自动当饱和。

## 6. 极限质量：先可积，再大尺度cutoff回收

对0≤chi_R≤1的光滑紧支撑径向cutoff，弱L2收敛给∫sigma chi_R≤M。取chi_R↑1，得sigma∈L1且∫sigma≤M，先付清可积性。

坐标Cauchy的真实缩放给||Lchi_R||∞=R^-1||Lchi||∞。u∈L1，所以

|<Lu,chi_R>|=|∫u Lchi_R|≤R^-1||u||1||Lchi||∞→0。

测试f=sigma−Lu并令R↑∞，得到∫sigma=M。这样排除sigma质量在weak极限中逃向无穷远；未假设u_t L1强收敛、统一support或uniform tail tightness。面积估计kappa|D|≤M随之成立。

最后f、sigma∈L2、u∈Dom(L)，可合法以u配对，得到连续能量恒等式E(u,u)+kappa||u||1=∫f u。

## 7. 非局部接触集的真实结论

令L_epsilon用每条坐标轴上|s|>epsilon的跳跃截断，总率2n/(pi epsilon)有限。其非负符号psi_epsilon≤psi且趋psi，所以u∈Dom(L)保证L_epsilon u→Lu于L2。

在D^c上u=0，L_epsilon u=(1/pi)Σ_j∫_(|s|>epsilon)u(x+s e_j)/s² ds≥0。非负积分随epsilon↓0增加；结合L2收敛子列，极限有限并等于Lu a.e.。因此Lu≥0于D^c，sigma=f+Lu≥f在那里，可能严格大于f。这是沿坐标轴的非局部尾/接触结算，不是Delta的Lu=0局部零集结论。

## 8. 范围与公开来源

证明实际使用n>1的transient s^-n时间积分，本文只登记n3/4。任意Markov半群未必满足该上界、频率尾控制或cutoff缩放；不能把本构造升格为任意Markov定理。有限矩阵容量反例另见计算批。初始输入限制f∈L∞且紧支撑；全L1比较在kernel传递文件最后处理。

公开定位：[Lin作者PDF《c2<4》](https://coolrmal.github.io/assets/papers/c2-less-than-4.pdf) Lemma2（PDF pp2–5）给二维坐标stable的变分构造，本文原文已读；其Lemma3（pp5–8）给比较机制。本文不是把原文的二维表述冒充已陈述n3/4，而是给上述独立离散时间/真实连续极限证明。[网页§3 Principle2](https://coolrmal.github.io/articles/two-partial-balayage-principles/)是所需接口的简述。此处未提供任何低成本核，理论闭合不产生新的U3/U4值。
