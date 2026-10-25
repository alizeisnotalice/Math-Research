# 可直接写入两篇正文的有限源比较证明

本说明是第二审新写的自足分析路径；不把024/026/031的审计标签作为证明。一般核的CND/无限源审阅另见audit.md。本路径用于最终正部核，故上界主定理可以不以Lévy–Khintchine为必要依赖。所有积分均为Lebesgue，L2 Fourier变换可取unitary；卷积核/有限测度的变换取\(\widehat K(\xi)=\int e^{-ix\cdot\xi}K(x)dx\)，故卷积乘子没有额外\((2\pi)\)因子。

## A. 真实模型和障碍分解

取有限非零向量\(v\)及正权\(c_v\)，其线性张成是\(\mathbb R^n\)，这里\(n=3,4\)。设
\[
 \psi(\xi)=\sum_v c_v|v\cdot\xi|,\qquad L=-\psi(D),\quad P_t=e^{tL}.
\]
独立一维标准Cauchy变量之和\(X_t=\sum_vvC^v_{c_vt}\)给真实正、自伴、保质量的卷积半群。选任一方向基B，前n项具有严格正密度；其余项为概率卷积。因此\(p_t>0\)、\(p_t\in C_0\)、\(\|p_t\|_\infty\le Ct^{-n}\)。另有\(p_t(x)=t^{-n}p_1(x/t)\)、正时间L1连续及\(a|\xi|\le\psi(\xi)\le b|\xi|\)。从而\(\mathrm{Dom}_2L=H^1\)、形式域\(H^{1/2}\)。本项目选body基，使\(C=1/(4\pi^3)\)或\(1/(16\pi^4)\)，增加三维权1/3的face方向只作额外概率卷积。

对\(\phi\in C_c^\infty\)，
\[
 L\phi(x)=\frac1{2\pi}\sum_vc_v\int_{\mathbb R}\frac{\phi(x+sv)+\phi(x-sv)-2\phi(x)}{s^2}\,ds.
\]
二差分由\(\min\{s^2|v|^2\|D^2\phi\|,4\|\phi\|\}\)控制，故\(L\phi\in L^1\cap C_0\cap L^\infty\)，且大cutoff满足\(\|L\chi_R\|_\infty=R^{-1}\|L\chi\|_\infty\)。

**障碍引理。** 对非负有界紧支撑f及\(\kappa>0\)，存在\(u\ge0\)、\(u\in L^1\cap L^\infty\cap H^1\)及\(\sigma\in L^1\cap L^\infty\)，使
\[
 f=\sigma-Lu,\quad 0\le\sigma\le\kappa,\quad \sigma=\kappa\text{ a.e. on }D=\{u>0\},\quad \int\sigma=\int f.
\]
特别\(\kappa|D|\le\|f\|_1\)。

证明：记\(F=\|f\|_\infty\)、\(M_f=\|f\|_1\)。对\(0<t\le1\)，令
\[
 V_t=t\sum_{j\ge0}P_{jt}f,\quad
 u_{t,0}=0,\quad u_{t,k+1}=(P_tu_{t,k}+t(f-\kappa))_+.
\]
\(\|P_sf\|_\infty\le\min(F,CM_fs^{-n})\)保证固定t的V_t有界，\(V_t-tf\in C_0\)，且\(P_tV_t=V_t-tf\)。迭代单调趋于固定点\(u_t\)，满足\(0\le u_t\le(V_t-\kappa t)_+\)。这先保证每个固定t的u_t有essential紧支撑、属于L1及L2；不要求V_t可积或u_t支撑一致。

令\(L_t=(P_t-I)/t\)、\(\sigma_t=f+L_tu_t\)。固定点两分支给\(0\le\sigma_t\le\kappa\)及\(\sigma_t=\kappa\)于\(u_t>0\)。由u_t可积及半群保质量得\(\int\sigma_t=M_f\)。与u_t配对得
\[
 \langle u_t,-L_tu_t\rangle+\kappa\|u_t\|_1=\int fu_t.
\]
递减函数右Riemann和给
\[
 \|V_t\|_\infty\le F+\frac n{n-1}(CM_f)^{1/n}F^{1-1/n}=:C_f.
\]
故\(\|u_t\|_1\le M_fC_f/\kappa\)、\(\|u_t\|_2^2\le M_fC_f^2/\kappa\)、\(\|\sigma_t\|_2^2\le\kappa M_f\)、\(\|L_tu_t\|_2\le\|f\|_2+\sqrt{\kappa M_f}=:C_2\)。

写\(\psi_t=(1-e^{-t\psi})/t\)。当\(t\le1/Q\)时，\(\psi\ge Q\)处\(\psi_t\ge(1-e^{-1})Q\)，因而\(\int_{\psi\ge Q}|\widehat u_t|^2\le C_2^2/((1-e^{-1})^2Q^2)\)。低频投影的H1界和局部Rellich，结合该尾估计，给子列\(u_t\to u\)于L2_loc及a.e.，同时\(u_t\rightharpoonup u\)、\(\sigma_t\rightharpoonup\sigma\)于L2。Fatou给u的L1/L∞界。\(L_t\phi\to L\phi\)于L2允许传递分解；Fourier得到\(\psi\widehat u\in L2\)，即u属于H1。在u>0的a.e.点最终u_t>0，故\(\sigma_t=\kappa\)最终成立；对D上的有界支撑测试用DCT与弱收敛得到饱和。

最后先用\(\chi_R\uparrow1\)测试弱极限，得\(\sigma\in L1\)、\(\int\sigma\le M_f\)。再由\(u\in L1\)及大cutoff估计，\(|\langle Lu,\chi_R\rangle|\le R^{-1}\|u\|_1\|L\chi\|_\infty\to0\)，从\(f=\sigma-Lu\)回收\(\int\sigma=M_f\)。这排除质量逃逸，不使用u_t的L1强收敛或H1嵌入连续函数。

## B. 正部势的全空间有限源

令\(G(x)=\int_0^\infty p_t(x)dt\)为扩展非负函数，\(c>0\)，\(H=(G-c)_+\)。写\(G=U_T+W_T\)，其中\(\|U_T\|_1=T\)、\(\|W_T\|_\infty\le CT^{1-n}/(n-1)\)。取\(T=(2C/((n-1)\lambda))^{1/(n-1)}\)，得\(|\{G>\lambda\}|\le A\lambda^{-n/(n-1)}\)。这给G局部L1及H全空间L1，包括方向线邻域及全部远尾。

对\(\phi\in C_c^\infty\)，
\[
 \int G|L\phi|\le\int_0^\infty\min(\|L\phi\|_\infty,Ct^{-n}\|L\phi\|_1)dt<\infty.
\]
半群导数和\(P_T\phi(0)\to0\)于是给\(LG=-\delta_0\)。由正Tonelli与半群律，\(P_tG=\int_t^\infty p_sds\)，故a.e.有\(G=P_tG+U_t\)。保正保常数及正部1-Lipschitz性质给
\[
 P_tH\ge(P_tG-c)_+\ge H-U_t.
\]
测试非负\(\phi\in C_c^\infty\)，用自伴性除以t，并令\(t\downarrow0\)。左边由H L1和\((P_t\phi-\phi)/t\to L\phi\)于L∞，右边由\(U_t/t\rightharpoonup\delta_0\)，得到\(LH+\delta_0\ge0\)于整个空间分布。因此\(\mu=LH+\delta_0\)是正Radon测度。此处没有删去方向线或level界面，也不假设H具有紧支撑。

用\(\chi_R\uparrow1\)及H L1，
\[
 \int\chi_Rd\mu=1+\int HL\chi_R\longrightarrow1,
\]
故\(\mu(\mathbb R^n)=1\)。令\(h=\min(G,c)\)，有\(\mu=-Lh\)且\(0\le h\le c\)。小cutoff的准确缩放\(\|L\chi_\epsilon\|_1=\epsilon^{n-1}\|L\chi\|_1\)给
\[
 0\le\int\chi_\epsilon d\mu\le c\epsilon^{n-1}\|L\chi\|_1\longrightarrow0,
\]
所以\(\mu(\{0\})=0\)。正源总质量与原点负原子均已结算。

## C. 有限源直接比较：同一个D用于所有r

**有限源比较引理。** 设\(K\ge0\)、\(K\in L1\)，且全空间\(LK=a(\mu-\delta_0)\)，其中\(a>0\)、\(\mu\)为正概率测度。令\(M_K=\int K>0\)。则
\[
 \lambda\big|\{x:\sup_{r>0}K_r*|f|(x)>\lambda\}\big|\le M_K\|f\|_1,
 \quad K_r(x)=r^{-n}K(x/r),\quad f\in L1.
\]

证明：\(LH\)及\(LK\)作为tempered分布的Fourier变换确为\(-\psi\widehat H\)、\(-\psi\widehat K\)：对Schwartz测试，由\(\psi|\widehat\phi|\in L1\)和K L1可直接Fubini验证。全空间有限测度等式因此给连续函数恒等式
\[
 -\psi(\xi)\widehat K(\xi)=a(\widehat\mu(\xi)-1).
\]
对任意\(u\in H1\)，两边乘\(\widehat u\)，用K L1及有限测度L2 Young界，得
\[
 K*Lu=a(\mu*u-u)\quad\text{in }L2.
\]
这一步避免不合法地把非紧支撑分布随意卷积。若\(\mu_r=(z\mapsto rz)_\#\mu\)，order1齐次给
\[
 LK_r=\frac a r(\mu_r-\delta_0),\qquad
 K_r*Lu=\frac a r(\mu_r*u-u)\quad\text{in }L2.
\]

先令f非负、有界紧支撑，取A节障碍分解，\(\kappa=\lambda/M_K\)。在同一个\(D^c=\{u=0\}\)上，对每个固定r，有限正测度卷积\(\mu_r*u\ge0\) a.e.（即使测度集中在方向线上也由\(dx\,d\mu_r\)-Fubini成立）。故
\[
 K_r*f=K_r*\sigma-K_r*Lu\le\kappa M_K=\lambda\quad\text{a.e. on }D^c.
\]
先合并正有理r的可数零集，再用
\[
 \|K_r*f-K_s*f\|_\infty\le\|f\|_\infty\|K_r-K_s\|_1\longrightarrow0
\]
扩至全部r。这证明严格坏集包含于D（差一个零集），因而\(\lambda|\mathrm{bad}|\le\lambda|D|\le M_K\|f\|_1\)。同一个D只依赖f及\(\lambda/M_K\)，从未随r重选。

一般f取\(f_j=\min(|f|,j)1_{[-j,j]^n}\uparrow|f|\)。非负卷积的MCT及\(\sup_r\sup_j=\sup_j\sup_r\)给极大函数单调增加，严格超水平集亦单调增加，输入质量趋\(\|f\|_1\)。这同时给全尺度sup的可测性。无需为一般L1输入先证明H1障碍分解。

本项目\(K=H/N\)、\(N=2^n(g_0-c)\)，故\(a=N^{-1}\)。原点无正原子的结论已证明，但有限源直接比较本身只用\(\mu\ge0\)、\(\mu\)有限及全空间源等式。

## D. 从连续g/I证书到常数

设\(Q=[-1,1]^n\)，有严格连续证书\(0<g_0\le\operatorname{essinf}_QG\)。取\(c=g_0/n\)、\(N=2^n(g_0-c)\)。则\(K\ge2^{-n}1_Q\) a.e.，所以所有中心轴平行方体共同半径平均均被\(K_r*|f|\)支配。

令\(q=n/(n-1)\)、\(J=\sum_{F\subset\partial Q}\int_FG(y)^q\,dS\)。G的齐次性与weak-Lq先保证\(J<\infty\)。逐face沿\(x=ry\)的Jacobian为\(r^{n-1}\)，得到
\[
 \int H=\frac{n-1}{n}c^{-1/(n-1)}J.
\]
故
\[
 C_n\le\int K=\frac{n^{1/(n-1)}J}{2^ng_0^q}.
\]
body、face各方向类均在cube群下保权，因此若\(I=\int_{[0,1]^{n-1}}G(1,y)^qdy\)，则\(J=n2^nI\)。认证\(I\le I_0\)给
\[
 C_n\le n^{n/(n-1)}I_0/g_0^{n/(n-1)}.
\]
因此三维平方接口为\(C_3^2\le27I_0^2/g_0^3\)，四维立方接口为\(C_4^3\le256I_0^3/g_0^4\)。g0/I0的具体连续认证不属于本文，不能以样本minimum或浮点积分代替。
