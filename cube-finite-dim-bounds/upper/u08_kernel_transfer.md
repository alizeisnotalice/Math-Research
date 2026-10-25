# U08：坐标Cauchy核比较，含无限Lévy源与全部尺度

2026-10-25，n3/4、Lebesgue，L=−Σ|D_j|。障碍分解使用u08_cauchy_obstacle.md，不套Delta H²/有限退出测度模型。

## 1. 准确目标与标准依赖

设K≥0、K(-x)=K(x)、K∈L1(R^n)，M=∫K>0。分布LK通过<K,Lphi>定义；要求LK≥0于R^n\{0}，不要求K光滑/有界，也不要求LK是有限测度。定义K_r(x)=r^-n K(x/r)。结论为

lambda |{sup_(r>0) K_r*|f|>lambda}|≤M||f||1，f∈L1。

若另有K≥2^-n1_[-1,1]^n，便给中心cube弱型上界M。没有这样的低成本核时，接口证明不会改善8/16。

唯一新增标准表示依赖是：连续实偶条件负定rho、rho(0)=0，有rho=xi^TQxi+∫(1−cos(xi·z))nu(dz)，Q≥0、nu对称且∫min(1,|z|²)nu<∞。所用表述已在[Gao作者论文《A Probabilistic Characterization of Negative Definite Functions》§5.1式(5.2)/(5.3)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7604715/)原文核到；不是声称重证该基础定理。下面正最小原则→条件负定直接验证，不需Courrège原文应用。

## 2. 原点cutoff付账，包含坐标轴奇性

Bphi=K*Lphi。若phi∈C_c∞、phi≥0、phi(0)=0，则grad phi(0)=0。取0≤eta≤1、eta=1于零点邻域，g_epsilon=eta(x/epsilon)phi，Taylor给

||g_epsilon||∞=O(epsilon²)，||D²g_epsilon||∞=O(1)。

真实坐标Cauchy有对称二差分表示

Lg(x)=(1/(2pi))Σ_j∫[g(x+s e_j)+g(x−s e_j)−2g(x)]/s² ds。

把每条轴的积分分为|s|≤epsilon、|s|>epsilon，分别用O(s²)、O(epsilon²)，得到||Lg_epsilon||∞≤C epsilon。因此<K,Lg_epsilon>→0，而phi−g_epsilon非负、支撑避0，假设给<LK,phi−g_epsilon>≥0。故Bphi(0)=<K,Lphi>≥0（偶性处理卷积反射）。这用的是每条坐标轴真实跳率，未漏掉轴上分布，也未只检查K的原点外点值。

Lphi∈C0且有界，K∈L1使Bphi∈C0；在固定测试支撑上由有限导数范数连续控制。平移性将零点结论移到任意x。

## 3. 直接条件负定符号

令k(xi)=∫e^(−i xi·x)K(x)dx（这是卷积乘子），rho=psi k，psi=Σ|xi_j|。k实偶连续，|k|≤M。任取有限xi_i、复系数z_i且Σz_i=0，置p(x)=Σz_i exp(i xi_i·x)，phi_R=chi(x/R)|p(x)|²，chi≥0、近0恒1。phi_R≥0、phi_R(0)=0，前节给<K,Lphi_R>≥0。

R≥1时phi_R及其二阶导全局一致有界，Cauchy二差分在|s|=1处分割给Lphi_R一致有界；同一小跳s²/大跳常数控制给逐点Lphi_R→L|p|²。K∈L1支配收敛，所以

−Σ_(i,j) z_i conj(z_j)rho(xi_i−xi_j)=<K,L|p|²>≥0。

因此rho连续条件负定、实偶且rho0=0；两点测试给rho≥0，从|k|≤M得0≤rho≤M psi≤M sqrt(n)|xi|。

调用标准Lévy–Khintchine表示，Q必须0：若某方向e有e^TQe>0，则rho(r e)≥r²e^TQe与上述线性增长矛盾。故rho=∫(1−cos(xi·z))nu(dz)，nu可能总质量无限。没有killing项（rho0=0），没有漂移（实偶对称）。不能把它写成nu−nu(R^n)delta0并把无穷大当常数。

## 4. 定义域与无限源的接触集极限

障碍定理给u∈Dom(L)。因rho≤M psi，rho u_hat∈L2，定义Bu为L2乘子−rho作用u。Fourier恒等式合法给

Bu=K*Lu于L2。

对delta>0，把nu截到|z|>delta为有限正测度nu_delta，令

B_delta u=nu_delta*u−nu_delta(R^n)u。

其符号rho_delta=∫_(|z|>delta)(1−cos(xi·z))nu(dz)非负、≤rho且逐点增加到rho。于是B_delta u→Bu于L2，靠rho u_hat∈L2支配收敛，而不是靠一致有穷跳率。

u≥0，u=0 a.e.于D^c，所以B_delta u=nu_delta*u≥0于D^c。有限测度卷积保正并由L2 Young界定义；即使nu集中在低维集合，Fubini(dx nu_delta)保证a.e.论证合法。取L2收敛子列，得Bu≥0 a.e.于D^c。无限源从未与u在零集外作未定义的“无穷大减无穷大”。

## 5. 同一坏集、所有共同尺度与一般L1

对有界紧支撑f≥0，取Cauchy障碍分解f=sigma−Lu，kappa=lambda/M，D={u>0}。于D^c，

K*f=K*sigma−Bu≤kappa M=lambda。

每个K_r也满足原点外L正性（order1齐次缩放），其乘子rho_r(xi)=psi(xi)k(r xi)=r^-1rho(r xi)，0≤rho_r≤M psi。因此同一个u属于全部B_r域，同一个sigma/D可用，不能给每个r重新选择坏集。

先对正有理r去掉可数并零集。对有界f，|K_r*f−K_s*f|≤||f||∞||K_r−K_s||1，尺度映射在L1连续，故在同一剩余集上扩到全部r>0。得到lambda|{sup_r K_r*f>lambda}|≤lambda|D|≤M||f||1。

一般f∈L1令f_m=min(|f|,m)1_[-m,m]^n↑|f|。非负核保证每个卷积及尺度sup单调增加，严格超水平集增加到目标集，输入质量也增加，常数M不变。若K≥2^-n1_cube，每尺度支配正是(2r)^-n cube积分，未改为矩形独立尺度。

## 6. 真实假设边界与恒定背景

u08分解使用真实全空间product-Cauchy半群的守恒、L∞收缩、s^-n transient上界和order1 Fourier控制；本比较还使用K偶、非负、L1、原点外分布正性。缺任一条件不能由有限矩阵正性替代。特别非局部零集Lu不是0，K*L的Levy源也不保证有限；这两处不能从已审Delta有限Green模型照搬。

Lebesgue进入对称Fourier/质量保守、Fubini(dx nu)、<u,Lchi_R>、sigma帽给K*sigma≤kappa∫K，以及kappa|D|≤∫sigma的Lebesgue接收体积。换w(x)dx后分母与卷积/生成元不再匹配，不推广到此前2^n坏背景。固定Lebesgue下原子输入平滑传递不被排除。

公开出处已按原文核：[Lin作者PDF Lemma3（PDF pp5–8，坐标stable二维）](https://coolrmal.github.io/assets/papers/c2-less-than-4.pdf)给相关比较路线；上述Gao论文§5.1核标准负定函数表示。本文针对n3/4 alpha1补上真实障碍极限及L2截断无限源链。没有新的低成本kernel实例或新U3/U4常数。
