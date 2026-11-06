# 第39轮：固定短时初始化与中间时间尺度反例

日期：2026-11-06。状态：固定时刻全空间初始化已证明；现有二次非线性势在接近锐系数时的直接全空间延拓有严格反例。一般锐对数系数仍未解决。本轮使用本机 Wolfram 15，未启动 MATLAB。

## 1. 账本与本轮选择

统一 \(\Phi_n(t)=t(1+\log^+t)^{n-1}\)，\(\mu=1+e^{-\mu}=\log B\)，\(\sigma=\sqrt\mu\)，\(\bar\Phi_G(z_n)=\sigma z_n/\sqrt n\)，\(d_n=\sigma\sqrt n z_n\)。一般严格账本不变：
\[
 \log Q_n=n\mu-d_n+\frac{9-\mu}{24}\log n+o(\log n),\qquad
 Q_n\le C_n\le\mathcal R_n,\quad
 \log C_n=n\mu-d_n+O(\log n).
\]
第29轮的 \(g=(\mu-1)/8\) 是反馈上界模型相对 Q_n 的对数级收益下界，不能登记为真实 C_n 下界。

第37–38轮已将剩余问题定位为初始化及中间时段。故本轮先证明固定短时初值可由有限密度校正支配，再检查该形状能否直接延拓；同时对真实 source 权重与 receiver 比例作联合扫描。结论是：固定初值可以完成，但在 \(t\sim\sqrt{\log N}\) 必须修改尾部形状，不能仅令固定膨胀系数趋于1。

记 \(\mathcal U_{t,N}=T^t\phi_N\)，\(T=\sup_{c\ge0}K_c\)，\(K_ch(r)=Eh(r-Z_c)\)，\(Z_c=\operatorname{Exp}(1)+c\operatorname{Bern}(p_c)\)、\(p_c=(1+e^c)^{-1}\)；\(\phi_N(r)=e^{-r}\)（r>0）、\(1-r/N\)（r≤0）。匹配 a=a_N→mu，\(m(a)=EZ_a\)、\(b=m(a)-a\)、\(\delta=\mu-m(a)=O(b^2)\)、\(b^2\asymp\log N/N\)。参考 \(v_t^a=K_a^t\phi_N\)、密度 f_t^a，\(\ell=2a-m(a)\)。精确二次 Poisson 势为
\[
 \frac{J_t^a}{f_t^a}=\frac{p_a}8
 \left\{\frac{(r-t\ell)^2}{t}+\frac{C_t}{t}(H_t-1)\right\}\ge0.
 \tag{1}
\]
H_t 是调和数，C_t 是第36轮已给出的条件桥方差量。

## 2. 固定短时全空间初始化（严格）

**命题。** 给定固定整数 k≥4、alpha>0，存在有限常数 s0、Lambda，独立于 N，使所有充分大匹配 N 和所有实数 r 都满足
\[
 \mathcal U_{k,N}(r)\le
 I_{k,N}(r):=v_k^a(r)\exp\{[\alpha J_k^a(r)+s_0f_k^a(r)]/v_k^a(r)\}
 +\Lambda kb^2/N.
 \tag{2}
\]
因此 \(\min(O_k^-,O_k^+,I_{k,N})\) 也支配真值。这里 s0 只是一项有限常数密度校正，不是新增 log N 源。

**证明。** 把 a 限于 mu 的一个预先固定的正紧邻域。取 \(\Lambda b^2\ge\delta\)，它由 \(\delta/b^2\) 一致有界保证。

r≤0 时增量非负，精确地
\(\mathcal U_{k,N}=1-r/N+k\mu/N\)、\(v_k^a=1-r/N+km(a)/N\)，而 f=J=0；(2) 的仿射余量吸收差值。

取固定 \(0<r_0<\min(\inf a/2,1/2)\)。对充分大 N，所有 x≤r0 有 \(\phi_N(x)\le1-x/N\)，故 \(\mathcal U_{k,N}\le L_k^\mu=1-r/N+k\mu/N\)。当 0<r<r0 时，仅零 Bernoulli 位移项参与
\[
 L_k^\mu-v_k^a=k\delta/N+D(r),\quad
 D=E\{[1-e^{-(r-S_k^a)}-(r-S_k^a)/N]1_{S_k^a<r}\}.
\]
由 integrand≤r−S 及该 Gamma 密度，
\[
 0\le D/f_k^a\le e^{r_0}\frac{r^2}{k(k+1)}.
 \tag{3}
\]
而 \(I-v\ge\alpha J+s_0 f+\Lambda kb^2/N\)，所以有限 s0 覆盖小正半轴；不能将均值差 kδ/N 除以趋零的 f。

在固定紧区间 [r0,R]，零 Bernoulli 项给 f 的一致严格正下界。又 r≥0 时 \(\phi_N(r-S)\le1+S/N\)，故真值≤1+kmu/N≤1+mu。增大有限 s0 即覆盖该紧区间。

最后处理右尾。有限 Gamma–二项混合对固定 k 给一致界
\[
 c e^{-r}r^{k-1}\le f_k^a(r),\quad
 c e^{-r}r^k\le v_k^a(r)\le C e^{-r}r^k\quad(r\ge R).
\]
由 (1)，\(\alpha J_k^a/v_k^a\ge\kappa r\)，选固定 0<kappa<1。这里正确增长是**线性 r**，不是 r²；指数势未必趋于无穷，不能靠常数上界收尾。取 \(\theta=1-\kappa/2\)。对 N≥1/theta，\(\phi_N(x)\le e^{-\theta x}\) 全局成立，且统一核密度≤2e^{-z}，所以全部反馈策略满足
\[
 \mathcal U_{k,N}(r)\le[2/(1-\theta)]^k e^{-\theta r}.
\]
另一方面 \(I\ge c r^k e^{-(1-\kappa)r}\)，以严格较慢的指数衰减支配该上界。先选足够大 R，再完成上述紧区间的 s0 选择，即得 (2)。所有常数对充分大 N 一致。证毕。

右尾的精确斜率为
\[
 \lim_{r\to\infty}\frac{J_k^a(r)}{r v_k^a(r)}
 =\frac{p_a}8\frac{k+H_k}{k+1}>0,
 \qquad \lim_{r\to\infty}r f_k^a/v_k^a=k.
 \tag{4}
\]
Wolfram 对 k=4,8,16 精确核验 (4) 的多项式首系数；a=mu 时斜率约 .03312553、.03242326、.03103934。小端式 (3) 的 leading coefficient 是 \((1-1/N)/[k(k+1)]\)。

**范围。** (2) 只给固定时刻初值。第36轮局部推导对密度常数的依赖、第37轮 N^gamma 起点以及中间时段全空间比较不能由它自动消除。下面证明，现有形状在接近锐系数时确实存在中间时段反例。

## 3. 接近1的固定膨胀系数不能直接全空间延拓（严格障碍）

令 \(\Lambda_c(\eta)=\log Ee^{\eta Z_c}\)、\(\Lambda_* =\sup_c\Lambda_c\)，\(I_*\) 为 \(\Lambda_*\) 的 Legendre 率函数。正小倾斜处最优动作分支光滑。设
\[
 y(\eta)=\Lambda_\mu'(\eta),\quad
 I_\mu(y)=\eta y-\Lambda_\mu(\eta),\quad
 D(\eta)=\frac{p_\mu}8(y-\mu)^2\eta(1-\eta),\quad
 \mathscr J_\alpha=I_\mu-\alpha D.
\]
\(\mathscr J_\alpha\) 是指数化势的率函数，和 Poisson 势 J_t 不同。

**命题。** 存在 epsilon0>0：对每个固定 \(\alpha=1+\varepsilon\in(1,1+\varepsilon_0)\)，把第36轮势（可加任意固定 s0）与第34轮两侧粗屏障直接取 min，不能成为从固定时间起对所有状态有效的 supersolution。反例可取
\[
 t_N=\lfloor\sqrt{\log N}\rfloor,\qquad r_N=t_Ny(c\varepsilon),
 \quad c>A_3/K_4,
 \tag{5}
\]
这里可固定 c=6；更一般地，先固定任一 c>A3/K4，再选依赖 c 的 epsilon0。其常数为
\[
 A_3=\frac{\mu(\mu-1)}8,
 \quad K_4=\frac{\mu(\mu-1)(3+5\mu)}{384}>0,
 \quad A_3/K_4=5.110556937\ldots.
\]
密度源限定为实际使用的收敛序列
\(s_t=s_0+B\sum_{j=t_*+1}^t\log^6j/j^{3/2}\)（或相同收敛性质的固定序列）；不把结论扩张至任意有界振荡源或修改率函数的新势。

### 3.1 率与 Hamiltonian 的共同正余量

第36轮精确四阶系数给
\[
 \Delta:=\mathscr J_{1+\varepsilon}-I_*
 =-\varepsilon A_3\eta^3+K_4\eta^4+O(\varepsilon\eta^4+\eta^5).
 \tag{6}
\]
具体地 \(I_\mu-I_*=A_3\eta^3+B_{gap}\eta^4+O(\eta^5)\)，
\(D=A_3\eta^3+B_D\eta^4+O(\eta^5)\)，其中
\(B_{gap}=\mu(\mu-1)(99-43\mu)/384\)、\(B_D=A_3(2-\mu)\)，两者之差正是 K4。

取 eta=c epsilon，得到
\(\Delta=\varepsilon^4c^3(cK_4-A_3)+O(\varepsilon^5)>0\)。若 prime 表示 y 导数，则
\[
 \mathscr H_\alpha:=\mathscr J_\alpha-y\mathscr J_\alpha'
                 +\Lambda_*(\mathscr J_\alpha')
 =\Delta+\tfrac12\Lambda_*''(I_*')(\Delta')^2+o((\Delta')^2)
 =\varepsilon^4c^3(cK_4-A_3)+o(\varepsilon^4)>0.
 \tag{7}
\]
此处 \(\mathscr J_\alpha'\in(0,1)\)、\(\Delta'=O(\varepsilon^3)\)。由此同时得到 inner 率高于粗屏率、以及其单步 Hamiltonian 符号错误。只比较率函数而不检查 (7) 不足以得出反例。

### 3.2 一步比值需要的 prefactor，而非仅大偏差率

这里 eta 是先固定的正小数，再令 t→∞。在它的固定紧邻域，以 nu、kappa3 表示倾斜增量的方差、三阶累积量。直接 Fourier 展开给
\[
 \frac{f_t(ty-z)}{f_t(ty)}=e^{\eta z}
 \left[1+\frac{\kappa_3z/(2\nu^2)-z^2/(2\nu)}t
 +O(\operatorname{polylog}(t)/t^2)\right],\quad |z|\le C\log t.
 \tag{8}
\]
这里给出所需展开的直接证明，避免把裸大偏差率当作一步估计。固定正 eta 紧区间和 a 的正紧邻域。令 psi(u) 是倾斜且中心化的单步特征函数，g_t 是其和的密度。指数分量给一致界：原点邻域 |psi(u)|≤exp(−c u²)，邻域外 |psi(u)|≤rho<1，且积分 |psi|² 一致有限。Fourier 反演为
\[
 g_t(-z)=\frac1{2\pi\sqrt t}\int_{\mathbb R}e^{iwz/\sqrt t}\psi(w/\sqrt t)^t\,dw.
\]
在 |w|≤t^(1/12)，解析 cumulant 展开给
\[
 t\log\psi(w/\sqrt t)=-\nu w^2/2-i\kappa_3w^3/(6\sqrt t)
 +\kappa_4w^4/(24t)+i\kappa_5w^5/(120t^{3/2})+O(w^6/t^2).
\]
把指数乘积展开，t^(−1/2)、t^(−3/2) 的项均为奇次 w 多项式，沿对称区间积分精确为零；余项由 Gaussian 包络及多项式矩控制。对 |z|≤C log t，得
\[
 g_t(-z)=\frac1{\sqrt{2\pi t\nu}}\left[1+\frac{A_0+\kappa_3z/(2\nu^2)-z^2/(2\nu)}t
 +O((1+|z|^8)/t^2)\right],
 \quad A_0=\kappa_4/(8\nu^2)-5\kappa_3^2/(24\nu^3).
\]
中间频段由 exp(−c w²) 控制，|u| 超出固定邻域后用 rho^(t−2) 积分 |psi|²；两者都是指数小量，且位移相位模为1。具体外频积分不超过 rho^(t−2)乘以积分|psi|²，后者由倾斜 Exp 因子统一控制。除以 z=0 的式子，A0/t 抵消；再用原分布与倾斜分布的精确关系，即得 (8)。所有余项在上述紧参数集一致。

同一 Fourier L1 界给 g_t 的全局 sup≤C/sqrt(t)，中心≥c/sqrt(t)，故全局原密度比≤C exp(eta z)。因此把 |z|>C log t 的 payoff 积分尾压到任意指定幂次是合法的；C 可随固定正 eta 紧区间变化，不要求 eta 随 t 趋零。这一限定正适用于先固定 epsilon 的反例。

辅助密度的 mgf 相对原 t 步 mgf 的乘子为
\(G_1(\theta)=(1-\theta)^{-2}\)、\(G_2(\theta)=e^{\theta a}/(1-p_a+p_ae^{\theta a})^2\)。它们在 saddle 邻域解析，且沿该固定倾斜竖线一致有界：后者的分母模至少 \(1-p_a-p_ae^{\eta a}>0\)，因为 eta<1。对同一 Fourier 积分乘以 G 并展开，中心比值为 G(eta)+O(1/t)；奇次项仍消失。这严格给出下文 C_t/t 的后验展开。

对 \(\phi_\infty=\min(1,e^{-z})\) 积分得
\[
 v_t/f_t=h_0+h_1/t+o(t^{-1}),\quad
 h_0=1/[\eta(1-\eta)],\quad
 h_1=\kappa_3h_0'/(2\nu^2)-h_0''/(2\nu).
 \tag{9}
\]
这里 h0 的 prime 为 eta 导数。精确辅助密度恒等式
\[
 f_{t,a}C_{t,a}=t g_{\Gamma(t+2)+a\mathrm{Bin}(t)}
 +a^2t p_a(1-p_a)g_{\Gamma(t)+a\mathrm{Bin}(t-2)}(\cdot-a)
\]
同样给 \(C_{t,a}/t=\nu_a+O(t^{-1})\)，先在 a 的紧邻域一致成立，再取匹配 a_N 的极限。置 \(A=p_\mu(y-\mu)^2/8\)、\(B_0=p_\mu/8\)，则 (1) 给
\[
 J_t/f_t=tA+B_0\nu\log t+B_0\nu(\gamma_E-1)+o(1).
\]
所以无限终端的 inner 具有参数一致 prefactor
\[
 \log I_t(ty)=-t\mathscr J_\alpha(y)
 +[-1/2+\alpha B_0\nu/h_0]\log t+F(y)+o(1),
 \tag{10}
\]
F 光滑（包括 h1 修正与 s_infty/h0）。这里 s_t 收敛至 s_infty；仅有有界而无收敛/差分控制的源不能直接使用此式。

对 (5) 的有限 N 候选，先对 a_N 使用以上一致展开，并由方差正对隐式倾斜参数作光滑调整；所有系数在 a 上光滑。在所选时间尺度，\(t_N|a_N-\mu|\to0\)，参数改变对 log(inner) 的影响为 o(1)。由同样积分有 call/f=eta^(−2)+O(1/t)，finite-call 的 v/f 修正为 O(1/N)，在指数中的影响为 O(t/N)；而仿射 floor \(\Lambda t_Nb_N^2/N\) 相对 \(e^{-t_N\mathscr J_\alpha}\) 趋零，故同一展开仍成立。把 t−1 处斜率展开到 \(y+(y-z)/(t-1)\)，(10) 严格推出每个固定 z≥0 的比值
\[
 I_{t_N-1,N}(r_N-z)/I_{t_N,N}(r_N)
 \longrightarrow\exp\{\mathscr J_\alpha-y\mathscr J_\alpha'+z\mathscr J_\alpha'\}.
 \tag{11}
\]

### 3.3 min 不能遮住反例

沿 (5)，Ominus→1，inner→0；Oplus 的率为 I*，由 Delta>0，它指数大于 inner。因此当前点及每个有界左移的前一时间点都选择 inner 分支。选择有限动作 c 使
\(\mathscr J_\alpha-y\mathscr J_\alpha'+\Lambda_c(\mathscr J_\alpha')>0\)，它由 (7) 的 supremum 定义保证。对非负积分先限于有界 z 区间，应用 Fatou，再扩区间，得
\[
 \liminf_N K_cW_{t_N-1,N}(r_N)/W_{t_N,N}(r_N)>1.
 \tag{12}
\]
这与全局单步 supersolution 矛盾。无需对无限左移尾作支配收敛。证毕。

该反例不否定第37轮 t≥N^gamma 的全空间比较：此处 t_N远小于N^gamma，floor 尚不足以遮蔽 inner；也不否定新的状态/时间相关率函数或额外尾部修补。

## 4. Mathematica 与下界联合优化

`initialization_checks39.wl` 实际精确核验固定 k 的小端 Beta 积分、右尾密度/payoff 首项及 (4)，k=4,8,16 全部通过。

`global_inflation_check39.wl` 实际以高精度求解最优动作及 Legendre saddle，epsilon=.01,.005,.002,.001，c=5,6,7。c=6 时 Delta/epsilon^4 约1.30405、1.50034、1.60653、1.64014，Hamiltonian 同号；c=5 时均负。所有 saddle 与候选导数都在 (0,1)。基准精度80及提升至160后，归档的24位显示结果一致；初次试算的嵌套求根精度警告已通过精确有理 epsilon 及统一参数精度处理，最终两次运行无该警告。严格结论来自 (6)–(12)，这些数值不是区间证书。

下界使用 n=64、b=32、33组/65锚点，128是求积节点。保留真实 Phi64 分母及接受接收锥并集。输入权重改为
\[
 \widetilde w_j=w_j e^{\eta(v_j/16)^2}
 \frac{\sum_i w_i}{\sum_iw_i e^{\eta(v_i/16)^2}},
\]
receiver/source 同用 \(c_j=c_0+h(v_j/16)^2\)。初始 h=.075,.15,.3、eta=−.2,0,.2、scale1.025，9组64/32；按实际 log(V/J) 选择前两项，在128/32分别调scale1.020、1.030，共4次复核。归一化保持上述组权重之和，并未把它冒充实际输入的 Orlicz 成本。

最佳 h=.15、eta=.2、scale1.020：log(V/J)=71.28079768609，minAvg=1.00080796668，全33组数值通过。仍低于同维标量参考71.71499208795约.43419440，未登记新认证下界。进展指标64logB−log(V/J)约10.54093，同维标量约10.10674；越小越好。负eta及h=.3未改善实际评分。

## 5. 复现、来源与下一轮

从项目根目录运行：

```sh
wolframscript -file code/initialization_checks39.wl /tmp/init39.json
wolframscript -file code/global_inflation_check39.wl /tmp/rate39.json
wolframscript -file code/global_inflation_check39.wl /tmp/rate39-high.json 160
wolframscript -file code/lower_reweight39.wl /tmp/lower39.json
```

新增脚本通过现有第38轮 lower kernel 的 library guard 复用定义，避免重跑上一轮17项扫描。无随机输入；数值扫描非认证。只归档本报告、必要脚本与最小摘要，更新当前总览；没有需要移除的冗余跟踪材料。

来源索引：本项目 r015 的指数目标优化、r016/029 的倾斜 Fourier 方法、r034 两侧粗屏障、r035/036 精确四阶系数与辅助密度恒等式、r037/038 晚时段比较及初始化障碍。本文固定时刻/中间时段机制结论的证明已在上文给出；不作未经文献核实的公开优先权声明。

下一轮两项动作：① 根据 (6)–(7) 修改正倾斜区域的率函数，使其在管内保留锐系数、在中间区满足正确 Hamiltonian 符号；不能继续仅增大有限密度常数。② 下界重新优化完整source权重或接收曲线形状，同时重算真实成本及并集体积；若始终低于标量，转向更具区分力的非对称/增长块结构。
