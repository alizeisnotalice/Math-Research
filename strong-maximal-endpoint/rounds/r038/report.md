# 第38轮：统一初始化面积、精确控制损失与移位障碍

日期：2026-11-06。状态：初始化面积适用范围得到严格扩展；证明全局参考移位的线性代价障碍，并给出保留反馈与损失的误差传播接口。一般锐对数系数尚未闭合。本轮使用本机 Wolfram 15，未启动 MATLAB。

## 1. 不变的一般账本

统一 \(\Phi_n(t)=t(1+\log^+t)^{n-1}\)，\(\mu=1+e^{-\mu}=\log B\)，\(\sigma=\sqrt\mu\)，\(\bar\Phi_G(z_n)=\sigma z_n/\sqrt n\)，\(d_n=\sigma\sqrt n z_n\)。严格结果仍为
\[
 \log Q_n=n\mu-d_n+\frac{9-\mu}{24}\log n+o(\log n),\quad
 Q_n\le C_n\le\mathcal R_n,\quad
 \log C_n=n\mu-d_n+O(\log n).
\]
第29轮的额外 \(g\log n\)、\(g=(\mu-1)/8\)，属于上界模型 \(\mathcal R_n\) 的下界，不能转写为真实 \(C_n\) 下界。

沿用 payoff \(\phi_N=e^{-r}\)（r>0）、\(1-r/N\)（r≤0），\(T=\sup_{c\ge0}K_c\)、\(K_ch(r)=Eh(r-Z_c)\)，\(Z_c=\operatorname{Exp}(1)+c\operatorname{Bern}(p_c)\)、\(p_c=(1+e^c)^{-1}\)。真值 \(\mathcal U_{t,N}=T^t\phi_N\)，固定动作参考 \(v_{t,N}^a=K_a^t\phi_N\)。用 \(m(a)=EZ_a\)、\(\nu(a)=\operatorname{Var}Z_a\)。第35–37轮匹配动作满足
\[
 |a-\mu|=O(\sqrt{\log N/N}),\quad
 \delta=\mu-m(a)=O(\log N/N),\quad b=m(a)-a,\quad b^2\asymp\log N/N.
\]
第37轮全局晚时段势为 \(W_t=\min(O_t,I_t)\)，其中 \(O_t=\min(O_t^-,O_t^+)\) 是全局粗屏障，\(I_t\) 是匹配非线性势。对固定 \(0<\gamma<1\)，从 \(t\ge\lceil N^\gamma\rceil\) 有全空间 \(TW_{t-1}\le W_t\)。仍需真正初始化，不可直接把面积误差按参考桥传播。

## 2. 匹配参考面积：从 gamma<1/2 扩展至所有 gamma<1（严格）

第37轮用逐步耦合得到 \(CT|a-\mu|\) 误差。本轮以 iid 和的 Gaussian 耦合替换它，证明对 \(T\le N\)、a 位于 mu 的固定紧邻域，
\[
 \boxed{\int(\mathcal U_{T,N}-v_{T,N}^a-T\delta/N)_+\,dr
 \le C\log(2+T)+C+T\delta+C\sqrt T|a-\mu|+CT/N.}
 \tag{1}
\]
因此 \(T=\lceil N^\gamma\rceil\) 时，对任意固定 \(0<\gamma<1\)，后三个显式余项均趋于0，而原逐步耦合只对 \(\gamma<1/2\) 给出这一结论。

**证明。** 令 \(G_T^c\sim N(Tm(c),T\nu(c))\)。Gamma–Bernoulli 增量在 c 的固定紧邻域具有统一有限四阶矩、方差下界1，标准化 iid 和的 Wasserstein–1 CLT 给出
\[
 W_1(S_T^c,G_T^c)\le C
 \tag{2}
\]
（未标准化尺度，常数对该紧集一致）。用同一个标准正态耦合，
\[
 W_1(G_T^\mu,G_T^a)
 \le T\delta+\sqrt{2T/\pi}\,|\sqrt\mu-\sqrt{\nu(a)}|
 \le T\delta+C\sqrt T|a-\mu|.
 \tag{3}
\]
故 \(W_1(S_T^\mu,S_T^a)\le C+T\delta+C\sqrt T|a-\mu|\)。局部光滑性直接由 \(\nu(c)=1+c^2p_c(1-p_c)\) 给出，且
\(m''(\mu)=-(\mu-1)/\mu\)、\(\nu(\mu)=\mu\)、\(\nu'(\mu)=\mu-1\)。

第37轮已由第21轮全状态 CDF 屏障、正 mgf 尾截断证明
\(\int(U_T^\infty-v_T^{\mu,\infty})_+\le C\log(2+T)\)。对 \(\phi_\infty(r)=\min(1,e^{-r})\)，加独立 Exp 后的 survival 面积差不超过 \(W_1\)，故以 (2)–(3) 替换第37轮逐步耦合项即可。

有限 N 部分保留精确分解 \(\phi_N(r-S)=\phi_\infty(r-S)+(S-r)_+/N\)。对任意反馈策略，中心鞅 Bernstein 尾与 \(S_T\le\mu T+M_T\) 给 call 正面积差 \(CT\)；固定参考的均值补偿用
\[
 C_{S_T^a}+T\delta\ge C_{S_T^a+T\delta},\qquad
 \int|C_X-C_Y|\le(\operatorname{Var}X+\operatorname{Var}Y)/2
\]
（X,Y 同均值），仍仅贡献 \(CT/N\)。这正是第37轮 call 模块，不需新的二阶 transport 假设。合并得到 (1)。

匹配尺度给
\(T\delta=O(N^{\gamma-1}\log N)\)、\(\sqrt T|a-\mu|=O(N^{(\gamma-1)/2}\sqrt{\log N})\)，均对固定 gamma<1 趋零。证毕。

这是真实适用范围改进，但仍是积分结论，没有点态初始化或锐端点上界的含义。

## 3. 裸全局参考移位必须付出线性代价（严格障碍）

不能直接用 \(O(\log T)\) 的固定空间移位把真值全实轴压在固定参考之下。

**命题。** 固定整数 T≥1、有限 N≥T、a≥0。若无任何加法尾补丁地有
\[
 \mathcal U_{T,N}(r)\le v_{T,N}^a(r-H)\quad\text{对所有实数 r},
 \tag{4}
\]
则必有
\[
 \boxed{H\ge T\log(1+e^{-a}).}
 \tag{5}
\]
对 a=mu 即 \(H\ge T\log\mu\)，排除 \(H=O(\log T)\) 的全局方案。

**证明。** 对任何固定常值动作 c，Gamma–二项卷积在 r→∞ 的最高项为
\[
 v_{T,N}^c(r)
 =\frac{e^{-r}r^T}{T!}
 \left(\frac2{1+e^{-c}}\right)^T[1+O_{T,c,N}(r^{-1})].
 \tag{6}
\]
条件于 Bernoulli 位移和 B，指数支贡献 \(e^{-r+B}(r-B)^T/T!\)。CDF 支与有限 N 的 call 支仅有最高 T−1 次多项式；call 中的 T 次项精确抵消。再取 \(Ee^B=[2/(1+e^{-c})]^T\) 得 (6)。由真值≥任意常值参考，比较 (4) 两侧的 (6) 首项，再令固定 c→∞，得
\(e^H[2/(1+e^{-a})]^T\ge2^T\)，即 (5)。先 r→∞ 再 c→∞，没有交换非一致极限。证毕。

此障碍不排除“管内参考移位、管外使用 O”的拼接，也不排除加法尾补丁。在第37轮晚时段外管，H=O(log T) 的移位仍留在稍窄外管，因此 \(W_T(r-H)=O_T(r-H)\ge O_T(r)\) 可用；真正缺口在管内的相对比较及常数依赖。

特别地，外管宽度 A 满足 A² 与 1/gamma 同阶时，初始化移位 \(C_A\log T\) 的系数可能依赖 A²；不能据 \(\log T=\gamma\log N\) 就宣称令 gamma→0 会消除全部代价。第21轮 Gaussian 误差指数的提高也会改变初始化常数。

## 4. 保留反馈损失的精确传播接口（严格）

记晚时段每个动作的非负损失
\[
 R_t^c(r)=W_t(r)-K_cW_{t-1}(r)\ge0.
\]
沿任意可预测反向策略、\(r_t=q-\sum_{j=t+1}^N Z_j\)，条件期望逐项相消，给
\[
 W_N(q)-E_\pi W_T(r_T)=\sum_{t=T+1}^N E_\pi R_t^{c_t}(r_t).
 \tag{7}
\]
因此初始化误差的收益必须与同一策略的基准损失一起结算。

令 \(E_T=(\mathcal U_{T,N}-W_T)_+\)，定义
\[
 e_T=E_T,\qquad
 e_t(r)=\sup_c[K_ce_{t-1}(r)-R_t^c(r)]_+
       =[T(W_{t-1}+e_{t-1})(r)-W_t(r)]_+.
 \tag{8}
\]
这是满足所有 \(K_ce_{t-1}\le e_t+R_t^c\) 的最小非负递推。由 (7) 或逐步比较，
\[
 \mathcal U_{N,N}(q)\le W_N(q)+e_N(q).
 \tag{9}
\]
下一真正所需引理是对**实际结构误差**证明 \(e_N(q)\lesssim f_N^a(q)\log T\)。第37轮已否定无损失的任意 L1 传播；(8) 没有把控制核替换为参考桥，但也尚未完成 (9) 中 e_N 的锐估计。

### 4.1 左外侧的定量损失

在晚时段左外管，实际选择分支为 \(W_t=O_t^-\)。理由：置 \(x=t\mu-r>0\)，Jensen 给 \(O_t^+\ge e^{x/N}\)，而 \(O_t^-=1+x/N+B_t^\Gamma\)。若 r≤rho_N，B=0；否则 r−rho_N≤mu t−N+1≤(mu−1+o(1))t，Gamma 小尾给 \(B\le Ce^{-ct}\)。另一方面 \(e^{x/N}-1-x/N\ge x^2/(2N^2)\gtrsim t\log t/N^2\)，在 t≥N^gamma 上压过该指数小量。故 O^-≤O^+，结合第37轮外管拼接得到所述分支。

Gamma 修补非减，且共同耦合下 Z_c≥Exp(1)，所以
\[
 R_t^c(r)\ge O_t^-(r)-K_cO_{t-1}^-(r)
 \ge[\mu-m(c)]/N.
 \tag{10}
\]
这给所有反馈策略的左外管累计均值损失杀死项。它只在该外管使用，不能替代管内损失估计。

### 4.2 可验证的相对粗屏缺口条件

由于 T 为正线性动作的上确界，对固定 \(0\le\lambda\le1\)，\((1-\lambda)W_t+\lambda O_t\) 也是 supersolution。若能验证初始化
\[
 E_T\le\lambda(O_T-W_T),
 \tag{11}
\]
则归纳给 \(e_t\le\lambda(O_t-W_t)\)。精确代数为
\[
 K_c[\lambda(O_{t-1}-W_{t-1})]-R_{W,t}^c
 =\lambda(O_t-W_t)-(1-\lambda)R_{W,t}^c-\lambda R_{O,t}^c.
\]
此条件保留控制损失，可作为具体认证目标。临界点粗屏缺口至多 O(d_N/N)，故欲保持 o(f_N log N)，一个充分尺度为 \(\lambda=o(\log^{3/2}N/\sqrt N)\)。仅有 (1) 不会推出 (11)，本轮没有声称该所需 lambda 已获证明。

## 5. Mathematica 实际运行与证据边界

`area_checks38.wl`：精确核验均值二阶导、方差一阶导、仿射损失及凸缺口恒等式，尾多项式 T=1,2,4,8 的次数/首系数全部通过。对 N=2^16,2^24,2^32,2^40、gamma=.75，以60位及100位独立求根计算匹配动作。共同 Gaussian 耦合项分别约 .0439655、.0270975、.0152588、.00835094；旧逐步尺度 T|a−mu| 则约20.97、113.36、553.28、2562.95。精度提升差在50位以下。Gaussian 数值项不包含 (2) 的统一 O(1) CLT 常数，不能将该小数直接当作全部 W1 误差。

`initialization_probe38.wl`：实际从终端开始执行完整 adaptive payoff Bellman 递推，N=10^8、t至128，h=.05/.025；在 t=16,64,128 比较固定 mu 参考的指数/双侧校正、alpha=1.1、beta=0/.01。归档管宽 A=3 的24组摘要（其余嵌套小管不另存）。beta 源为探索性 log^4 可求和项，与证明使用的匹配参数、仿射余量和 log^6 大起点项不同。

h=.025、t=16、beta=0 时，指数/双侧的最大初值正缺口约 .000413400/.000551133，最大“正缺口/粗屏缺口”约 .00379156/.00684684。beta=.01 后，这个管内正缺口降到机器舍入量级；t=64,128 的这些探针也未见正缺口。该结果仅提示有限时间初始化值得认证，不是点态证明、渐近 lambda 界或全状态证书。分段线性插值误差、机器尾下溢及粗屏缺口≤1e-10处的比值截断都未认证；最终脚本记录递推与候选的单调性检查，两种网格全部为0违例。

下界 `lower_scan38.wl`：n=64、b=32、33组/65锚点；128指求积节点。source 与 receiver 使用同一组 \(c_j=c_0+h(v_j/16)^2\)，同时扫描 h=−.10,−.025,.025,.075,.15 及 source scale=1.025,1.030,1.035，共15组64/32求积，按实际 Phi64 分母与接受锥并集评分挑选，再对两项作128/32复核。最佳 h=.15、scale1.025 的评分由71.27603027776变为71.27603061284，全部33组通过数值平均，minAvg=1.00193923725。临界候选 h=.025、scale1.030 的 minAvg=1.00005026029、评分71.26996558308。

最佳值略高于旧曲线候选，仍低于同维标量参考71.71499208795。评分为 log(V/J)，最佳候选的缺口64logB−log(V/J)约10.54570012，同维标量参考缺口约10.10673865（越小越好）。所有结果为浮点诊断，未增加认证下界；并未据固定维数33模板试验宣称实现渐近增长模板机制。

## 6. 复现、来源与下一轮

从项目根目录执行：

```sh
wolframscript -file code/area_checks38.wl /tmp/area38.json
wolframscript -file code/initialization_probe38.wl /tmp/init38.json
wolframscript -file code/lower_scan38.wl /tmp/lower38.json
```

代码调用现有第32、35、36轮必要库与第30轮 seed，不依赖仓库外路径；无随机输入。本轮仅归档最终报告、3个脚本和对应最小 JSON 摘要，以及更新当前总览。MATLAB 未启动，Wolfram CLI 批次结束退出。

来源索引：本项目 r021 全状态 Gaussian CDF 屏障、r029 倾斜 Fourier 方法、r032 payoff 与预算接口、r034 粗屏障、r036 局部非线性势、r037 晚时段全空间比较与有限 call 面积证明。统一 W1-CLT 使用 Nathan Ross, *Fundamentals of Stein’s method*, Probability Surveys 8 (2011), 210–293, DOI10.1214/11-PS182，[原始预印本 Theorem3.2，第10页](https://arxiv.org/pdf/1109.1880#page=10)，有限四阶矩版本；紧参数族矩界与方差下界确保常数一致。本轮无新的公开优先权声明。

下一轮两项动作：① 对实际初始化误差验证 (11) 或直接控制 (8)，优先检查固定短时校正的可认证初始化及保留管宽常数的倾斜比较；不能以裸全局平移或面积传播替代。② 在已全33组数值通过的 h=.15 附近调整 source 权重、scale及比例曲线，按实际评分选择；只对真正胜过同维标量的候选启动完整下界认证。
