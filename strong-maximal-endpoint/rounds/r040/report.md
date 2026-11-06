# 第40轮：正侧再进入的局部剪切、率修补与屏障交接

日期：2026-11-06。状态：证明正侧局部剪切结构、修补后的率函数符号及其固定时刻初始化，并确定固定 mu 模型的交接尺度；有限时间全动作中间带仍开放。实际计算使用本机 Wolfram 15，未启动 MATLAB。

## 1. 一般账本与本轮结果

统一 \(\Phi_n(t)=t(1+\log^+t)^{n-1}\)，\(\mu=1+e^{-\mu}=\log B\)，\(\sigma=\sqrt\mu\)，\(\bar\Phi_G(z_n)=\sigma z_n/\sqrt n\)，\(d_n=\sigma\sqrt n z_n\)。严格结果仍为
\[
 \log Q_n=n\mu-d_n+\frac{9-\mu}{24}\log n+o(\log n),\qquad
 Q_n\le C_n\le\mathcal R_n,\quad \log C_n=n\mu-d_n+O(\log n).
\]
第29轮 g=(mu−1)/8 的额外收益属于反馈上界模型，不能登记为真实 C_n 下界。本轮没有改变一般上下界。

上一轮证明旧二次势在固定小正倾斜区会重新低于粗屏障，导致中间时间反例。本轮得到两种处理方式：一是利用内部分支已高于粗屏障的带，强制在更右侧保留粗屏障，排除再进入；二是修改正侧率函数。第一种结构更简单，且数值计算未显示第二种修补能普遍降低有限时间残差。因此下一步优先推进局部剪切留下的中间带，不把正确率符号当作全局单步证明。

记 \(T=\sup_{c\ge0}K_c\)、\(K_ch(r)=Eh(r-Z_c)\)，\(Z_c=\operatorname{Exp}(1)+c\operatorname{Bern}(p_c)\)，\(p_c=(1+e^c)^{-1}\)。终端 \(\phi_N=e^{-r}\)（r>0）、\(1-r/N\)（r≤0），真值 \(\mathcal U_{t,N}=T^t\phi_N\)。匹配 a_N→mu，\(b=m(a)-a\)、\(\ell=2a-m(a)\)、\(b^2\asymp\log N/N\)。旧势为
\[
 I_{t,N}=v_t^a\exp[(\alpha J_t^a+s_tf_t^a)/v_t^a]+\Lambda tb^2/N,
 \quad \alpha=1+\varepsilon>1,
\]
其中 s_t 是固定有界收敛非负密度系数；精确 Poisson 公式见 r036。全局粗屏障为 \(O_t=\min(O_t^-,O_t^+)\)。

## 2. 正侧局部剪切排除再进入（严格结构引理）

令 \(\Lambda_c(\eta)=\log Ee^{\eta Z_c}\)、\(\Lambda_*=\sup_c\Lambda_c\)，\(I_*\) 为其 Legendre 率函数。以 \(y(\eta)=\Lambda_\mu'(\eta)\) 参数化正侧，r/t≤mu 时把 eta 记为0。

旧势的固定 mu 率为 \(\mathscr J_\alpha=I_\mu-\alpha D\)，其中
\[
 D=\frac{p_\mu}8(y-\mu)^2\eta(1-\eta),\quad
 A_3=\frac{\mu(\mu-1)}8,
 \quad K_4=\frac{\mu(\mu-1)(3+5\mu)}{384}.
\]
r039 给 \(\mathscr J_{1+\varepsilon}-I_*=-\varepsilon A_3\eta^3+K_4\eta^4+O(\varepsilon\eta^4+\eta^5)\)。先取固定 d>0 充分小，再取 epsilon 充分小，使整个带
\[
 d\varepsilon\le\eta\le2d\varepsilon
 \quad\Longrightarrow\quad
 \mathscr J_\alpha-I_*\le-c\varepsilon\eta^3<0.
 \tag{1}
\]
这与上一轮反例 eta=c epsilon、c>A3/K4 所处区域不同。

对所有充分大匹配 N，以及 \(t_0(\varepsilon)\le t\le N\)，在该固定带上一致有
\[
 I_{t,N}(ty(\eta))\ge O_t^+(ty(\eta)).
 \tag{2}
\]
证明使用 r039 已证明的固定倾斜 Fourier prefactor：log(inner)/t=−mathscrJ+O(log t/t)，而 Oplus=e^(−tI*)。匹配 a_N 的率连续性将 (1) 的严格余量保留；无需错误地假设所有 t≤N 都有 t|a_N−mu|→0。有限 call 对 log(inner) 的修正 O(t/N)，即率修正 O(1/N)，同样可由先选充分大 N 吸收。正 floor 只会增加 inner。源系数非负且有界。

定义局部剪切
\[
 W^{loc}_{t,N}(r)=
 \begin{cases}
 \min(O_t(r),I_{t,N}(r)),&\eta(r/t)\le2d\varepsilon,\\
 O_t(r),&\eta(r/t)>2d\varepsilon.
 \end{cases}
 \tag{3}
\]
由 (2)，整个 eta≥d epsilon 的区域都已是粗屏障：在带内 min 本来就选择 O，带外强制选择 O。接口附近无值跳跃。

由正性及全局 \(TO_{t-1}\le O_t\)，当前 eta≥d epsilon 时立即有
\(TW^{loc}_{t-1}\le O_t=W^{loc}_t\)。在 eta<d epsilon 的当前点，所有可达前态都是 r−z≤r；时间减一最多使正倾斜参数增加 O(1/t)，因此对充分大 t 不会到达 eta>2d epsilon 的强制区。故对所有可达前态仍有 \(W^{loc}_{t-1}\le I_{t-1}\)。**任何已经证明的当前点内分支不等式 \(TI_{t-1}\le I_t\)，都可与全局 O 比较同时使用，转移到 (3)**。

此引理移除了正侧固定倾斜尾部再进入，无须修改内部势的公式。它保留 r036 自然管的单步比较及 r039 固定 k 初始化（因为只将旧 min 向上放宽）。但它没有证明自然管与 eta<d epsilon 之间的正中偏差单步比较，也没有证明左侧中间剪切；这些仍是全局闭合的实质缺口。固定带起效时间可很大，不能把 (2) 当作 t=16 或128 已成立的数值事实。

## 3. 率函数修补：可用的替代方案（严格率结论）

设
\[
 E(\eta)=I_\mu(y)-I_*(y),\quad F=E-D=K_4\eta^4+O(\eta^5).
\]
对正侧给旧指数势的 log 加 tF；负侧延为0。因 F=O(eta^4)，该延拓为 C³。新的率函数是
\[
 \mathscr J_{new}=I_*-\varepsilon D.
 \tag{4}
\]
令 prime 表示 y 导数。由 Legendre 等式及 Taylor 的积分余项，
\[
 \mathscr H_{new}=\mathscr J_{new}-y\mathscr J_{new}'+\Lambda_*(\mathscr J_{new}')
 =-\varepsilon D+O(\varepsilon^2(D')^2)
 =-\varepsilon D+O(\varepsilon^2\eta^4)
 \le-c\varepsilon\eta^3<0
 \tag{5}
\]
在小固定 eta 区间内一致成立，只需 epsilon·eta 足够小。这消除了 r039 的特定率反例。自然边界 eta=A sqrt(log t/t) 处，tF=O(log²t/t)=o(1)，不会改变临界左侧的 g log N 候选系数。

对每个固定正 eta，修补 inner/Oplus 的指数因子为 exp(t epsilon D)，故粗屏障最终会被选择；它不再发生 r039 的固定正倾斜再进入。这里只是率与点态分支选择结论，尚非有限时间全动作比较。

### 3.1 修补后仍可作固定时刻初始化

该修改不会破坏 r039 的固定 k≥4 初始化存在性。近 r=0 和左侧根本未修改；紧区间可增大有限 s0 补偿 tF 的负部。右尾需要检查真实斜率，不能只凭小 eta 展开。

当 eta↑1，\(M_\mu(\eta)\sim(2/\mu)/(1-\eta)\)，\(M_*(\eta)\sim2/(1-\eta)\)。后者由 Bernoulli 因子≤2、取 c=(1−eta)^(−1/2) 达到极限2得到。Legendre 的常数项可用严格夹逼：M*≤2/(1−eta) 给 I*≥y−1−log y−log2；固定任意 c 的 I_c 给相反 limsup，再令 c→∞。因此
\[
 I_\mu(y)-I_*(y)\to\log\mu,\quad D/y\to p_\mu/8,
 \quad F/y\to-p_\mu/8.
\]
固定 k 的旧 Poisson 尾斜率 \(q_k=(p_\mu/8)(k+H_k)/(k+1)\)。修补后相对 e^(−r) 的增强率为
\(\alpha q_k-p_\mu/8>0\)，对 alpha≥1、k≥4 成立，并在 a_N→mu 后保留统一正余量。继续用 theta 指数矩屏障，按 r039 的小端、紧区间、右尾三段证明，即得修补势的固定 k 全空间初始化。它依旧不是从固定 k 以后全空间单步不等式。

## 4. 屏障交接尺度的精确常数（严格受限模型）

此节仅取固定 mu、无限终端 phi_infty=min(1,e^(−r))，不把其位移估计冒充所有匹配 N 的 transport 结论。令 s_t 有界，epsilon>0 固定，
\[
 \eta=C(\log t/t)^{1/3},\quad y=\Lambda_\mu'(\eta),\quad
 \widetilde V_t=v_t\exp[(\alpha J_t+s_tf_t)/v_t+tF(\eta)].
\]
则
\[
 \boxed{\log\frac{\widetilde V_t(ty)}{e^{-tI_*(y)}}
 =\left(\varepsilon A_3C^3-\frac16\right)\log t
 -\frac13\log\log t-\log C-\frac12\log(2\pi\mu)+o(1).}
 \tag{6}
\]
临界常数为 \(C_{crit}=(6\varepsilon A_3)^{-1/3}\)。C<Ccrit 时比值趋零，C>Ccrit 时趋无穷，等号时仍因负的 loglog 项趋零。1/6 来自 eta^(−1) payoff prefactor，不能误写成只看 t^(−1/2) 的1/2。

**所需统一 Mills 主项及证明。** 倾斜中心密度记 g_t。r039 的一致 Fourier 分割给 g_t(0)≥c/sqrt(t)、sup|g_t''|≤C t^(−3/2)。在中心导数积分中，Gaussian 的奇次项精确为0，保留 cubic cumulant 后得到 g_t'(0)=O(t^(−3/2))；绝对 Fourier 估计本身不够。全实轴 Taylor 因而给
\[
 |g_t(u)/g_t(0)-1|\le C(|u|+u^2)/t.
\]
对 phi_infty(z)e^(eta z) 积分，得到在 0<eta≤delta、带显示奇异余项的一致界
\[
 v_t/f_t=h_0\left[1+O\left(\frac1{t\eta^2}+\frac1{t\eta}+\frac1t\right)\right],
 \qquad h_0=[\eta(1-\eta)]^{-1}.
 \tag{7}
\]
不声称在 eta=0 有有限 h0。这个证明不用将主质量截于 O(log t) 的错误做法。

在 (6) 的尺度，均匀 saddle 密度给 f_t=e^(−tI_mu)/sqrt(2pi t mu)·(1+o(1))，辅助密度恒等式给 C_t/t=O(1)。因此
\(J_t/v_t=tD+O(eta log t)+O(eta)=tD+o(1)\)，Mills 误差的相位影响为 tD/(t eta²)=O(eta)。又 s_t f/v=O(eta)、tF=O(t eta⁴)=o(1)。合并 log v、t epsilon D=epsilon A3 C³ log t+o(1) 及 −log eta 即得 (6)。此处只需 (7) 的主项精度；它远弱于全动作单步 transport 要求。

## 5. 有限时间缺口与实际单步探针

固定 eta 的 r039 截断不能直接推广到自然边界 eta0=A sqrt(log t/t)：
\(e^{-eta_0 C\log t}\to1\)，主质量位移约1/eta0，远大于 log t。正态代理已经给
\[
 \int_0^\infty e^{-\eta u-u^2/(2\nu t)}du
 =\eta^{-1}-(\nu t\eta^3)^{-1}+3/(\nu^2t^2\eta^5)-\cdots.
]
所以 h0+h1/t 的余项不是对 eta 统一的 O(t^(−2))。(7) 解决了交接尺度所需主项，但仍需对 exact-reference normalized ratio 的时间、状态、动作导数逐项控制，才能把 generic 1/t 改善为 eta/t。t约N时也不能把 a_N 直接换成 mu；必须保留匹配参数与有限 call 项。尚未得到覆盖整个中间带的严格单步比较。

`corrected_probe40.wl` 实际使用完整 PWL 单步 MAX 核，N=256，t=16,64,128、alpha1.01、探索源 beta=0/.01，h=.05/.025，比较旧势和加 tF 的势，两者均与左右粗屏障取 min，共24组。另 t128、h=.0125 细化4组，并以同一细 rate 表在 h=.025 单独复核4组，区分状态网格与 rate 表变化。这些探针固定参考 mu，未实现全部匹配参数及仿射余量；所有候选单调性检查均无违例。

beta=0、h=.025 的全域最大残差约：t16为2.688e-5（两者相同），t64为7.142e-7（相同），t128旧势1.273e-7、修补势1.554e-7。h=.0125时后两值分别1.23609e-7、1.23678e-7，差别显著缩小。同 h=.025 下单独细化 rate 表，最大残差改变不足4e-15，表明该次粗网格增幅主要不是 rate 表误差。由此不能把 coarse-grid 残差增幅直接解释为严格数学反例；数值结果也不支持“率修补会普遍改善有限时间核”。beta=.01 后这些采样正残差均在浮点舍入量级。该源是探索性 log^4 可求和项，不等同于证明中的匹配 log^6 大起点源。

rate 表以50/80位计算，求根 AccuracyGoal=PrecisionGoal=30，eta网格 .004/.002、三次插值；追加4个非网格 saddle 直接求根检查，细表 F 绝对插值误差最大约1.41e-15、Istar约4.37e-12（只覆盖这4点）。最终求根未出现收敛失败消息。机器混合小尾仍可能下溢；PWL、rate 插值、动作根及状态范围未经区间认证。上述均为诊断。

`asymptotic_checks40.wl`：80及160基准精度实际验证8组新率符号，归档显示的率及Mills数值逐位一致，并核验 Mills 展开、截断质量极限1及 Bernoulli 尾因子2。epsilon=.01、eta=.06 时 F≈1.135377164e-7，Hnew≈−1.004879811e-7。正态 Mills 代理在 t=10^6,10^9,10^12 的缩放误差约−.78212695、−.78218361、−.78218791，趋于−1/mu≈−.78218829；该代理不是 Gamma–二项实际值。

## 6. 下界改换曲线几何（数值；未认证）

使用 n=64、b=32、33组/65锚点，128指求积节点。保留 r039 最佳权重归一化、h=.15、权重参数eta_w=.2、scale1.020，同步将 source 和 receiver 坐标改为
\(d_j\mapsto d_j+q(v_j/16)^p\)，p=2或4，q=±.08；q=0仅算一次，共5个64/32求积。按真实 Phi64 分母及接受锥并集评分，选两项128/32复核。

|几何|128/32 log(V/J)|接受组数|
|---|---:|---:|
|负二次 q=−.08|71.29433994|31/33|
|负四次 q=−.08|71.28834203|32/33|

负二次相对旧基线提高约.01354：接受锥并集体积增加约1.526%，超过成本的小幅增加，但同时失去两组接收点。未全锚通过不等于剪去失败锥后的数值见证无效；它仍未认证，也仍低于同维标量71.71499209约.42065。主指标64logB−log(V/J)约10.52739，同维标量10.10674，越小越好。正向形变未改善评分。此结果揭示本地方向的体积—成本取舍，没有实现一般渐近收益。

## 7. 复现、来源与下一轮

从项目根目录运行：

```sh
wolframscript -file code/asymptotic_checks40.wl /tmp/asym40.json
wolframscript -file code/asymptotic_checks40.wl /tmp/asym40-high.json 160
wolframscript -file code/corrected_probe40.wl /tmp/probe40.json
wolframscript -file code/corrected_probe40.wl /tmp/refine40.json refine
wolframscript -file code/corrected_probe40.wl /tmp/mesh40.json meshcheck
wolframscript -file code/lower_geometry40.wl /tmp/lower40.json
```

新增 r039 rate 库模式避免重复旧扫描，默认旧脚本行为不变。新代码仅依赖仓库内必要库和 r030 seed，无随机输入。本轮归档报告、3个必要脚本和最小数值摘要，更新当前总览；不归档工作草稿或重算缓存。

来源索引：r015/r034 的指数矩优化与粗屏障，r036 的精确 Poisson 势，r039 的四阶系数、固定 k 初始化及 Fourier 余项；本轮各新机制的证明见上文。没有新增公开优先权声明。

下一轮两项动作：① 以 (7) 为起点证明匹配参数下的归一化时间/状态/动作导数估计，闭合正侧剩余中间带，并与左侧剪切一同审查；局部化可避免处理已被粗屏接管的远正尾。② 下界在几何获益方向重新调整输入或改用非对称块结构，继续结算真实分母和接受并集；只有胜过同维标量的候选才进入完整下界认证。
