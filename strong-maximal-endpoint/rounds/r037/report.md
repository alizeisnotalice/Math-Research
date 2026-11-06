# 第37轮：晚时段全空间比较与初始化面积

日期：2026-11-06。状态：本轮两个上界引理及面积传播障碍已证明并经独立审查；一般锐对数系数尚未闭合。实际运行以本机 Wolfram 15 为主，本轮未启动 MATLAB。

## 1. 账本与本轮增量

沿用 \(\Phi_n(t)=t(1+\log^+t)^{n-1}\)。令
\(\mu=1+e^{-\mu}=\log B\)、\(\sigma=\sqrt\mu\)，并以
\(\bar\Phi_G(z_n)=\sigma z_n/\sqrt n\)、\(d_n=\sigma\sqrt n z_n\) 定义临界缺口。目前严格结论仍为
\[
 \log Q_n=n\mu-d_n+\frac{9-\mu}{24}\log n+o(\log n),\qquad
 Q_n\le C_n\le\mathcal R_n,\quad
 \log C_n=n\mu-d_n+O(\log n).
\]
第29轮的 \(g=(\mu-1)/8\) 改善是反馈上界模型 \(\mathcal R_n\) 的下界，不是真实 \(C_n\) 的新下界。本轮没有改变这份一般账本。

本轮推进：第36轮的自然时间管局部非线性势，与第34轮两侧全局屏障取最小后，在 \(t\ge N^\gamma\) 的全部状态上满足 Bellman 单步比较。另证明此初始层真值相对补偿参考的正缺口面积为 \(O(\log t)\)。这两个结论尚不能直接拼成端点上界：仍需把面积误差传至终端的有效估计。

## 2. 记号与晚时段全空间比较（严格）

终端与算子为
\[
 \phi_N(r)=\begin{cases}e^{-r},&r>0,\\1-r/N,&r\le0,\end{cases}
 \quad K_ch(r)=E h(r-Z_c),\quad T=\sup_{c\ge0}K_c,
 \quad Z_c=\operatorname{Exp}(1)+c\operatorname{Bern}((1+e^c)^{-1}).
\]
记真值 \(\mathcal U_{t,N}=T^t\phi_N\)，固定动作参考 \(v_t^a=K_a^t\phi_N\)，
\(m(a)=EZ_a\)。使用第35–36轮匹配动作 \(a=a_N\)：
\[
 b=m(a)-a,\qquad b^2\sim\frac{\mu\log N}{4N},\qquad
 \delta=\mu-m(a)=O(b^2),\qquad |a-\mu|=O(\sqrt{\log N/N}).
\]
第36轮证明，对任意固定管宽 \(A\)、\(\alpha>1\)，可选择固定常数和充分大的固定起点，使
\[
 I_t=v_t^a\exp(U_t/v_t^a)+\Lambda tb^2/N,\quad
 U_t=\alpha J_t+s_tf_t^a,\quad
 s_t=B_0\sum_{j>t_*}^t\frac{\log^6j}{j^{3/2}},
 \tag{1}
\]
在当前状态 \(|r-tm(a)|\le A\sqrt{t\log t}\) 满足 \(TI_{t-1}\le I_t\)。其中 \(J_t\ge0\) 为第36轮精确 Poisson 校正，\(s_t\le1\)。
第34轮的全局 MAX 屏障记为 \(O_t^-,O_t^+\)，置
\[
 O_t=\min(O_t^-,O_t^+),\qquad W_t=\min(O_t,I_t).
 \tag{2}
\]
最小值保持 supersolution 的理由是正性：若 \(TW_{t-1}\le TO_{t-1}\le O_t\) 且 \(TW_{t-1}\le TI_{t-1}\le I_t\)，即可取最小值；不需要 \(T\) 与最小值交换。

**命题。** 对任意固定 \(\gamma\in(0,1)\)，适当依次选择固定 \(A,B_0,\Lambda,t_*\)，所有充分大 \(N\) 和整数 \(\lceil N^\gamma\rceil\le t\le N\) 满足
\[
 TW_{t-1}(r)\le W_t(r)\qquad(r\in\mathbb R).
 \tag{3}
\]

**证明。** 选定 \(K>2/\gamma-1\)，再增大固定 \(A\)。由第34轮的终端精确分解与统一 Chernoff 界，在左侧外部
\[
 O_t^--v_t^a\le t\delta/N+P_a\{S_t<r\}+B_t^\Gamma
 \le t\delta/N+Ct^{-K}+Ce^{-ct}.
 \tag{4}
\]
这里删去的参考 Gamma 修补项符号有利。右侧外部有
\(O_t^+\le Ct^{-K}\)。匹配中心与 \(t\mu\) 的差满足
\(t\delta=O(t\log N/N)=o(\sqrt{t\log t})\)，在 \(N^\gamma\le t\le N\) 一致成立，故第34轮右尾界可用。

增大 \(\Lambda\)，使 \(\Lambda b^2-\delta\ge b^2\)。由于
\[
 \frac{t^{-K}}{tb^2/N}
 \le C\frac{N^{2-\gamma(K+1)}}{\log N}\longrightarrow0,
 \tag{5}
\]
(4) 的两种尾项都由仿射余量吸收。因此左外部 \(I_t\ge O_t^-\)，右外部 \(I_t\ge O_t^+\)，两边都得到 \(W_t=O_t\)。在外部，直接以全局 \(O\) 比较；在内部，以 (1) 的局部比较及 \(W\le I,O\) 同时比较，即得 (3)。固定时刻差一的位移由严格多项式余量吸收。证毕。

这一证明没有对前一状态所在区域作限制，也没有使用“所有控制很少进入外带”的错误假设。与第35轮固定 \(\mu\)、无匹配仿射余量的反例不同，此处 \(b_N^2>0\) 是外部吸收的关键。

## 3. 初始化正缺口面积（严格）

对 \(1\ll T\le N\)，证明
\[
 \int_{\mathbb R}\bigl(\mathcal U_{T,N}-v_{T,N}^{a}-T\delta/N\bigr)_+dr
 \le C\log T+CT|a-\mu|+CT/N.
 \tag{6}
\]
常数对匹配动作的固定紧邻域一致。取 \(T=\lceil N^\gamma\rceil\)、\(0<\gamma<1/2\)，后两项为 \(o(1)\)，主项至多 \(C\gamma\log N+O(1)\)。

### 3.1 有界终端部分

令 \(\phi_\infty=\min(1,e^{-r})\)。第21轮全状态 CDF 屏障给每个可预测策略
\[
 F_{S_T^\pi}(q)\ge F_{G_T}(q)-2T^{-4},\quad
 G_T\sim N(T\mu+H_T,\mu(T+t_0)),\quad H_T,t_0=O(\log T).
 \tag{7}
\]
若 \(Y\sim\operatorname{Exp}(1)\) 独立，\(E\phi_\infty(r-S)=P(S+Y\ge r)\)，故真值上界为 \(P(G_T+Y\ge r)+2T^{-4}\)。不能将常数误差在无穷区间直接积分：取 \(\theta=1/2\)，统一条件 mgf 给
\[
 U_T^\infty(r)\le2M_*(1/2)^T e^{-r/2}.
\]
截于 \(R_T=2[T\log M_*(1/2)+5\log T]=O(T)\)，外尾积分为 \(O(T^{-5})\)，内部常数误差为 \(O(T^{-3})\)。负半轴两个有界终端均为 1。

一维 survival 面积距离等于 \(W_1\)，所以
\[
 \int(U_T^\infty-v_T^{\mu,\infty})_+
 \le W_1(G_T,S_T^\mu)+O(T^{-3})=O(\log T).
 \tag{8}
\]
最后一步：共同 Gaussian 耦合处理 \(H_T,t_0\)，代价 \(O(\log T)\)；iid 参考和的未标准化 Gaussian \(W_1\) 距离为 \(O(1)\)。这里可直接采用 Ross 的 Theorem 3.2（有限四阶矩版本）；本增量具有全部矩。不同固定动作之间用共同 Exp、共同 uniform Bernoulli 耦合，有 \(W_1(S_T^\mu,S_T^a)\le CT|a-\mu|\)，不需要动作间一阶随机序。

### 3.2 有限终端 call 部分

精确地
\[
 \phi_N(r-S)=\phi_\infty(r-S)+(S-r)_+/N.
\]
令 \(M_T=S_T-\sum_i m(c_i)\)、\(m=\mu T\)。统一条件指数矩给中心鞅 Bernstein 尾
\[
 P_\pi(|M_T|>x)\le C e^{-c\min(x^2/T,x)},\qquad S_T\le m+M_T.
\]
对 \(C_X(r)=E(X-r)_+\)，分别积分左右尾并用 \(EM_T=0\)，得到
\[
 \sup_\pi C_{S_T^\pi}(r)\le(m-r)_++H_T^*(r-m),\quad
 H_T^*(s)=C\int_{|s|}^\infty e^{-c\min(x^2/T,x)}dx,
 \quad \int H_T^*=O(T).
 \tag{9}
\]
参考 call 由 Jensen 大于确定均值 call，故其正面积误差为 \(O(T)\)。此外
\(C_{S_T^a}+T\delta\ge C_{S_T^a+T\delta}\)，后者与 \(S_T^\mu\) 同均值。任意同均值有限方差变量满足
\[
 \int|C_X-C_{X'}|\le\frac{\operatorname{Var}X+\operatorname{Var}X'}2,
 \tag{10}
\]
因为 \(C_X\ge C_{EX}\)，且 \(\int(C_X-C_{EX})=\operatorname{Var}X/2\)。两个参考方差均为 \(O(T)\)。合并 (8)–(10)，并对和的上确界用次可加性，即得 (6)。不需额外 \(W_2\) 定理。

### 3.3 尚缺什么

真值已知 \(\mathcal U_{T,N}\le O_T\)，而 \(I_T\ge v_T^a+T\delta/N\)。因此
\[
 E_T=(\mathcal U_{T,N}-W_T)_+\le
 (\mathcal U_{T,N}-v_T^a-T\delta/N)_+,
 \quad \int E_T\le C\gamma\log N+O(1).
\]
但 (3) 需要点态初值；面积界不能替代它。所需新引理应将这个特定误差传播到临界点，代价为 \(O(f_N\int E_T)\) 或同等有效形式。尚未证明此引理。不能把最优控制传输换成参考桥，也不能宣称所有控制密度都被 \(Cf_N\) 支配。固定短时的 \(4T\) 平移在 \(T=N^\gamma\) 时也不是可忽略代价。

### 3.4 任意面积传播的明确反例（严格）

令 \(T=\lceil N^\gamma\rceil\)、\(0<\gamma<1/2\)，\(y_0=Tm(a)\)，取固定 \(e_0>0\)、\(M=\log T\) 及
\[
 E_T(r)=e_0\mathbf1_{[y_0,y_0+M/e_0]}(r).
\]
它高度有界、面积为 M、支持在 T 自然管内。设 \(q=N\mu-d_N+O(\log N)\)、\(k=N-T\)。选择常值动作 \(c_N\) 满足
\[
 m(c_N)=(q-y_0)/k=\mu-[d_N-T\delta+O(\log N)]/k.
\]
该右端小于 \(\mu\) 且趋于 \(\mu\)；\(m\) 在 \(\mu\) 有严格二次最大值，故存在 \(c_N\to\mu\)，且
\(|c_N-\mu|\asymp(\log N/N)^{1/4}\)。于是 \(q-y_0\) 正是未来常值和的均值。该三角阵固定参数卷积有一致局部 CLT（参数在紧邻域，Gamma 平滑且具有一致指数矩），长度 O(log T) 的中心窗口密度一致至少 \(c/\sqrt N\)。因此
\[
 T^k E_T(q)\ge K_{c_N}^k E_T(q)\ge cM/\sqrt N,
 \qquad f_N^a(q)\asymp\log N/N.
\]
两者比值至少 \(c\sqrt N/\log N\to\infty\)。故不存在适用于任意这类非负误差的 \(T^kE_T(q)\le C f_N^a(q)\|E_T\|_1\)。一致局部 CLT 可按本项目第16轮 Fourier 方法：原点特征函数二阶展开在参数紧集一致，离原点模严格小于1，Exp 因子提供可积高次幂；中心 o(sqrt N) 窗口用 Fourier 反演即得。

这不反驳特定真值误差的补偿传播：在 \(T^k(W_T+E_T)-T^kW_T\) 中，同一偏移控制还牺牲基准 payoff。反例说明必须同时结算这项损失，不能仅用算子次可加性后单独估计 \(T^kE_T\)。

## 4. 双侧有界校正与 Mathematica 诊断（数值）

为减少左侧截断的非线性误差，测试
\[
 \mathcal B_L(v,U)=\frac{L}{1+(L-v)v^{-1}\exp[-UL/(v(L-v))]},\quad 0<v<L,\ U\ge0.
 \tag{11}
\]
它满足 \(v\le\mathcal B_L\le L\)、初值 \(v\)、初导数 1、初二阶导数 \((L-2v)/(v(L-v))\)。右侧趋于指数化校正，左侧趋于剩余缺口的指数衰减。Wolfram 精确符号核验这些恒等式及补对称关系。

在临界点，\(L-v\asymp\sqrt{\log N/N}\)、\(U=O(\log^2N/N)\)，Taylor 误差相对 \(f_N\) 为 \(O(\log^{5/2}N/\sqrt N)=o(1)\)。这仅说明临界价格合理，**未证明 (11) 是全局 supersolution**。

实际执行 `logistic_probe37.wl`：N=256，t=16,64,128，alpha=1.1，beta=0,.01，网格 h=.05,.025，双侧变换与截断指数变换共24组；另 N=10^8、同三时刻、beta=.01、h=.025共6组作短时诊断。beta 对应探索性 \(\sum\log^4j/j^{3/2}\)，不同于证明采用的固定大起点 log^6 误差项。无随机输入。

h=.025、beta=0 时，双侧变换最大正残差依次约 1.410e-5、1.484e-7、2.098e-8；指数变换约 2.151e-5、6.238e-7、1.496e-7。beta=.01 时，N=256 的正残差不超过 2.23e-16；N=10^8 不超过 1.78e-15。所有测试曲线单调。N=256 双侧核心最大残差依次约 −1.022e-6、−7.458e-8、−1.599e-8。

这些探针使用固定参考 mu、左侧屏障和探索性校正，未实现定理中的全部匹配参数及仿射余量；它们用于比较变换行为。上述为机器精度、分段线性 Bellman 核诊断，不是区间认证；极小概率可能下溢。库的6处50位独立公式核对最大绝对差约 1.11e-14。N=10^8 初次运行的 Exp[-N] 下溢已记录；新增精确分支 r≤N−1 时 Gamma 修补为0，避免无必要计算该指数。数学分支依据 \(\rho_N\in(N-1,N)\)，不改变结果。

## 5. 真实 mixed source/receiver 下界诊断

使用第30轮 seed：**n=64、b=32**，33组/65锚点；128是求积节点而非维数。每锚点独立 c_j，source likelihood 与 receiver 都重新生成，分母是实际 \(\Phi_{64}\) 成本，并按通过阈值的接收锥结算并集体积。

两模板 c_low=c0、c_high=c0+.25 的三个数值探针：

|模式|最小平均|接受组见证评分|
|---|---:|---:|
|uniform-low|0.999206|70.88094348|
|outer-high|0.997219|70.88167355|
|inner-high|0.999148|70.71817961|

评分为 log(V)−log(J)，其中 V 是接受锥并集体积，J 是实际 Orlicz 分母；进展缺口为 n log B 减该评分。同维固定 c0 标量参考 log 评分为71.71499209。所有锚点未同时通过不等于 pruned 接收集合无效，但三项都不竞争，也均未认证。固定两模板还受第28轮 ≤2Q_n 障碍，不能改变渐近 log n 系数。另测试33组变化曲线 c_j=c0+.25(v_j/16)^2 及反向曲线，评分分别为71.24294779、71.13579289，接受30、29组，最小平均分别为 .99751165、.99991139，均未超过同维标量参考。完整摘要见 mixedratio JSON；所有候选仍为浮点数值见证，未认证。

## 6. 复现、来源、下一轮

从项目根目录运行：

```sh
wolframscript -file code/extension_checks37.wl /tmp/extension37.json
wolframscript -file code/logistic_probe37.wl /tmp/logistic37.json
wolframscript -file code/logistic_probe37.wl /tmp/initial37.json initial
wolframscript -file code/lower_mixedratio37.wl /tmp/mixed37.json
```

附加50位求根检查在 N=4096、65536、1048576、16777216 验证反例中的动作均值方程，残差在50位精度内为零；动作偏移除以第四根尺度分别约3.001、2.904、2.879、2.883。这只是尺度复核，不承担一致局部CLT证明。

精确核验、机器扫描及下界数值摘要分别保存在本目录 JSON；代码共享第35–36轮必要卷积/核库。符号计算只是对应恒等式的复核，严格论证为上文。

来源索引：本项目 r021 全局 CDF 屏障、r028 固定锚定模板障碍、r034 外侧屏障、r035 匹配参数与 r036 局部非线性定理。外部使用 Nathan Ross, *Fundamentals of Stein’s method*, Probability Surveys 8 (2011), 210–293, DOI 10.1214/11-PS182；[原始预印本 Theorem 3.2，第10页](https://arxiv.org/pdf/1109.1880#page=10) 已核验，其有限四阶矩 Wasserstein 界足以给未标准化 iid 和 O(1) 误差。本轮不提出公开优先权声明。

下一轮最有区分力的两项动作：

1. 研究 (6) 的特殊正缺口在晚时段非线性比较中的传播，保留均值损失/Doob 权重；先证明或否定所需 L1 到临界点的比较，不以参考桥替代控制核。
2. 下界从固定模板数转为随组数增长的 c_j 耦合，继续联动 source、receiver 与精确并集体积；只有胜过同维标量参考的候选才进入高精度及定向认证。
