# 候选：当前bin宽下的Lipschitz方差帽

本文件在先运行 `variance_probe.py` 后给出数学推导。候选只涉及source class核在连续uniform query-bin内的方差；不改本轮冻结source、driver、payload、CI或主结论。

## 先运行的小验证

`CANDIDATE_VARIANCE_RESULTS.json` 保存32个Fraction真实box case及24个common-D真实核数值case，总.0701秒。

Fraction小源为K4/R4、w4/v8、D0、真实box宽δ=1/128，query-bin全宽h=δ/6；分别用p=.5的正背景源和p1纯峰源，表0/3，plus/minus class，预定进入/退出ramp、class边界和N0控制。每个核按所有breakpoints分片，精确积分线性函数及其平方，得到真实连续uniform-bin variance。32例都满足新帽；其中10例非零affine ramp达到新帽等号。这些不是38轮两点隐藏law，后者不能使用1/12常数。

common-D数值核使用当前冻结Source的R4/R32、表0/末表、L7/6，并在每class的N_partial=0/1/2各取预定first bin。按breakpoints分片24点Gauss计算variance，24例均低于新帽，最大实际variance/新帽=.9938832。N0例仍有非零background方差，证实不能把N0当方差0。普通FP及Gauss数值不是完整interval证书；D0有理部分与下面一般推导分开。

## 1. 双独立uniform差平方恒等式

设X、X'是同一个全宽h区间上的独立uniform变量，Z=z(X)。若z在这个bin上具有Lipschitz常数L_loc，则

`Var Z = (1/2) E[(z(X)−z(X'))²] ≤ (L_loc²/2) E[(X−X')²] = L_loc² h²/12`。

等号可由affine z达到。这对一维连续uniform hidden query成立；对两点端点law，差平方期望不同，variance可达L_loc²h²/4，不能沿用1/12。

## 2. 当前range公式确实等于L_loc h

固定表、class、bin和macro L。source class核是窗口积分除L，绝对连续；其a.e.导数是两端class density值之差。

1. 常数背景的class overlap导数绝对值≤1，贡献≤(1−p)/L。
2. common D-mask贡献导数绝对值≤(p/K)2Dsup/L，即使端点经过class边界，窗口积分仍连续/绝对连续，该a.e.帽给Lipschitz帽。
3. 每个peak的endpoint slope绝对值≤1/δ。bin内曾产生非零slope的atom必在某window endpoint的δ+h候选区内；source的N_partial含所有这样的候选。class内部gap证书保证每端≤1，总N≤2。因此peak贡献≤(p/K)N/(δL)。

故整个bin上可取

`L_loc=[(1−p)+(p/K)N/δ+(p/K)2Dsup]/L`。

当前h=1/(6wKv)=δ/6，所以冻结源码的

`range_cap=[(1−p)h+(p/K)N min(1,h/δ)+(p/K)2Dsup h]/L`

恰好就是L_loc h。于是source中原先用于relative variance的 `range_cap²/4` 可以解析地改为 `range_cap²/12`，改进因子1/3。

此推导需要h≤δ。若h>δ，原range中的min会饱和，range_cap/h不再是同一个Lipschitz帽；不能一般把任意振幅平方的1/4改成1/12。N_partial也必须覆盖整bin所有slope区域，而不是只看中心是否在峰内。当前endpoint候选算法满足此前提；N0仍保留全部background项。

## 3. 从单核到产品与分层variance

均值下包保持不变，因此每个坐标class relative variance帽可由旧r_i替换为r_i/3。表产品帽新值为

`R_new=∏_i(1+r_i/3)−1 ≤ (1/3)[∏_i(1+r_i)−1]=(1/3)R_old`。

证明只需展开正系数多项式：k阶项由3^{-k}≤1/3控制。重复bin counts也只是重复的n个坐标。**不能反向用** `(1+R_old)^{1/3}−1` 作上帽；该式一般是下方向的比较。

相同的表均值上包u_b下，分层Minkowski帽满足

`V_new ≤ (Σ_b u_b sqrt(R_old,b/3))² = V_old/3`。

没有新增R或monomial数费用。fixed0、joint、prior均可用此共同因子，其他PB均值/BE/π上下包不变。

## 4. 保存conditional lower的确定变换

设旧某L的Cantelli lower为b=g²/(V_old+g²)，g=(mean_lower−λ)_+；g0时b0。用V_old/3仍为合法variance帽，因此新lower至少

`φ_{1/3}(b)=3b/(1+2b)`。

一般c∈(0,1]时 `φ_c(b)=b/[c+(1−c)b]`。φ单调，故对共同尺度lower max有

`max_L φ_c(b_L)=φ_c(max_L b_L)`。

所以保存的conditional lower B(N,h)足以给一个更强的确定lower，不需重新求核或重采样。当前payload保存的是贡献x=λW B与logD；在理想实数算术下可恢复B=x exp(logD−h)，再变换。未来实现应使用稳定log域/方向正确的包围，不能假装underflow零能恢复丢失信息；保守保留0仍可作较弱lower。

由于新lower仍受真实条件概率/均值支配，旧source-only Markov贡献帽依然可用。此φ变换不直接给保存的union/Markov upper同样改进；upper要保留每L分解或另证明新映射，不能盲改。

本轮只给候选和小验证，**未应用φ到任何pilot payload，未生成新CI**。下一轮若采用，应冻结后处理与统计报告口径，再用现有payload确定重算；原504项CI继续原样保留。普通FP证据边界与原均值/BE误差前提仍须携带。
