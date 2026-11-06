# 第43轮：同一势函数的全空间拼接与固定平移初始化

本文件以本轮proof.md两条早期左定理及已归档r036、r037、r040、r041为输入。结论是离散Bellman的全空间上界势与初始化；尚不在本文件中把它换算成反馈端点常数的最终锐渐近公式。

## 1. 一个势、一组常数

记alpha=1+epsilon>1，匹配a_N、b、delta、ell同proof.md。统一使用

s_t=B∑_(j=t₀+1)^t log^6j/j^(3/2)，s_(t₀)=0，0≤s_t≤1，
I_t=v_t exp[(alpha J_t+s_t f_t)/v_t]+Lambda tb²/N，O_t=min(O_t^−,O_t^+)。

源在t≤t₀置0。取Lambda足够大，使Lambda b²−delta≥c_*b²，且满足r036、r037所需下界；进一步增大Lambda只在一步差中增加有利的−Lambda b²/N。

按以下顺序选常数，避免密度初始化常数的循环：

1. 先将源的统一上界固定为1，选本轮early左带及r041正带的下管常数A_−、A_+。取gamma=1/3；取旧晚左floor引理的A_late。
2. 取本轮带外C_out，再取本轮C₀，使C₀ inf_a Lambda_a''(0)>2C_out。取自然管宽A_nat大于2 sup_a Lambda_a''(0) max(A_−,A_+)和2A_late/sqrt(gamma)，留严格重叠。
3. 用r036选源系数B及floor下界；用r040选足够小的右clip常数d，使其接口上界严格落在r041正带上界内。取Lambda满足全部下界。所有这些量都固定、不随N增长。
4. 最后取整数t₀足够大，满足各局部定理、相邻时间的缓冲，并令B∑_(j>t₀)log^6j/j^(3/2)≤1。增加t₀只减少该源，不违反统一上界。随后取N充分大使匹配偏移和固定clip严格余量成立。

令eta_mu(r/t)是Lambda_mu'(eta)=r/t的非负解；r/t≤mu时置0。定义全实轴候选

W_t(r)=min(O_t(r),I_t(r))，eta_mu(r/t)≤2d epsilon；
W_t(r)=O_t(r)，eta_mu(r/t)>2d epsilon。                       (1)

r040给接口带d epsilon≤eta_mu≤2d epsilon上I_t≥O_t^+，故W_t=O_t。**远强制区不要求O≤I**，那里正是允许旧内分支再进入而主动弃用它。

## 2. 所有空间/时间区域的一步核查

始终W_(t−1)≤O_(t−1)，故K_cW_(t−1)≤O_t对所有c成立。只需核查实际选择内分支的点。

### 2.1 右侧及自然管

当前eta_mu≥d epsilon时，由接口带或强制定义W_t=O_t，比较立即成立。当前eta_mu<d epsilon时，全部可达前态r−z≤r；时间从t变到t−1只令该固定紧区间的正倾斜增加O(1/t)。取t₀大，所有前态eta_mu<2d epsilon，因此

W_(t−1)(r−z)≤I_(t−1)(r−z)，z≥0。                            (2)

不把(2)错误地延到当前已经在强制区的点。

自然管|r−tm(a)|≤A_nat sqrt(t log t)内用r036的全动作I比较，配合(2)和粗屏障得K_cW_(t−1)≤min(O_t,I_t)=W_t。正侧管外但eta_mu<d epsilon时，eta_a的下界落入r041正中偏差带；其上界因固定clip的严格缩窄仍在r041范围。匹配偏移|a−mu|+b=O(sqrt(log N/N))≤C sqrt(log t/t)（N≥t）保证其假设，故同理成立。

### 2.2 早期左侧：t≤sqrt(N)

带外r≤tm(a)−C_out t^(2/3)由本轮定理2有I_t≥O_t^−≥O_t，当前W_t=O_t。

该带外与自然管之间，由严格常数重叠，当前点的负倾斜在本轮定理1范围内。此时r<rho_N，所有前态也小于rho_N，且不可能到达正侧强制区。于是

W_(t−1)≤min(O_(t−1)^−,I_(t−1))
         =min(L_(t−1)^mu,I_(t−1))。

本轮定理1给其核平均≤min(L_t^mu,I_t)。再同时用K_cW_(t−1)≤O_t^+，得到所需三分支最小值。自然管的固定缓冲处理两端重叠；没有在接触环带留下未经验证的一步。

### 2.3 晚期左侧：t≥N^gamma

这里只复用r037/r038已证明的晚floor机制，不作为新增结果。对r≤tm(a)−A_late sqrt(t log N)，精确仿射恒等式给

O_t^−−I_t≤P_a(S_t<r)+B_t(r)−(Lambda b²−delta)t/N。

统一负Chernoff界令概率≤N^(−M)，M可由A_late增大到任意所需常数；B_t≤Ce^(−ct)。因为t≥N^gamma、b²≈log N/N，二者均小于c_*tb²/N，故I_t≥O_t^−、W_t=O_t。该点无需r<rho_N：A_t≤P(S_t<r)在整个左侧仍是上界。

剩下的左侧点满足tm(a)−r≤A_late sqrt(t log N)≤(A_late/sqrt(gamma))sqrt(t log t)，已在r036自然管内。晚区与早区N^gamma≤t≤sqrt(N)严格重叠；每层使用同一个I、源、floor和O。正侧仍由2.1处理。

所有t₀<t≤N、全部实r和c≥0至此均覆盖。因此

K_cW_(t−1)(r)≤W_t(r)，T W_(t−1)≤W_t。                     (3)

## 3. 整个屏障平移后的固定短时初始化

不要求未平移的零密度源I_(t₀)直接支配真值。首先O_t是已有全局粗上界，所以O_t≥mathcalU_t≥v_t；I_t≥v_t，强制区也用O。故

W_t(r)≥v_t^a(r)，所有r。

对任意可预测策略，单步密度≤2e^(−z)给

P_pi(S_k>x)≤2^k P(Gamma(k,1)>x)。

Gamma整数阶hazard递增，且P(Gamma(k,1)>4k)≤(4e^(−3))^k<2^(−k)。故r033的全策略随机支配加强为对任意固定a≥0均有

S_k^pi≤st Gamma(k,1)+4k≤st S_k^a+4k。

phi_N非增，所以取固定h=4t₀有

mathcalU_(t₀,N)(r)≤v_(t₀,N)^a(r−h)≤W_(t₀,N)(r−h)。          (4)

定义**整个**平移屏障widehat W_t(r)=W_t(r−h)。平移与每个K_c交换，(3)原样成立于widehat W；由(4)迭代得

mathcalU_(N,N)(r)≤W_(N,N)(r−h)，h=4t₀固定。                 (5)

不能把(4)偷换成v(r−h)≤W(r)。在r≤0处平移产生的h/N由平移屏障自身承担，不要求b²/N量级的floor去支付它。

## 4. 临界点固定平移的代价

设q_N=N ell_N=N mu−d_N+O(log N)，r=q_N+O(log N)。这里eta_mu=0，故(5)右侧≤I_N(r−h)。临界参考密度f_N(q_N)≈log N/N是既有r022的结论。负倾斜lambda_N≈sqrt(log N/N)，Nlambda_N→∞，有限终端Mills/辅助密度导数给v_N≥c>0、|v_N'|≤C f_N(q_N)。精确Poisson公式在r−Nell=O(log N)上给

J_N/f_N=O(log N)，J_N'/f_N=o(1)，f_N'/f_N=O(lambda_N)。

所以U_N'=o(f_N)、U_N/v_N=o(1)，由

V_N'=exp(U_N/v_N)[(1−U_N/v_N)v_N'+U_N']

得到|I_N'|≤C f_N(q_N)，floor空间导数为0。于是固定h只导致

|I_N(r−h)−I_N(r)|≤C_h f_N(q_N)。                            (6)

它不产生新的log N级阈值位移。后续仍须把(5)–(6)代入实际反馈端点与任意输入归约，逐项计算锐系数；本文件不以势函数拼接替代这一步结算。
