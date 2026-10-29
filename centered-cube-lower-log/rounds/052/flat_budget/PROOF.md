# 方差平坦尺度族上的幅度混合比较

先行计算为本目录probe.py的12项实际参数有理预算；随后check_order.py又实际完成64项共同排序有理检验，并给出无共同排序时弱泛函由1增至3/2的有限模型反例。表内有限n的D=.03/1明确是待满足的方差平坦假设，未冒充有限n核认证。以下结果为自含条件定理；随后说明第45轮接受源的渐近参数满足它。

## 1. 设定与有限常数

在概率空间(T,dx)上，对有限共同尺度族Lambda（基数J≥1），给定实函数d_L，满足∫d_L=0、|d_L|≤h。令
q_L=n∫d_L²，且|q_L−v|≤D，q_L≤qbar。
各坐标取乘积概率测度。对|alpha|≤A、Ah<1，定义
T_{alpha,L}(x)=prod_{i=1}^n(1+alpha d_L(x_i))。
它是非负质量1密度。真实周期立方体应用中d_L=A_L rho−1；只在源1+alpha(rho−1)合法时把T称为实际产品来源的输出。

令nu为[-A,A]上的任意固定概率测度，允许依赖n；所有尺度与查询共用它。定义M_L=∫T_{alpha,L}dnu。假设选择epsilon>A h qbar，置
 t=epsilon−A h qbar>0，
 E=t²/[2((1+Ah)h²qbar+h²t/3)]，
 delta=A²(D+epsilon)/2+n(Ah)³/[3(1−Ah)]。
记W为基准概率空间上的弱(1,1)泛函，使用严格越阈集合。

若nu仅支撑非负alpha，或仅支撑非正alpha，则
 W(max_L M_L) ≤ exp(2delta)∫W(max_L T_{alpha,L})dnu +2J² exp(−E)。 (1)
对一般含两种符号的nu，右端积分前系数改为2exp(2delta)，最后的2J²exp(−E)保持不变。特别可将积分界为其支撑中单一产品弱泛函的上确界；不假设该上确界有O(sqrt n)上界。

## 2. 以输出质量控制坏集

令V_H=sum_i d_H(x_i)²，S_L=sum_i d_L(x_i)，坏集
 B=union_{H∈Lambda}{|V_H−q_H|>epsilon}。
对任一alpha,L，用密度T_{alpha,L}定义乘积概率Q_{alpha,L}。在该律下
 |E_Q V_H−q_H|=|alpha n∫d_L d_H²|≤Ah q_H≤Ah qbar，
 sum_i Var_Q(d_H(x_i)²)≤n E_Q d_H⁴≤(1+Ah)h² qbar。
每个独立随机变量d_H²位于[0,h²]。由Bernstein界，
 Q_{alpha,L}(|V_H−q_H|>epsilon)≤2exp(−E)。
为避免外部尾界假设，所用Bernstein形式可直接由中心变量Z_i、|Z_i|≤b=h²的幂级数推导：E|Z_i|^k≤b^{k−2}E Z_i²，k!≥2·3^{k−2}（k≥2），故0<s<3/b时log E exp(s sum Z_i)≤s² sigma²/[2(1−sb/3)]。取s=t/(sigma²+bt/3)，再对负尾使用−Z_i，得到上式。

对H作union bound，再由Tonelli及max≤sum得到
 ∫_B max_L M_L dx ≤ sum_L∫Q_{alpha,L}(B)dnu≤2J² exp(−E)。 (2)
因此即便所选lambda很大，坏集的弱型费用仍由(2)支付，而不是把极小Haar概率直接当作输出质量。未调用待证极大弱型界。

## 3. 好集上的共同排序

由log(1+u)的幂级数，当|u|≤Ah<1时，
 |log(1+u)−u+u²/2|≤|u|³/[3(1−Ah)]。
所以在B的补集上，同时对所有alpha,L，
 |log T_{alpha,L}−[alpha S_L−alpha²v/2]|≤delta。 (3)

若alpha≥0，则理想模型exp(alpha S_L−alpha²v/2)在所有alpha下都由同一S*=max_L S_L取最大。令
 F(s)=∫exp(alpha s−alpha²v/2)dnu(alpha)。
F连续非减。由(3)，在好集上max_L M_L≤exp(delta)F(S*)，且对每个alpha有exp(alpha S*−alpha²v/2)≤exp(delta)max_L T_{alpha,L}。

需要的共同排序引理是：对有限测度mu、有限有界随机量Z、非负连续非减函数phi_a，F=∫phi_a dnu，并假设F在Z的有限闭值域上有限连续，
 W_mu(F(Z))≤∫W_mu(phi_a(Z))dnu。 (4)
在当前紧参数及有界Z下F连续、积分有界。证明：给lambda及正测度越阈事件H={F(Z)>lambda}，取t=ess inf_H Z；由连续性lambda≤F(t)。又phi_a(Z)≥phi_a(t)在H上成立，令阈值从下趋近phi_a(t)，得W_mu(phi_a(Z))≥phi_a(t)mu(H)。积分即得lambda mu(H)≤∫W_mu(phi_a(Z))dnu；再取sup。零质量事件无影响，phi_a恒常（alpha=0）亦包含在内。

以mu为基准测度限制到好集，在(4)中用phi_alpha(s)=exp(alpha s−alpha²v/2)，再两次使用(3)，得到
 W_good(max_L M_L)≤exp(2delta)∫W(max_L T_{alpha,L})dnu。
好集和坏集在同一个lambda下分开，(2)支付后者，得到(1)。负alpha时改用−min_L S_L，完全相同。

对混合符号，先把nu按alpha≥0和alpha<0分开。好集上max(M_++M_−)≤max M_++max M_−；用两个阈值lambda/2的事件union bound，会真实支付因子2。对两份有限子测度分别使用上述共同排序引理后，积分权重相加为1。坏集仍直接对完整M使用(2)，所以不再多付一次坏集预算。这是有明确阈值费用的比较，未声称弱L1函数相加具有常数1三角不等式。

## 4. 接受源渐近参数的回代

第45轮接受源在安全局部尺度上有|d_L|≤C/K，q_L=v+O_w(1)，n=Theta_p(K²v)，且K/v^r→∞对每个固定r成立。取任何固定A、D足够大，qbar=nh²=O(v)，epsilon=1/v。则
 Ah qbar=O(v/K)=o(epsilon)，
 E≥c K²/v³，
 n(Ah)³=O(v/K)=o(1)。
若J≤C'K^u v^r、u,r固定，则2J²e^(−E)=o(1)，delta=A²D/2+o(1)。于是固定符号幅度混合满足
 W(max_L M_L)≤[exp(A²D)+o(1)] sup_alpha W(max_L T_{alpha,L})+o(1)， (5)
一般符号将系数改为2exp(A²D)+o(1)。常数不依赖尺度数、prior原子数或其最小权重。不从这里推导单产品的O(sqrt n)上界；(5)仅说明此幅度混合相对同尺度最佳产品没有无界的对数放大。

若原网格含不安全节点，只能另行支付这些节点的质量1 Markov费用；使用第50轮单bank的坏节点计数，在间距≥a/(Kv)的I内有限网格上，节点数为poly(v)=o(sqrt n)。因此完整该网格的比较多一个o(sqrt n)，不是无条件o(1)。

第51轮全部计数的一致比较把任何非负计数门权归约为这样的合法幅度产品混合；合法soft split保证|tau|统一有界，故alpha区间统一紧。结合(5)，在这些指定网格内，任意计数门权的输出弱泛函≤常数倍最佳合法倾斜产品+o(sqrt n)。第51轮的乘性o(1)保留为乘性误差，未在缺少速率时把它乘J后称为o(sqrt n)。

周期模型向Lebesgue可积截断输入的转移仍只使用第50轮同一有限网格桥梁；J为上述多项式时可取足够大的整数截断，使边界费用o(1)。本证明没有估计完整连续所有尺度的C_n；没有得到新的下界，也没有证明原单产品构造已经最优。
