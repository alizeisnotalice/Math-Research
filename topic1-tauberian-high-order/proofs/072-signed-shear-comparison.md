# 正负共同剪切的真实九阶比较：保留全源交叉退款

2026-11-05。冻结送审候选，待独立审查。本文讨论指定长见证及共同孔质量的真实有限输入比较，不识别C_n的绝对九阶系数，不加入微孔。n=3是首个非平凡情形；以下证明固定每个n≥3。

## 1. 参数、算子与待证比较

记u=t^(1/n)、v=2u、a=(2u)^n、g=(1+2u)⁻¹。κ(u)在固定正紧区间，Λ(u)有界；不要求解析。置
\[
 e=\kappa u^3,\quad
 T_+=(1-e)I+eJ,\quad T_-=(1+e)I-eJ,\quad
 w=ug+\Lambda u^6,\quad \rho=1-2w,\quad C=[-w,w]^n,
 \quad M=a\rho^n.
\]
两符号始终使用同一w、ρ、M及每个源点的原规则半径
\[
 s=w-\|x\|_\infty,\qquad r(x)=\rho+s.
\]
令b=(2w)^(n−2)，
\[
 F(s)=\frac a b[(\rho+s)^n-\rho^n],\quad
 f=F'(0)=na\rho^{n-1}/b,\quad
 h=e^2/f,\quad\delta=e/f,\quad\lambda=e/w,\quad
 V_0=b e^3/f.
                                                               \tag{1}
\]
固定参数范围后，小u时
\[
 f\asymp u^2,\quad \delta\asymp u,\quad \lambda\asymp u^2,
 \quad h=e\delta\asymp u^4,\quad V_0\asymp a u^5.
\]
特别V₀δλ∼a u⁸，**不能丢**；而V₀δ²λ、V₀λ²、V₀h/ρ均O(a u⁹)。

对任意可测删除集V⊂C，定义完整实际坏集
\[
 R_\pm(V;x)=\frac1b\sum_{k\in\mathbb Z^n\setminus\{0\}}
 |(T_\pm k+(C\setminus V))\cap Q(x,\rho+s)|,\qquad
 \mathcal B_\pm(V)=\partial C\cup\{R_\pm(V;x)\ge F(s)\}.
                                                               \tag{2}
\]
U₀^±=B_±(∅)，U_(j+1)^±=B_±(U_j^±)。设计M不按测试V重算。

核心候选为：对每个固定j≥0，
\[
 \big||U_j^-|-|U_j^+|\big|\le C_{n,j} a u^9.             \tag{3}
\]
故若某个固定偶层J≥2的正剪切参数有真实安全余量
\[
 D-|U_J^+|\ge c\,a u^8,\qquad D=(2w)^n-M,\quad c>0,     \tag{4}
\]
负剪切同一参数、同一质量的U_J^-也安全。两有限输入可取得
\[
 \mathscr R_-(u)-\mathscr R_+(u)
 =\left[\frac{4\binom n3\kappa_0^3}{2^9}+o(1)\right]av^9,
 \qquad\kappa(u)\longrightarrow\kappa_0.
                                                               \tag{5}
\]
若κ₀=8n/3，系数为4 binom(n,3)n³/27，n=3时为4。

下面证明(3)时完整保留V₀δλ的cross-pair退款，并证明其正负相等。不能仅用det展开或沿用旧balanced-even清单。

## 2. 全整数邻孔清单及真实浅层

若原Q(x,r(x))与T_±k+C有正交量，必有
||T_±k||∞≤||x||∞+r+w=1。另一方面
||T_±k||∞≥[1−(n−1)e]||k||∞。取e小，使2[1−(n−1)e]>1，所有||k||∞≥2的整数均被排除，不是有限枚举的外推。

对k∈{−1,0,1}^n记S=Σk_i。
正剪切T_+k=(1−e)k+eS1。S≠0时，选与S同号的非零坐标，中心模≥1，只有零体积接触。正交量只来自S=0，支撑大小j≥2为偶数，支撑中心为±(1−e)，其余中心0。

负剪切T_-k=(1+e)k−eS1。若k含正负非零坐标，要使正坐标模<1须S>1，使负坐标模<1须S<−1，不可能。故只有同号非零坐标、支撑j≥2：
\[
 k=\pm1_A,\quad |A|=j,\qquad
 (T_-k)_i=\begin{cases}
 \pm[1-(j-1)e],&i\in A,\\
 \mp je,&i\notin A.
 \end{cases}                                          \tag{6}
\]
支撑1只有边界接触。负剪切新增所有奇支撑，首个是j=3；n=3有6个j2邻孔及2个j3邻孔。

上述清单也是“某原半径可有正cap”的完整清单：中心在开单位盒时x=0、r=1−w给每轴正交长。但深部可能同时捕获多个目标，不能先将R写成一个浅角帽。

每目标支撑交长≤C_n e，非支撑全宽进入，因为2w+ne<ρ。归一化单目标质量≤C_n e^j/w^(j−2)。有限清单给R_±≤C_n e²，所以所有内点坏点s≤C_n h≪w−ne。浅层内，目标支撑方向必须与源最近面符号一致，否则该方向穿透≤C_n e−w+s<0。

在最近面符号σ及全局最浅轴q的图中，写d_l=s+y_l，y_q=0，其他y_l∈[0,w]。源纤维上端
H_top(y)=w−max_l y_l，雅可比1。一个支撑j的活动目标穿透为
\[
 ((c_j e-y_l)_+)_{l\in A},\qquad
 c_j=1\text{（正平衡）},\quad c_j=j-1\text{（负同号）}, \tag{7}
\]
不依s；非支撑捕获目标完整原宽。因此真实坏高度是
min{H_top,F⁻¹(全部目标帽质量之和)}。所有目标来自不同周期孔、彼此不交，故质量先在同一个原见证内相加。本文保留全部图及真实clip。

## 3. 去掉高支撑及容量曲率的已付误差

固定一个支撑j、一个目标点。Tonelli积分其可见源底面：
q∈支撑时底面积≤C_n e^(j−1)w^(n−j)，q不在支撑时≤C_n e^j w^(n−j−1)。除归一化b与f，得到该目标对坏集L¹差的小常数至多
\[
 C_n\delta\lambda^{j-2}+C_n\delta\lambda^{j-1}.
\]
有限目标求和为C_nδ<1/2，对正／负、任意测试V,W都成立；真实F′≥f和clip的1-Lipschitz性没有放大该常数。由此U_j体积≤C_nV₀，所有输出高度≤C_nh。

完整名义帽同样积分：支撑j的坏体积上界为
C_n V₀ λ^(2j−4)（q不在支撑时再多λ）。所以两符号全部j≥3目标的合计≤C_n V₀λ²=O(a u⁹)。负j3条不能宣称不存在，只是其实际预算在这里已有这个上界。

又F(s)=f s[1+O_n(s/ρ)]，0≤s≤C_nh。因此对任意测试V，真实F⁻¹与线性s=R/f的输出体积差≤C_n V₀h/ρ=O(a u⁹)。

定义辅助算子A_±：只保留j2目标，容量用f s，但**保留全H_top clip**。将其迭代记V_j^±。从相同空集出发，以上逐步误差及C_nδ小常数给每个固定j
\[
 |U_j^\pm\triangle V_j^\pm|
 \le C_{n,j}V_0(\lambda^2+h/\rho)=O_{n,j}(a u^9).       \tag{8}
\]
辅助算子仅为误差比较工具，不把它自身当真实见证条件；最终构造使用完整U_j^-。

## 4. 两轴算子与unclipped体积核

全文E、F一律是**无序坐标轴对**，不是有向角。轴对E={p,r}的pair算子A_±^E保留其两个相反方向j2目标、全部符号图及全局最浅轴图。同轴对两反角之间的首退款可达a u⁶，已完整收入isolated-pair账本；不能将它另列作a u⁸的cross。只有E≠F才是后文的两边forest。
正剪切最近面符号满足σ_pσ_r=−1；负剪切满足σ_pσ_r=+1。两轴临界穿透均为(e−y_l)_+。负pair非临界中心虽为∓2e，但全部非临界原宽都在见证内，故不改其local target cap算子。

记孤立pair迭代V_j^{E,±}。单pair的两个符号算子经反射一个轴精确共轭，含全局最浅轴不在pair的图、原ρ、f和H_top。因此
\[
 |V_j^{E,+}|=|V_j^{E,-}|\quad\text{对全部固定j}.          \tag{9}
\]

为写共同体积账本，先将各图R/f高度作不裁顶的代数延伸，再扣回clip。此延伸不是对C外源点声称有见证。固定一条有向target pair，令d₁,d₂为目标点到相应两面对端的物理深度。全图Tonelli给unclipped退款核
\[
 K_E(z)=\boldsymbol1_{\{d_1\le e\}}(e-d_2)_+
       +\boldsymbol1_{\{d_2\le e\}}(e-d_1)_+
       +\frac{n-2}{w}(e-d_1)_+(e-d_2)_+.               \tag{10}
\]
每个无向pair把其两个面对角的核相加。q在pair内给前两项；q在pair外时，两支撑底变量长度之积，其他最近面符号求和并除b，每个外轴贡献1/w，故末项n−2。特别n=3的外最浅轴项不能省略。

写
Q_E(V)=f⁻¹∫_C K_E1_V。K_E≤C_ne，故|Q_E(V)−Q_E(W)|≤C_nδ|V△W|。在各固定符号／最浅图内，K_E沿s在跳平面外的导数≤C_n，跳幅≤C_ne；跳仅在某d_i=e。

设N_E为单pair完整名义unclipped体积，C_E(V)为其源顶端裁剪量。joint名义N=ΣN_E，且有精确式
\[
 |A_\pm(V)|=N-\sum_EQ_E(V)-\mathsf C(V),\qquad
 |A_\pm^E(V)|=N_E-Q_E(V)-\mathsf C_E(V).                \tag{11}
\]
全源joint高度依目标总和，不是假设各pair的坏点集合互不交。

## 5. 实际顶端clip的费用

所有unclipped高度≤C_nh。clip发生必有某自由底坐标y_l≥w−C_nh。任何涉及该坐标的pair帽均为零；其余pair帽及删除退款与此自由坐标无关。沿该轴积分，裁剪带宽≤C_nh，相对完整自由宽w为O(h/w)。所以对|V|≤C_nV₀，
\[
 0\le\mathsf C(0)-\mathsf C(V)
 \le C_n(h/w)\sum_EQ_E(V)
 \le C_nV_0\delta(h/w)=C_nV_0\delta^2\lambda.           \tag{12}
\]
孤立pair同理。这保留了名义clip本身可能为V₀h/w=V₀δλ∼a u⁸的九阶费用，没有把它扔进误差。

joint名义clip与孤立名义clip之和的差，只有至少两个pair同时活跃时发生。它们占至少三条不同坐标面，至少两个底变量在O(e)带；另有一个不参与的自由坐标在O(h)中分面带。故差≤C_n h²e²w^(n−4)=C_nV₀δλ²。n=3没有这种额外自由轴，此差为零。若最浅轴不在活跃坐标并中，需更多O(e)底变量，只更小。

因此对各实际迭代j≥1，
\[
 |\mathsf C(V_{j-1})-\sum_E\mathsf C_E(V_{j-1}^E)|
 \le C_nV_0(\delta^2\lambda+\delta\lambda^2)
 =O(a u^9).                                           \tag{13}
\]
这些是同一次source clip账本，不额外消费任何删除质量。

## 6. joint层与孤立层的比较：先保leading cross价

以下略去±。在共同最近面图中，记joint高度H_j、孤立高度H_j^E。令S_j为高度min(H_top,Σ_EH_j^E)的真实可测层；它包含每个孤立层。单pair高度≤h，仅在对应两轴最近面深≤e+C_nh处可非零。

设D_j=|V_j△S_j|。对j=0，
min(t,Σx_E)=min(t,Σmin(t,x_E))，故D₀=0。min的1-Lipschitz及(10)给
\[
 D_j\le C_n\delta D_{j-1}
       +\sum_EQ_E(S_{j-1}\setminus V_{j-1}^E).
\]
后一集合进入K_E支撑后，另一个pair也活跃，所以至少三坐标近面；高度O(h)，体积≤C_n h e²w^(n−3)=C_nV₀λ。故
\[
 D_j\le C_n\delta D_{j-1}+C_nV_0\delta\lambda,\qquad
 D_j\le C_{n,j}V_0\delta\lambda=O(a u^8).               \tag{14}
\]
这个量可能真在九阶，不能丢；但将它再次放进体积退款Q时另乘δ，才是O(a u⁹)。

还需weighted加性，而非将S_j当各孤立层的union。对每个E，
\[
 Q_E(S_j)=\sum_FQ_E(V_j^F)+O_n(V_0\delta^2\lambda+
                                      V_0\delta\lambda^2). \tag{15}
\]
证明：沿每个物理s纤维，ΣH_F以堆叠方式代入积分。常数权精确相加。远离K_E跳平面，逐层堆叠误差≤(C_n/f)Σ_(F<G)H_FH_G；只有至少两个pair活跃的三面底区贡献，其底面积≤C_ne²w^(n−3)，给C_nh²e²w^(n−3)/f=C_nV₀δ²λ。若一条d_i=e跳平面被纤维穿过，相关底变量在宽O(h)的带，另至少一个变量在O(e)带；跳幅≤C_ne/f，乘高度O(h)及底带O(h e w^(n−3))仍给同一量。若目标pair面与最近面符号相反，K_E在该浅层为零，不额外产生跳。ΣH的顶端截断只在多pair加自由中分面带，§5的同一几何给更小的V₀δλ²上界。这证明(15)。

孤立pair还满足各固定j
||H_j^F−H₀^F||∞≤C_nhδ。因为其U₀在任意pair cap中的质量≤C_n h e b：pair内最浅轴给薄臂，pair外最浅轴给额外h薄层，后者再小e/w；除b及f即得。与不同E的cap相交至少三面，因此
\[
 Q_E(V_j^F)-Q_E(V_0^F)
       =O_n(V_0\delta^2\lambda),\qquad E\ne F.           \tag{16}
\]
这里V₀^F是完整n维孤立pair首层；包括最浅轴在pair外的条，绝不是把它换成独立二维plane。

将(11)减去各孤立pair账本，用(13)、(14)再进Q、(15)、(16)，对每个j≥1得到
\[
 |V_j|=\sum_E|V_j^E|-
          \sum_{E\ne F}Q_E(V_0^F)
       +O_{n,j}\!\left(V_0\delta^2\lambda+V_0\delta\lambda^2\right).
                                                               \tag{17}
\]
**cross和本身为V₀δλ∼a u⁸，完整保留。** 本公式不是把joint U₀当union，也未删pair外最浅轴。j=0另由名义线性加性与§5给Σ孤立体积加O(V₀δλ²)。

## 7. 两边forest共轭：保留的交叉项正负精确相同

设E、F为两个不同轴对。它们在坐标轴顶点构成至多两条边，必为forest。存在对角反射D=diag(ξ_l)，ξ_l∈{±1}，同时使每条边两端ξ乘积为−1。共享轴例如E=12、F=13，可取(ξ₁,ξ₂,ξ₃)=(1,−1,−1)；不相交边分别选一端反射。

该反射并不把完整T_+晶格变为T_-。实际使用的是§4已经证明的局部pair帽：非临界方向精确吞入全目标宽，∓2e的中心错位完全退出帽算子；剩下两临界坐标的穿透相同，仅最近面符号积由−1变+1。反射保持C、s、f、顶端clip、体积和全部全局最浅轴图。因此可对这个给定(E,F)同时得到
\[
 K_E^+(Dz)=K_E^-(z),\qquad
 D(V_0^{F,-})=V_0^{F,+}
\]
（同一有向面对角可相应重编号）。所以
\[
 Q_E^+(V_0^{F,+})=Q_E^-(V_0^{F,-}),\qquad E\ne F.       \tag{18}
\]
包括F外最浅轴、K的(n−2)/w项、所有面符号。可对每个两边项分别选反射，不要求一个反射同时满足整张complete graph。

联立(9)、(17)、(18)，得对每固定j≥1
\[
 \big||V_j^-|-|V_j^+|\big|
 \le C_{n,j}V_0(\delta^2\lambda+\delta\lambda^2)=O(a u^9).
\]
j=0用名义线性和clip同样成立。再用(8)恢复完整F及所有高支撑，得到(3)。

三条边的triangle符号积不被坐标反射改变，所以本证明没有把所有joint高阶相互作用都宣称正负同价。它仅证明影响坏体积到a u⁸的项，要么是一个完整孤立pair（精确共轭），要么是保留的两个pair交叉退款（forest精确共轭）。剩余三面／曲率／odd-support误差以a u⁹付款，落到最终halo至少十阶。

## 8. 相同设计质量的两个真实有限输入

取共同κ,Λ，并假设正偶层J有(4)。由(3)，
D−|U_J^-|≥c a u⁸−C a u⁹>0。两符号都有足量开可用集C\U_J^±，其中切取质量恰为M的紧孔P_±。例如先取有限内部紧盒面积超过M，再连续分位切取；无须孔形彼此相同。

反单调性给U₁⊂U_J⊂U₀及B(U_J)⊂U_J（J偶）。每个x∈P_±的同一个原半径r=ρ+s，含完整本胞实际孔（r≥ρ>2w），且
\[
 M+\sum_{k\ne0}|(T_\pm k+P_\pm)\cap Q(x,r)|
 \le M+bR_\pm(U_J^\pm;x)
 <M+bF(s)=ar^n.                                      \tag{19}
\]
这是原严格阈值t=a/2^n，全部孔点含边界均成立。完整清单包含负odd支撑，不是只证明pair测试安全。

令Ω_N^±=T_±([-N,N]^n)，保留全部边界切孔，
E_N^±=Ω_N^±\∪_k(T_±k+P_±)。周期积分精确给
|E_N^±|=(2N)^n(d_±−M)>0。r<1、||(T_±)⁻¹||∞<2，所有内层孔点的原见证保留；边界孔总质量O(aN^(n−1))。输入密度点几乎处处进入halo。

Ω_N^±外任一中心、任意正半径的中心方体，分离半空间与中心对称性给其在Ω_N^±内占比≤1/2，故小t时不在严格halo。因此完整极大算子比值为
\[
 \mathscr R_\pm(u,N)=
 \frac{|\mathcal H_n(1-t;E_N^\pm)|}{|E_N^\pm|}
 =\frac{d_\pm}{d_\pm-M}+O_n(a/N).                    \tag{20}
\]
N=ceil(u⁻¹²)使误差o(av⁹)。所有参数为每个充分小实u而非子列；取孔、严格见证和有限化都真实。

## 9. det triangle收益及八阶一致性

特征值直接给
d_+=(1−e)^(n−1)[1+(n−1)e]，
d_-=(1+e)^(n−1)[1−(n−1)e]。
二阶同为−binom(n,2)e²，三阶分别±2binom(n,3)e³，所以
\[
 d_--d_+=-4\binom n3e^3+O_n(e^5).                    \tag{21}
\]
相同M时是精确差
\[
 \frac{d_-}{d_--M}-\frac{d_+}{d_+-M}
 =\frac{M(d_+-d_-)}{(d_--M)(d_+-M)}.
\]
M/a→1、d_±→1，e=κu³，结合(20)即得(5)。n=3时d_+−d_-=4e³恰为等式；κ₀=8，九阶比值差为4av⁹+o(av⁹)。

### 非空的八阶匹配参数供应

现有正剪切真实宏观构造给以下明确参数（本稿沿用其**实际正安全预算**，不靠“已通过”状态）：
\[
 \kappa=\frac{8n}{3}g^2+\frac{320n}{81}g^3u,\quad B_n=\binom n2,
\]
\[
 \Lambda=\frac{B_n\kappa^3}{8n^2g}
 -\frac{5B_n\kappa^4}{96n^3g^2}u+
 \left[\frac{13B_n\kappa^5}{512n^4g^3}
       +\frac{B_n(n-2)\kappa^4}{32n^2g^2}\right]u^2+c_nu^3.
                                                               \tag{22}
\]
其来源为paper/sections/eighth.tex的宏观有限构造、真实两次反馈质量计算，亦见R060实际源§6。该计算保留全局最浅轴在pair外的三面名义量，给
D−|U₂^+|=a u⁸[n c_n g^(n−2)+O_n(1)]。
先在预定紧区间取统一误差，再选大固定c_n，最后小u，因此可取≥c a u⁸；(4)不是空的附加假设。正输入是已给的八阶宏观下界族：
R_+=P_(7,n)+J_n av⁸+O_n(av⁹)。

同一(22)的负剪切U₂^-因(3)安全，给相同的已知八阶系数，同时相对于该正输入多出(5)。(22)中的九阶安全常数影响两个实际比值的共同九阶部分，不影响其差。R060任何保留c a u⁸余量的改进偶层也可用同一比较。微孔未加入本稿；不能把其他轮的微孔常数直接叠加而省去其保护。

所以这是正／负真实输入族之间的九阶改进，不是C_n绝对九阶系数或新全输入上界。若要比较到十阶，则本稿支付掉的F曲率、负j3目标和三边共同源项必须重新精确计算；当前O(a u⁹)坏面积恰在该处变得不足。

## 10. 实算与证明边界

数学证明是§2全整数清单、§3实际全chart预算、§4–6完整核及weighted堆叠误差、§7局部forest共轭和§8真实有限化。有限枚举或符号det不替代这些节点。

已实际读取Terra的Mathematica 15.0.0 ARM：
scratch/r072-computation/determinant.wl/log；
neighbors_batch.wl/log/json；
shared_pair.wl/log及其scope-correction；
full_charts_batch.wl/log、full_charts_results.wl及metadata。
其中n=2,3,4的228个非零{-1,0,1}整数，正balanced/负unipolar分类零mismatch；全整数完整性另由§2解析界证明。det通式残差0，n3 κ=8的固定孔质量归一化det效应为4，但它自身不证明安全。

最早shared_pair积分只取plane extension，曾遗漏target最浅轴在outgoing pair外的同阶条，也将d₂≥s限制延伸为0；历史log保留，scope-correction明确该系数不能代表完整3D cross退款。本稿从未用它替代(18)。

后续完整单相容incoming12/outgoing13的三图为：
q1: y₃∈[0,e]，s≤e(e−y₃)/f，d₁=s，d₂∈[s,e]；
q3: y₁∈[0,e]，s≤e(e−y₁)/f，d₁=s+y₁，d₂∈[s,e]；
q2: y₁,y₃∈[0,e]，s≤(e−y₁)(e−y₃)/f，d₁=s+y₁，d₂=s。
全δ<1、e<w/4保证d₁,d₃≤e，顶端clip不活跃。三个图×K三项共9个独立精确积分给
\[
 \int K_{12}\mathbf1_{V_0^{13}}
 =e^4\left[\frac{7(3+\lambda)}{12}\delta
       -\frac{77+40\lambda}{72}\delta^2
       +\frac{20+21\lambda}{80}\delta^3-\frac{\lambda}{20}\delta^4\right].
\]
这是完整single-F线性图的回归，仍不把它当joint U₀。另以独立Fraction对2组参数×3图×3核项共18项逐项核对WL精确值，全相等。主证明(18)不需该系数；它只帮助核对offmin图和反射方向。精确式、实际输出及哈希另存manifest。所有失败初核及scope修正保留。

本稿尚为送审候选：不签微孔相加、所有safe孔形最优、九阶上确界或一般矩阵符号优化。可得的结构解释是：源预算到当前精度只见单边和两边forest，两个符号精确同价；体积密度却先在triangle的e³项区分，给明确真实比较。更高阶源预算也将开始看triangle，不能无限延长“符号免费”的断言。


