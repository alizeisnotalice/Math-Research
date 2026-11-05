# 第21轮：一般上界的对数级余项

日期：2026-11-06（Asia/Shanghai）。状态：一般输入上界；使用精确非线性Gaussian移位、空间矩比较和热半群积分。

## 1. 主要结论

令mu=log B、sigma=sqrt(mu)、L=log n，以标准正态上尾定义barPhi_G(z_n)=sigma z_n/sqrt n，d_n=sigma sqrt n z_n。
存在固定K_0使充分大n时
\[
\boxed{\log C_n\le\log\mathcal R_n\le n\mu-d_n+K_0\log n.}\tag{1}
\]
结合第19轮的实际输入下界，
\[
 n\mu-d_n+\frac{9-\mu}{24}\log n+o(\log n)
 \le\log C_n\le n\mu-d_n+O(\log n).\tag{2}
\]
一般上下界的对数差距由第20轮的O((log n)²)进一步降为O(log n)。最优log n系数及一般乘法渐近等价仍未闭合。
本轮只将第20轮屏障的全带统一曲率补偿替换为状态相关补偿；一般输入桥仍使用第11轮已证明的C_n≤mathcalR_n。

## 2. 固定时域与光滑截断

固定大整数N，L=log N≥1，取A=10，t_0=C_0 L。参数K≥1及C_0在下文依次选取，均不依赖N。
令a=A sqrt L，定义非负偶函数psi_L(z)：|z|≤2a时为z²；|z|≥3a时为44a²/7；中间令v=|z|/a−2，并取
\[
 \psi_L(z)=a^2P(v),\quad
 P(v)=4+v^2+4v-10v^4+8v^5+v^6-12v^7/7.
\]
P'(v)=2(v+2)(1−v)³(1+3v+6v²)≥0，端点前三导数匹配，所以psi_L为C³，且
\[
 0\le\psi_L\le C_A L,\quad
 |\psi_L'|\le C_A\sqrt L,\quad
 |\psi_L''|\le C_A,\quad |\psi_L'''|\le C_A/\sqrt L.\tag{3}
\]
对t∈[t_0,t_0+N]及x∈R，置
\[
 H(t,x)=K\left[\psi_L\left(\frac{x}{\sigma\sqrt t}\right)+2\log(t/t_0)\right],\quad
 y=x-H(t,x),\quad W(t,x)=\Phi_G(y/(\sigma\sqrt t)).\tag{4}
\]
H非负且全局≤C_{A,K}L。继续增大C_0可使|H_x|≤1/4、|H_xx|≤C_AK/t、|H_xxx|≤C_AK/(t^(3/2)sqrt L)。因此W_x>0；不以线性展开Phi−Hf替代W。

## 3. 精确热残差

记F(t,y)=Phi_G(y/(sigma sqrt t))，f(t,x)=F_y(t,x−H(t,x))，以及mathscrL=partial_t−mu partial_xx/2。
链式法则精确给
\[
 \frac{\mathscr L W}{f}
 =-(H_x-H_x^2/2)y/t+\mu H_{xx}/2-H_t.\tag{5}
\]
在|z|≤2A sqrt L、z=x/(sigma sqrt t)内，psi_L=z²，所以(5)精确等于
\[
 -\frac Kt-\frac{Kx^2}{\mu t^2}
 +\frac{2KxH}{\mu t^2}+\frac{2K^2x^2(x-H)}{\mu^2t^3}.\tag{6}
\]
后两项除以K(1+z²)/t的绝对值，由
\[
 C_A\{K/\sqrt{C_0}+K^2/C_0\}
\]
控制。这里使用log(t/t_0)/sqrt t≤2/(e sqrt t_0)、log(t/t_0)/t≤1/(e t_0)。
选择固定C_0足够大，得到整个较宽内带
\[
\boxed{\mathscr L W(t,x)\le-\tfrac12K\frac{1+z^2}{t}f(t,x).}\tag{7}
\]
本估计按移位后的密度f计量；不要求H|x|/t小，也不把移位密度与未移位密度比较。
全域由(3)、(5)、sup_y |y|F_y=O(1)、sup_y F_y=O(t^-1/2)得|mathscrL W|≤M_{A,K}，其常数独立于N。

## 4. 统一空间矩比较

沿用第20轮的一步变量Z_r，D=mu−E Z_r≥0，D≤mu−1，|Var Z_r−mu|≤C_m sqrt D，密度≤2e^-z。它与参考G~Normal(mu,mu)均有统一加权五阶矩。
在t≥t_0、当前|z|≤A sqrt L附近，固定时间的空间导数满足
\[
 W_x\asymp f,\qquad
 |W_{xx}/W_x|\le C_A\{(1+|z|)/\sqrt t+K/t\}.\tag{8}
\]
对任意空间位移u，利用y(x+u)−y(x)的全局Lipschitz界及(3)，
\[
 \frac{|W_{xxx}(t,x+u)|}{f(t,x)}
 \le\frac{C_A}{t}(1+z^2+u^2/t)e^{\eta_0|u|},\qquad\eta_0=1/8,\tag{9}
\]
只要C_0≥C_A(1+K)²足够大。常数C_A可选为与K及C_0进一步增大无关。
证明(9)时使用W_xxx=F_yyy y_x³+3F_yy y_x y_xx+F_y y_xxx，以及
F_y(t,y(x+u))/F_y(t,y(x))≤exp(|y(x)||y(x+u)−y(x)|/(mu t))。
内带|y|/(sigma sqrt t)≤C_A(1+|z|)，于是指数率可用C_0压到eta_0；K/sqrt t和K²/t产生的系数同样被C_0≥C_A(1+K)²吸收。

在时间t−1对Y=mu−Z_r及Y_0=mu−G作二阶Taylor并用(9)控制三阶积分余项。
主项差为D W_x+(Var Z_r+D²−mu)W_xx/2。按第20轮同一完成平方，先吸收D²项，再优化sqrt D项，由(8)—(9)得到
\[
 E_rW(t-1,x+\mu-Z_r)
 \ge P_1W(t-1,\cdot)(x)-C_*\frac{1+z^2}{t}f(t,x),\tag{10}
\]
P_s为方差mu s的Gaussian热半群。时间t−1与t处的f可比；这是下一节短时间对数导数界的直接特例。
C_*可取独立于K，只要求C_0先满足随K的上述充分大条件。这样后续先选K再增大C_0无循环。

## 5. 热半群积分给离散正余量

由直接微分P_sW(t−s)，有精确恒等式
\[
 P_1W(t-1,\cdot)(x)-W(t,x)
 =-\int_0^1P_s[\mathscr LW(t-s,\cdot)](x)\,ds.\tag{11}
\]
W及相关导数有界，(3)确保此微分和积分交换合法。
对当前|z|≤A sqrt L，较宽的2A带提供缓冲：Gaussian一步逃出|xi|≤2A sigma sqrt(t−s)L^(1/2)的概率≤C exp(−c_A t L)，对0≤s≤1一致。
在其补集，(7)使mathscrL W非正；不能丢掉该带内的|U|>1部分，因为其概率不是小量。
为取得严格正下界，只取|U|≤1这一固定正概率子事件。此时xi=x+sigma sqrt s U，以下量一致可比：
\[
 1+\xi^2/(\mu(t-s))\asymp1+z^2,
 \qquad f(t-s,\xi)\asymp f(t,x).\tag{12}
\]
后者可由短时间对数导数直接检查：在内带
|partial_x log f|≤C_A/sqrt C_0，
|partial_t log f|≤C_A\{C_0^-1+K C_0^-3/2\}。
这些界按移位密度计算，常数不含L；x只移动O(1)、t只变化≤1。扩大C_0后，可比常数固定不依赖K。

外逃区以第3节的全域界M_{A,K}处理。因为t≥C_0 L，M_{A,K}exp(−c_A t L)≤N^-5/4对充分大N成立。
故(7)、(11)—(12)给固定c_*>0下的离散收益
\[
 P_1W(t-1)(x)-W(t,x)
 \ge c_*K\frac{1+z^2}{t}f(t,x)-\tfrac14N^{-5}.\tag{13}
\]
先选K≥2C_*/c_*+2，再增大C_0使第2—5节的所有条件成立。结合(10)，当前中央带得到
\[
 E_rW(t-1,x+\mu-Z_r)\ge W(t,x)-\tfrac14N^{-5}.\tag{14}
\]

## 6. 外带和初值

令B_0=4sigma sqrt C_0，离散状态x=q−mu k−B_0 L，t=k+t_0。定义W_k(q)=W(t,x)、epsilon=N^-4、delta=N^-5及V_k=W_k−epsilon−k delta。
左外带z<−A sqrt L，由H≥0，W_k≤Phi_G(z)≤N^-50≤delta，而E W_(k−1)≥0。
右外带z>A sqrt L，选择C_0使H≤C_{A,K}L≤x/4。
在Z_r≤mu+x/4事件上，上一步移位后的坐标至少x/2，因此
\[
 E_rW_{k-1}(q-Z_r)\ge1-2e^{-x/4}-e^{-A^2L/8}\ge1-\delta
\]
对充分大N成立；这里x≥sigma A sqrt C_0 L。选C_0使第一个尾指数至少8，第二个指数为12.5。
结合(14)，全实轴及所有r均有E_r W_(k−1)(q−Z_r)≥W_k(q)−delta。

q<0时初值x=q−B_0L且H≥0，故W_0(q)≤Phi_G(−4sqrt L)≤epsilon；q≥0时W_0≤1。
所以V_0≤1_{q≥0}，且E_rV_(k−1)(q−Z_r)≥V_k(q)。如第20轮，条件反向迭代允许策略依赖完整过去。
又因H≤C_{A,K}L，最终得到固定C_H下的全状态分布界
\[
\boxed{P_\pi(S_N\le q)\ge
 \Phi_G\left(\frac{q-N\mu-C_HL}{\sigma\sqrt{N+t_0}}\right)-2N^{-4}
 \quad(q\in\mathbb R).}\tag{15}
\]

## 7. 一般输入账本与最后比较

第20轮对CDF下界的有限区间积分可以原样使用：先积分到T=(mu+1)N，再用一致正mgf处理无穷远尾。
由于(15)的Gaussian均值只增加O(L)，对s=Nmu−d、d=d_N+O(L)，得到
\[
 E_\pi(S_N-s)_+\le d+O(L).\tag{16}
\]
Gaussian负部尾积分为O(sqrt L)，已被O(L)吸收。没有假设原自适应总和的方差。

设可行非减轮廓h有J≤e^-s、E_pi h(S_N)>1。将第20轮中段截断长度改为b_N=sqrt N L。
模数反解给(h−1)_+≤(1+O(b_N/N))(t−s)/N于中段，远尾低振幅贡献≤exp(−cL²)，高振幅贡献≤exp(−c'N)。
由(16)及基准1分解，得到
\[
 P_\pi(S_N\le s-2L)\le d/N+O(L^{3/2}/N).\tag{17}
\]
选择s=Nmu−d_N+K_0L，K_0>C_H+3，(15)在q=s−2L给距离至多d_N−L、方差至少mu N。
因此
\[
 P_\pi(S_N\le s-2L)\ge
 \overline\Phi_G\left(z_N-\frac L{\sigma\sqrt N}\right)-2N^{-4}
 \ge d_N/N+cL^2/N.\tag{18}
\]
因为phi(z_N)/(sigma sqrt N)~L/N。(17)至多d_N/N+O(L^(3/2)/N)，与(18)矛盾。
故log mathcalR_N≤s；第11轮一般输入桥给(1)。

## 8. 本机计算与证明边界

Wolfram15实际核验(5)—(6)的精确非线性热残差，符号残差为0；光滑截断P的六个C³匹配残差均为0。
MATLAB R2026a对K=100、C_0=10^10，L=10,100,1000,10^6及t/t_0=1,2,100,10^6扫描精确热残差；|z|≤10sqrt L的归一化正余量最小约0.96478。
此外MATLAB使用64/128点Gauss–Laguerre积分实际一步核，用四点Gauss–Legendre计算稳定的CDF差；L=10,20,40，z=−sqrt L,0,sqrt L，并优化动作c及端点。
最小归一化一步余量约2.51e-10，64/128点差至多4e-13；均为浮点诊断，非统一证书。诊断参数不被宣称为存在性证明中的认证常数。
Wolfram另用80位自适应积分对L=10的三个固定mu动作独立复算实际一步差，检验稳定差值实现。
复现：MATLAB `addpath('code'); sm_shift_barrier_checks;`；Wolfram `wolframscript -file code/shift_barrier_checks.wl`。
三个固定mu动作的MATLAB/Wolfram归一化一步余量差至多1.67e-13；MATLAB Code Analyzer无报告。算法确定性，无随机种子；严格结论由第2—7节分析证明给出，浮点积分不承担全状态或全策略保证。

## 9. 差距与下一轮

一般log C_n−nmu+d_n目前被两个正的常数倍log n夹住；其最优系数未知，乘法渐近等价仍未闭合。
下一轮两项动作：① 保留状态相关残差的具体系数，研究能否确定反馈上界模型的log n项，而非只取充分大的K；② 优化真实输入的平台软段并与上界模型分开登记，判断现有下界系数是否可再提高。
本轮额外热半群及空间导数估计已直接证明，基础来源沿用sources.md；无新增优先权主张，无冗余跟踪材料需清理。
