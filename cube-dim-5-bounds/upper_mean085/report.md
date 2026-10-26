# upper_mean085：正纤维的方向二阶上界与真实盒精确最大坐标矩

本轮从正产品 Cauchy 密度及固定纤维平移独立证明：
\[
\boxed{D_u^2p_t(x)\le {9\over8t^2}\|u\|_\infty^2p_t(x)
\quad(x,u\in\mathbb R^5,t>0).}\tag{1}
\]
对空间积分实际使用的面内位移 \(u_1=0\)，常数进一步改进为
\[
\boxed{D_u^2p_t(x)\le {25\over24t^2}\|u\|_\infty^2p_t(x).}\tag{2}
\]
因此任意正时间段 J=[a,b]（a>0，b 可为∞）、固定 anchor c、真实四尾盒的均值满足
\[
\boxed{\operatorname{mean}_{B}B_J
\le m_J\left[1+{\|\mu\|_\infty\over a}
+{25\over48}{M_2R\over a^2}\right],}
\quad B_J(x)=\int_Jp_t(x)dt,\quad m_J=B_J(c),
\tag{3}
\]
其中 \(\mu=\mathbb E_B(x-c)\)、\(M_2=\mathbb E_B\|x-c\|_\infty^2\)，R 是整个 anchor-to-box 线段族的合法 band 比值上界。这里面内位移的首坐标固定为零；四尾向量的无穷范数与五维位移一致。

可直接复用的 `moment.py` 以 Fraction 返回任意有理盒与任意固定 anchor 的 M2、μ、旧 η 及均值系数，明确支持偏心、边界裁切、anchor 在盒外和零宽度坐标。本批的 16 个冻结盒、10 个有理 lift、7 个符号恒等式全部通过；独立条件最大坐标积分给相同 M2。零个 H、梯度、Green、Arb 或旧基线重算。

对四尾对称半宽 h 盒，\(M_2=2h^2/3\)，故新二阶均值系数为 \(25h^2/72\)，旧系数为 \(7h^2/9\)，新/旧恰为 25/56。更一般地，在同一 R 和同一一次项下，新二阶系数对任何盒均不超过旧系数的 20/21。没有新的完整 I、g、C5 或 Cn 结论；这些是供既有正 band 重新包络的解析和几何组件。

## 1. 合法纤维平移与最小 l1 提升

保持固定模型：16 条未单位化 \(v=(1,\pm1,\pm1,\pm1,\pm1)\)，共同时间 Cauchy 密度 \(s_t(z)=t/[\pi(t^2+z^2)]\)，Lebesgue 及逆 Fourier 因子 \((2\pi)^{-5}\)。令 B 为 5×16 方向矩阵，K=ker B，A=Bᵀ/16。075 已明确证明
\[
BB^T=16I_5,\quad J_B=\sqrt{\det BB^T}=1024,
\quad p_t(x)={1\over1024}\int_K f_t(Ax+w)dw,
\quad f_t(z)=\prod_{j=1}^{16}s_t(z_j).\tag{4}
\]
该纤维公式处处有限、正且连续；可由选取行列式 16 的五个 body 基列，先积分余下十一概率变量证明，而非仅引用 a.e. 版本。

对任意固定 a∈R16，只要 Ba=u，就有 a−Au∈K，平移 w↦w+s(a−Au) 保持 K 上 Euclidean Hausdorff 测度。因此
\[
p_t(x+su)={1\over1024}\int_Kf_t(Ax+w+sa)dw.\tag{5}
\]
提升不必是最小 Euclidean 范数右逆 Au；固定 a 的纤维内平移是合法的。

设 r=||u||∞>0。因为 conv(±V) 正是五维 cube，u/r 可写成 32 个 cube 顶点的正凸组合。一个明确构造是对 e∈{±1}5 取
\[
\lambda_e=\prod_{k=1}^5{1+e_ku_k/r\over2},\qquad
a_v=r(\lambda_v-\lambda_{-v})\quad(v\in V).
\tag{6}
\]
于是 Ba=u、\(\sum|a_v|\le r\)。另一方面
\(||Ba||∞\le\sum|a_v|\)，故实际上 \(\|a\|_1=r\)。u=0 时取 a=0。因为每个 body 向量首坐标为 1，还恰有 \(\sum a_v=u_1\)。若 u1=0，正负系数总量各为 r/2，所以
\[
|a_v|\le r/2\quad\text{对每个 v}.\tag{7}
\]
这里是完整的正凸组合/lift 证明，不是把若干顶点方向的二阶界直接延伸到任意方向，也没有假设 Hessian 为正半定。

## 2. 有界 score 支配与方向二阶常数

令 \(\ell_t=\log s_t\)。其导数为
\[
\ell_t'(z)=-{2z\over t^2+z^2},\quad
\ell_t''(z)={2(z^2-t^2)\over(t^2+z^2)^2},
\quad |\ell_t'|\le1/t,\quad |\ell_t''|\le2/t^2.
\]
对固定 a 和有界 s 区间，score 界给
\(f_t(z+sa)\le e^{|s|\|a\|_1/t}f_t(z)\)。零、一、二阶导由同一原正纤维密度的常数倍支配。对 t≥a_time>0，常数可统一为该正时间下端；075 的 \(p_t\le1/(16\pi^5t^5)\) 又支付无限时间尾。因此可以在 (5) 及所有正 band 内作两次空间求导。没有要求 Cauchy 变量存在任何一阶或二阶绝对矩。

直接的通用证明使用
\[
{s_t''(z)\over s_t(z)}={1\over t^2}{6q-2\over(1+q)^2}\le{9\over8t^2},
\quad q=z^2/t^2,
\]
\[
{9\over8}-{6q-2\over(1+q)^2}={(3q-5)^2\over8(1+q)^2}\ge0.
\tag{8}
\]
对角项使用 (8)，交叉 score 乘积使用 \(1/t^2\)，给
\[
{D_a^2f_t\over f_t}\le t^{-2}
\left[\|a\|_1^2+\tfrac18\|a\|_2^2\right]
\le {9r^2\over8t^2}.
\]
积分即证明 (1)。

为了得到更强的面内常数，令 \(b_i=|a_i|/r\)，\(\sum b_i=1\)。记无量纲 score \(S_i=t\ell_t'(z_i)\)、log 二阶 \(T_i=t^2\ell_t''(z_i)\)。加权 Cauchy–Schwarz 给
\[
{t^2D_a^2f_t\over r^2f_t}
=\left(\sum_i(a_i/r)S_i\right)^2+\sum_i b_i^2T_i
\le\sum_i b_i(S_i^2+b_iT_i).\tag{9}
\]
对任意 β≥0、q≥0，右侧单坐标函数是
\[
g_\beta(q)={(4+2\beta)q-2\beta\over(1+q)^2}
\le K(\beta):={(\beta+2)^2\over4(\beta+1)},
\tag{10}
\]
因为其差恰为
\[
K(\beta)-g_\beta(q)
={[(\beta+2)q-(3\beta+2)]^2\over4(\beta+1)(1+q)^2}.
\]
\(K'(\beta)=\beta(\beta+2)/[4(\beta+1)^2]\ge0\)。一般 b_i≤1 给 K(1)=9/8；面内 (7) 给 b_i≤1/2，所以每项至多 K(1/2)=25/24。结合 (9) 和正纤维积分，得到 (2)。

该 25/24 在上述统一产品密度估计中可达到：取两个相反尾方向的 lift 系数 ±1/2，对应 u=(0,1,1,1,1)，两个 score 符号与系数匹配，且两坐标 z²/t²=7/5。此时加权 Cauchy 与 (10) 均为等号。这不是关于积分后真实 p_t 最优常数的断言。

注意 (1)–(2) 是二阶导数的上界，允许二阶导数为负；不能称 p_t 凸或其 Hessian 正半定。

## 3. 正 band 到真实盒均值

固定 J=[a,b]，a>0。正性及积分求导给
\[
|D_uB_J(c)|\le\|u\|_\infty m_J/a,\qquad
D_u^2B_J(x)\le {25\over24a^2}\|u\|_\infty^2B_J(x)
\quad(u_1=0).
\tag{11}
\]
第一式由 (6) 的最小 l1 lift 和一阶 score 得到，亦与075相同。

R 只须满足 \(B_J(c+su)\le Rm_J\)，对所有真实盒位移 u 和 s∈[0,1]。现有合法选择包括
\[
R=\min\left\{\exp(h_\infty/a),
\exp\left(10\operatorname{asinh}{\sqrt{B_Q^2}\over2\sqrt5a}\right)\right\},
\quad h_\infty=\max_B\|u\|_\infty,
\quad B_Q^2=\max_Bu^TQu.
\tag{12}
\]
这里使用旧基 Q=I4/4−J4/36 的已证平移比及075的全 body 比值；实际应消费两者已有的向上端再取 min。即使 anchor 不在盒内，\(\|su\|∞\le\|u\|∞\)、\((su)^TQ(su)\le u^TQu\)，仍支付整个线段。不能在 anchor 变更后继续使用按旧 anchor 计算的 h∞、BQ² 或 R。

Taylor 的准确积分余项给
\[
B_J(c+u)=m_J+D_uB_J(c)+\int_0^1(1-s)D_u^2B_J(c+su)ds
\le m_J+D_uB_J(c)+{25\over48a^2}\|u\|_\infty^2Rm_J.
\]
先对盒平均一次项，得到 DμB_J(c)，再在均值方向使用 (11)，恰为 (3)。若先平均 |u| 会丢失对称抵消；非对称盒也不能误设 μ=0。上端 b=∞ 的尾段已由正时间支配覆盖，没有空间求导引起的时间端点项。

对每个已有正 band，安全消费形式是
\[
C_J=\min\{R_U,\;1+\|\mu\|∞/a+(25/48)M_2R_U/a^2,\;C_J^{\rm old}\},
\quad \operatorname{mean}H_\epsilon\le\sum_Jm_{J,U}C_J.
\tag{13}
\]
可以继续与合法的全局 cap 截紧。range 的 L、U 不因均值改进而自动改变；最后仍须用 \(z^{5/4}\) 在 [L,U] 的递增凸弦支付真实盒积分，不能反用 Jensen。

## 4. 与旧 η 项的统一比较

对四尾 u，Q=I4/4−J4/36。固定任一坐标 j，令其余三维向量
\(y_k=u_k-u_j/6\)。准确平方完成给
\[
u^TQu={5\over24}u_j^2
+y^T(I_3/4-J_3/36)y\ge {5\over24}u_j^2.
\tag{14}
\]
剩余矩阵在常数方向特征值为 1/6，其余为 1/4，均正。选取 |u_j|=||u||∞，再平均，得到
\[
\eta=\mathbb E u^TQu\ge(5/24)M_2,\qquad
{25\over48}M_2\le{5\over2}\eta
={20\over21}\,{21\over8}\eta.\tag{15}
\]
所以在同一 R 与一次项下，新二阶项对任意有理或实盒都不劣于旧项；这不是只从16个案例外推的结论。079/080 已采用相同的一次项 ||μ||∞；072 的旧 \(\|N^{-1}\mu\|_1\) 至少为 ||μ||∞，故同样可安全取小。

四尾独立均匀对称盒 [−h,h]4 的 max 绝对值分布为 \(P(M\le r)=(r/h)^4\)，0≤r≤h，所以 M2=2h²/3。由 (3) 新系数为 25h²/72，旧 η=8h²/27 给旧系数7h²/9。两者之比25/56，与本机符号资格一致。

## 5. M2 的精确分段算法与退化盒

令每个坐标独立均匀于 [lo_j,hi_j]。正宽度坐标的绝对值分布为
\[
F_j(r)={\operatorname{length}([lo_j,hi_j]\cap[-r,r])\over hi_j-lo_j},
\]
零宽度坐标则是固定值 c_j，CDF 因子为 \(1_{r\ge|c_j|}\)，不做除零。设 R0=max_j max(|lo_j|,|hi_j|)，则
\[
\boxed{M_2=\int_0^{R_0}2r\left[1-\prod_jF_j(r)\right]dr.}\tag{16}
\]
断点为0和所有 |lo_j|、|hi_j|；在每个开区间，正宽度坐标的交长为确定的有理仿射多项式，或为0，固定坐标指示为0或1。乘积次数不超过活动维数，逐系数作有理原函数即可。端点可能的跳跃不影响 Lebesgue 积分。

`max_abs_square_moment(offsets)` 实现 (16)，可选返回每片 CDF 系数与准确贡献。`box_statistics(bounds,anchor_tail)` 先计算真实偏移，返回准确 μ、M2、η、h∞、新旧二阶系数与普通四维体积。无须 anchor 是盒中心或在盒内；边界 clip 由实际 bounds 决定。零宽度时还返回活动维数和活动支撑测度：这是低维独立均匀概率的矩；若目标是四维积分，该盒普通体积为零，应贡献零，不能用低维支撑测度替代四维体积。

独立 `verify.py` 不导入生产的 moment 函数。它把事件按“哪一个活动坐标具有最大绝对值”划分，对该坐标 x 积分
\(x^2/\mathrm{width}_i\) 乘其余活动坐标的 CDF，并加入最大固定坐标压过所有活动坐标的原子贡献。连续活动坐标的并列最大事件为零测；固定坐标并列已统一为一个最大阈值。所有子积分用 SymPy 有理多项式准确积分。此路径与 (16) 的生存函数积分独立，全部16个结果一致。

## 6. 有限资格、实际运行与可调用接口

登记冻结16个案例：079 ordinary00、low00、nearS100、crossing00 的实际偏心 anchor；080 已保存 root 20、8816、52327、52359 的原 anchor；另8个对称、单侧边界、一般非对称、anchor在外、混合退化、全零点、固定各向异性和仅一活动坐标的精确案例。只提取这些旧几何，不运行52360个全层盒，不读取其 H 作为新数学输入。079保存的 offsets、μ、η、h∞、体积逐项准确一致。

另有10个预选有理 u，准确构造 (6) 的 lift 并核 Ba=u、l1=r，以及 u1=0 时逐系数上限r/2。SymPy1.14.0 实际验证7个平方/导数恒等式为零，包括 (8)、(10)、β单调性、两个等号点及 (14)。它们是有限代数核验；全体实参数的符号/正性与纤维推导由本报告证明，不把有限向量检查冒充全域证明。

生产 PID73058、session55601，先登记并armed，发句柄后运行。记录墙钟13.699097459秒，明确包含输入加载和等待启动换行的时间；实际 exit0。独立验证 PID73189，0.798330292秒，exit0。没有程序 wallcap/alarm，符合用户取消限时的当前指令；有限case前缀逐项保存，手动恢复要求源/登记/输入 hash 不变，完成批拒绝重跑。本轮未发生中断恢复。

最小调用方式（所有数字传 Fraction 或准确有理字符串）：
```python
from moment import box_statistics, band_mean_factor
stats = box_statistics(bounds4, anchor_tail4)
factor = band_mean_factor(a, existing_outward_R, stats)
band_mean_upper = existing_positive_mass_upper * factor
```
该函数只接收已经合法的 R≥1 和正时间下端 a>0，不评估 exp/asinh，不读取任何 band，也不生产 H。它返回 range、新 Taylor、旧 Taylor 的最小合法因子。旧质量上端必须来自同一 anchor、同一模型、同一 band；不能混用不同 ε 或时间端点。

完整来源 hash 在 registration 中；结果、独立核验、最小代码及本报告由 manifest 锁定。075/076、旧072/079/080及其它代理文件未修改。仅新增本目录，未提交或推送，未使用 MATLAB。新均值接口已可供082另行登记的旧正band重新包络使用；本批没有运行该全局消费，也没有新 I、g 或 C5 证书。
