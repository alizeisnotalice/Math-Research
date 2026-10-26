# 固定独立四子集的 rank≥4 aggregate 引理

本轮只替换 076 §3/§9 第 2 项的逐 pattern 常数求和，不重做其 rank3 in-plane、
quotient、signed layer 或全局 Green 证明。主审保留最终证明验收。下文结论自含；
没有计算实际 body32 的常数、全体余项曲率 Γ、g 证书或 Cn。

## 1. 对象和固定选择

设 A=(a_1,…,a_m) 是张成 R⁶ 的有限列阵，正谱权已吸收到实际物理列 a_i。
重复/平行列允许。证明也允许零列（它只能留在未解变量中）；076 的非零列假设
当然是其特例。令

\[
c_t(s)=\frac{t}{\pi(t^2+s^2)},\qquad 0<t\le1.
\]

固定 R>0，取偶 C² 函数 χ，满足 0≤χ≤1、χ=1 于 |s|≤R、χ=0 于 |s|≥2R，
且 |χ′|≤2/R、|χ″|≤6/R²。取 η=1−χ，记

\[
s_t=c_t\chi,\qquad \ell_t=c_t\eta.
\]

R、χ、列编号及以下选择均在任何物理微分之前固定，不依赖 x、t 或微分方向。
pattern L⊂[m] 表示 L 内用 ℓ_t，外用 s_t。其正 pushforward 密度记 p_{t,L}。
以下只考虑

\[
p_t^{\ge4}(x)=\sum_{\operatorname{rank}A_L\ge4}p_{t,L}(x).
\tag{1}
\]

这是平滑参数分解的 rank≥4 部分，不是对 x 事件做条件期望。

令 𝒥 为所有实际独立四列的索引子集，按递增四元组的字典序排序。每个
rank(A_L)≥4 的 L 唯一分给

\[
j(L)=\min\{J\in\mathcal J:J\subset L\},\qquad
\mathcal L_J=\{L:j(L)=J\}.
\tag{2}
\]

对每个 J，另固定 [m]\J 中字典序第一个二元组 K(J)，使

\[
N_J=[A_J\ A_{K(J)}]\quad\text{可逆}.
\tag{3}
\]

这里 J 的四列在前、K 的两列在后，组内各按编号递增。独立集合可在张成 R⁶
的实际列中扩充为基，因此 K 总存在。**同一 J 类的所有 assigned patterns
始终使用同一个 N_J**，包括随后用于扩大非负和的其他 patterns。

## 2. 四条尾与两条 completion 的导数包络

由 Cauchy 密度的显式导数及乘积法则可得

\[
\begin{aligned}
|\ell_t^{(r)}(s)|&\le A_r t,\quad r=0,1,2,\\
(A_0,A_1,A_2)&=\pi^{-1}(R^{-2},4R^{-3},20R^{-4}).
\end{aligned}\tag{4}
\]

三个导数的支撑均在 |s|≥R。具体地，那里
|c_t|≤t/(πR²)、|c_t′|≤2t/(πR³)、|c_t″|≤6t/(πR⁴)，
于是一阶常数为 2+2=4，二阶常数为 6+2·2·2+6=20。
在 |s|<R，ℓ_t 及其前两导数为零；边界由 C² 拼接覆盖。

completion 不能把“导数之和”当作“绝对导数之和”。定义真正的正包络

\[
S_r(t,s)=|s_t^{(r)}(s)|+|\ell_t^{(r)}(s)|.
\]

零阶正好 S_0=c_t；一、二阶只使用不等式

\[
\begin{aligned}
S_1&\le |c_t'|+2c_t|\chi'|,\\
S_2&\le |c_t''|+4|c_t'||\chi'|+2c_t|\chi''|.
\end{aligned}\tag{5}
\]

全空间有 |c_t′|≤1/(πt²)、|c_t″|≤2/(πt³)、c_t≤1/(πt)。因 t≤1，

\[
S_r(t,s)\le C_r t^{-1-r},\qquad
(C_0,C_1,C_2)=\pi^{-1}\bigl(1,1+4/R,2+8/R+12/R^2\bigr).
\tag{6}
\]

所以零阶 completion 不需额外乘 2；高阶保守支付了 χ/η 导数分别取绝对值的成本。
式 (5) 不使用不存在的逐点恒等式 S_r=|c_t^{(r)}|。

## 3. 先微分和取绝对值，再扩大 pattern 和

固定 J，置 B=N_J^{-1}，N 内的已解标量为

\[
w_N=B(x-A_{N^c}z),\qquad z\in\mathbb R^{m-6}.
\]

对任何 L∈ℒ_J，有固定 Lebesgue 积分

\[
p_{t,L}(x)=\frac1{|\det N_J|}
\int\prod_{i\in N} f_{i,L,t}(w_i(x,z))
       \prod_{i\notin N}f_{i,L,t}(z_i)\,dz,
\tag{7}
\]

其中 i∈L 时 f=ℓ，否则 f=s。四个 J 因子始终为 ℓ。物理方向 u 的导数是
在六个已解标量上施加固定算子 Σ_{i=1}⁶(Bu)_i∂_{w_i}。
未解变量及其积分域不动。对固定 t>0，已解因子的 C² 导数全局有界，
未解因子由其可积 Cauchy 密度支配，因此微分和连续性由支配收敛成立。

对一次或两次物理微分按 **有序微分命中索引** 展开 Leibniz 项。
固定某个命中多重指标 α=(α_1,…,α_6)，|α|≤2；方向系数在同一 J 下与 L 无关。
先对每个项取绝对值，再把 assigned 集合扩大为全部 L⊃J。所有项非负，故允许

\[
\sum_{L\in\mathcal L_J}(\text{该绝对项})
\ \le\sum_{L\supset J}(\text{同一 }N_J\text{ 下的该绝对项}).
\tag{8}
\]

右边可以逐标量精确因式分解：

\[
\begin{split}
&\sum_{L\supset J}
 \prod_{i\in N}|f_{i,L,t}^{(\alpha_i)}(w_i)|
 \prod_{i\notin N}f_{i,L,t}(z_i)\\
&=\prod_{i\in J}|\ell_t^{(\alpha_i)}(w_i)|
  \prod_{i\in K(J)}S_{\alpha_i}(t,w_i)
  \prod_{i\notin N}c_t(z_i).
\end{split}\tag{9}
\]

未解变量未被微分，故 s+ℓ=c 精确消去其 χ/η 选择；两个 completion 的选择产生
S_r 而非导数之间的抵消。未解的 Cauchy 概率密度乘积积分恰为 1。若 m=6，
空乘积/零维积分也定义为 1。

式 (8) 的扩大可以让一个 pattern 在不同 J 的上界中出现多次；这只是非负
majorant 的有限过计数，不影响 (1) 原分配的唯一性。不能先把 rank 事件写成
依赖 x 的函数再微分，也不能在微分后改选 N_J。这里两种操作都没有发生。

## 4. 明示的 k=0,1,2 aggregate 常数

定义无 π 的比值

\[
a_1=4/R,\ a_2=20/R^2,\qquad
b_1=1+4/R,\ b_2=2+8/R+12/R^2,
\]

\[
h_0=1,\qquad h_1=\max(a_1,b_1),\qquad
h_2=\max(a_2,b_2,h_1^2),\qquad
q_J=\frac1{\pi^6R^8|\det N_J|}.
\tag{10}
\]

若一次 Leibniz 项的 k 次微分有 q 次命中四条 J 尾，则由 (4)、(6)、(9)，
其 t 次幂是

\[
t^{4}\,t^{-2-(k-q)}=t^{2-k+q}\le t^{2-k}.
\tag{11}
\]

每个尾无论命中 0、1、2 次仍贡献一个 O(t) 因子。比值乘积分别由 h_k 控制：
二阶的可能比值只有 a_2、b_2、a_1²、a_1b_1、b_1²，全部包含于 h_2。

对任意物理方向 u_1,…,u_k（k=0 时空乘积为 1），有

\[
\boxed{
|D_{u_1}\cdots D_{u_k}p_t^{\ge4}(x)|
\le t^{2-k}\sum_{J\in\mathcal J}q_Jh_k
             \prod_{r=1}^k\|N_J^{-1}u_r\|_1 .}
\tag{12}
\]

证明中有序命中系数的绝对值总和，正好是右侧各 ℓ¹ 范数的乘积，故没有遗漏
二阶混合项的 2 或多重指标重数。

任取已认证 Λ_J≥||N_J^{-1}||_{2→1}，则

\[
K_k=\sum_{J\in\mathcal J}q_Jh_k\Lambda_J^k
\quad\Longrightarrow\quad
\|D^kp_t^{\ge4}\|\le K_k t^{2-k}.
\tag{13}
\]

例如 Λ_J=Σ_{i,j}|(N_J^{-1})_{ij}| 是完全有理且安全的选择；也可认证更紧的范数。
这里的 sum 只有独立四子集及其固定 completion，**不含 2^m 个 pattern 常数**。
对 m=32，候选四元组至多 C(32,4)=35960；本轮没有实际枚举 body32 四子集。

## 5. 保留时间幂的统一物理 Hessian 矩阵

给定 0<T≤1，设

\[
R_{\ge4,T}(x)=\int_0^T p_t^{\ge4}(x)\,dt.
\]

(12) 的 k=0,1,2 分别由 t²、t、1 支配，均可在零端积分；同一全空间支配同时
给出 R_{≥4,T}∈C²(R⁶) 和有界导数。这比只证明分布意义余项更强，但仅针对该部分。

固定 u，记 z=N_J^{-1}u、X=Σ_{i=1}⁴|z_i|、Y=|z_5|+|z_6|。不丢弃 (11) 中
额外的 t 次幂，可得到每 J 的二阶绝对 majorant

\[
q_J\{t^2 U(z)+t V(z)+W(z)\},
\tag{14}
\]

其中

\[
\begin{aligned}
U&=a_2\sum_{i\le4}z_i^2+2a_1^2\sum_{i<j\le4}|z_i z_j|,\\
V&=2a_1b_1XY,\\
W&=b_2(z_5^2+z_6^2)+2b_1^2|z_5z_6|.
\end{aligned}
\]

把 (14) 在 (0,T] 积分，再统一使用 2|z_i z_j|≤z_i²+z_j²，定义

\[
d_{\rm tail}=\frac{T^3}{3}(a_2+3a_1^2)+a_1b_1T^2,
\qquad
d_{\rm comp}=T(b_2+b_1^2)+2a_1b_1T^2,
\tag{15}
\]

\[
\boxed{
M_4(R,T)=\sum_{J\in\mathcal J}q_JN_J^{-\mathsf T}
\operatorname{diag}(d_{\rm tail}I_4,d_{\rm comp}I_2)N_J^{-1}.}
\tag{16}
\]

这是在所有方向之前固定的一个 PSD 矩阵，且

\[
-M_4\preceq D^2R_{\ge4,T}(x)\preceq M_4\quad\text{对所有 }x.
\tag{17}
\]

其 Young 余量可显式写成正平方和。以 x_i=|z_i|，对每 J，矩阵对角 majorant
减去 (14) 的时间积分，等于

\[
\frac{a_1^2T^3}{3}\sum_{i<j\le4}(x_i-x_j)^2
+\frac{a_1b_1T^2}{2}\sum_{i\le4<j}(x_i-x_j)^2
+b_1^2T(x_5-x_6)^2\ge0.
\tag{18}
\]

因此 (16) 不是把逐方向最小值拼成矩阵。若要另选更优 completion，应先整体固定
另一套 J↦N_J 并重作证书；不能在 u、x 或 pattern 内偷偷切换。

一个直接的各向同性支付额为 Γ_4=tr M_4，即

\[
\Gamma_4=\frac1{\pi^6R^8}\sum_J\frac{
d_{\rm tail}\sum_{i\le4}\|(N_J^{-1})_{i,:}\|_2^2
+d_{\rm comp}\sum_{i=5}^6\|(N_J^{-1})_{i,:}\|_2^2}
{|\det N_J|}.
\tag{19}
\]

有理 A、R、T 使 π⁶M_4 和 π⁶Γ_4 全部有理；无需浮点优化。Γ_4 **不是完整
Green 的 Γ**：rank3 各种 remainder、signed layer、上尾等仍按 076 另行支付。
固定 face 参数化 E:R⁵→R⁶ 时使用 EᵀM_4E，不能省略物理坐标拉回。

## 6. 可执行证书接口与范围

aggregate.py 不含密度或 kernel 调用，也不在 import/CLI 时执行实际 frame：

- independent_four_bases(A)：按字典序枚举独立 J，再固定首个可逆 K，返回精确
  determinant、N 和 N^{-1}。每个 J 出现一次，不枚举 scalar patterns。
- scalar_constants(R,T)：返回 (10)、(15) 的精确有理常数。
- majorant_from_bases(records,R,T)：在完整已认证 records 下返回 π⁶M_4、π⁶tr M_4、
  π⁶K_0,π⁶K_1,π⁶K_2。该低层函数的“records 覆盖全部 J”是显式输入义务。
- certificate(A,R,T)：自行生成完整 records 后调用上述接口，关闭覆盖义务。

接口只认证给定实际列和 cutoff 的 rank≥4 部分。未来真实 body 参数必须来自同一
model 的实际物理列，给定可用 R、T，并独立支付其他项；本轮没有调用该真实实例。

重复/平行列及零列不会破坏论证：它们不能被错误地作为独立四列或 completion；
留在未解变量时，其 χ/η 正和仍恰为 c_t，积分为 1。不同 rank3 空间的重叠和
S2 接缝也不参与本引理；本结论是固定 cutoff 的全空间 rank≥4 majorant。

本轮仅登记四个新小符号配置：一个非轴 6×8 toy frame 的分配/基扩充，28 个
低阶命中模式的正和恒等式，一个固定 N 的物理 chain rule，以及时间幂/Young
平方和/PSD 拉回。它们是有限代数证据，不能替代以上全称证明；没有重跑 076 的
三项配置，没有实际 body32 Γ、full6D 密度、fullfan、g 或 Cn 计算。
