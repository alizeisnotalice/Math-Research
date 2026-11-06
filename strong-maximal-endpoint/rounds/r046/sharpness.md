# 非积角条带：精确预算、倾斜界及稳定核

日期：2026-11-07。沿用 [proof.md](proof.md) 的记号。

## 真正非积的预算条带

取两块 \(b=N/2,c=\mu\)，\(S=T_1+T_2,D=T_1-T_2\)，固定 \(\lambda>0,w=4\)。令
\[
 d_0=\lambda\sqrt{NL},\quad d_1=d_0+w\sqrt N,\quad
 s_0=\tau+\log(w\sqrt N/2),\quad
 H_N=1_{\{S\ge s_0,\ d_0\le D\le d_1\}}.
\]
充分大的 N 有 \(s_0>d_1\)，正象限约束不活跃。Jacobian 为1/2，故
\[
 J(H_N)=\frac12\int_{d_0}^{d_1}\int_{s_0}^{\infty}e^{-S}dS\,dD=e^{-\tau}.
 \tag{1}
\]
H_N 不分离变量，但仍是实际预算输入。

设 \(\Lambda(\eta)=-\log(1-\eta)+\log[(1+e^{(\eta-1)\mu})/(1+e^{-\mu})]\)，
\(\Psi(\theta)=b[\Lambda(\theta)+\Lambda(-\theta)]\)。
按 D 倾斜，取唯一正根
\(\Psi'(\theta)=d_0+2\sqrt N\)。则
\(\theta=(\lambda/\mu)\sqrt{L/N}+O(N^{-1/2})\)。

在该概率 Q 下，D 均值在条带中点。对 \(\eta<1\)，令
\(p_\eta=(1+e^{\mu(1-\eta)})^{-1}<1/2\)，有
\[
 \Lambda'''(\eta)=\frac2{(1-\eta)^3}
 +\mu^3p_\eta(1-p_\eta)(1-2p_\eta)>0.
\]
因此 \(\Lambda'\) 凸，\(E_QS\ge N\mu\)。
在 \(|\eta|\le1/4\) 用 \(\mu<13/10\) 得
\(\Lambda''(\eta)\le16/9+169/400=7921/3600<9/4\)。
所以 \(\operatorname{Var}_Q D,\operatorname{Var}_Q S\le9N/4\)。
当以下有限条件都成立：\(\theta\le1/4,s_0>d_1,d_N>\log(2\sqrt N)\)，且
\[
 q_N=1-\frac9{16}-\frac{9N}{4[d_N-\log(2\sqrt N)]^2}>0,
\]
Chebyshev 并界给 \(Q(H_N=1)\ge q_N\)，且 \(q_N\to7/16\)。
由反向似然与普通 Chernoff，得到可直接计算的界
\[
 q_N e^{\Psi(\theta)-\theta d_1}
 \le EH_N\le \min\{1,e^{\Psi(\theta_U)-\theta_Ud_0}\},
 \quad\theta_U=\frac\lambda\mu\sqrt{L/N}.\tag{2}
\]
上界还要求 \(\theta_U<1\)，对本轮全部数值点均成立。
两端 log 除以 L 后都趋于 \(-\lambda^2/(2\mu)\)，故
\[
 EH_N=N^{-\lambda^2/(2\mu)+o(1)}.\tag{3}
\]
(1)—(3) 给出非积的角向预算障碍，但 EH_N<1，不构成真实 C_N 的下界改进。

## 稳定密度与概率实现

对 \(t>0,c>0\)，解 \(t/c=bm(d)/d\)。
函数 \(m(d)/d=1/d+(1+e^d)^{-1}\) 严格递减，根位于
\([bc/t,2bc/t]\)。令 \(\rho=d/c,\eta=1-\rho\)。
倾斜后 \(\rho Z_c\) 精确同分布于 \(Z_d\)，故
\[
 \log f_{b,c}(t)=b\!\left[\log(c/d)+\log(1+e^{-d})-\log(1+e^{-c})\right]
 -\eta t+\log\rho+\log g_{b,d}(\rho t-bm(d)).\tag{4}
\]
这里 g 是中心化 b 重卷积密度。数值实现保留实际残差 offset；不把近似根当精确根。
由中心化特征函数进行归一化 Fourier 积分，避免直接幂和溢出。
若截断归一化频率 \(|v|\le M\)，绝对截尾误差有解析上界
\[
 \pi^{-1}(1+M^2/b)^{-(b-2)/2}\arctan(\sqrt b/M).
\]
它只控制积分尾，未认证梯形离散、求根和舍入的总误差。

尾概率通过同一倾斜 Fourier 积分乘 resolvent
\((a\pm iv/\sqrt b)^{-1}\)，其中 \(a=|\eta|/\rho\)。
正号用于正倾斜的上尾；负号用于负倾斜的下尾。
当前程序只在相应倾斜符号匹配、远离零的测试点调用，未实现均值处的零倾斜极限。
小维数另以有限 Gamma–二项混合独立核对。
条带概率则对第一块积分，第二块使用 GammaRegularized 的区间形式，避免两 CDF 相减。

代码：[密度核](../../code/block_kernel46.wl)、[象限尾概率](../../code/quadrant46.wl)、
[条带积分](../../code/stripe46.wl)、[解析界](../../code/stripe_bounds46.wl)、
[恒等式及加精度检查](../../code/verify46.wl)。本轮未使用随机数。
