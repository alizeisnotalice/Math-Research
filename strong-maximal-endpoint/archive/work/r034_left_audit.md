# r034：左端全控制屏障审计（2026-11-06）

## 结论

对 \(N\ge2\)，题设的
\[
W_t(r)=L_t(r)+B_t(r),\quad
L_t(r)=1-r/N+t\mu/N,\quad B_t=G^t b,
\]
确为所有可预测控制的全局 MAX supersolution，其中 \(G\) 是纯 Exp(1) 核。它给出了有效的终端修补和一个可显式计算的左端初值；不单独给出锐的长时域校正。

## 边界与递推

令
\[
\psi(r)=1-e^{-r}-r/N,qquad
\rho=N+W_0(-Ne^{-N}),
\qquad b(r)=(-\psi(r))_+\mathbf1_{r>0}.
\]
这里必须取 Lambert \(W\) 的主支。\(\psi\) 的两个非负零点为 \(0\) 和
\(\rho\)，且 \(\rho\in(N-1,N)\)。所以 \(b=0\) 于 \(( -\infty,\rho]\)，并在 \(r>\rho\) 为
\[
b(r)=e^{-r}+r/N-1.
\]
由 \(e^{-\rho}=1-\rho/N<1/N\)，有
\(0<b'(r)=1/N-e^{-r}<1/N\) 于 \(r>\rho\)。故 \(b\) 非减且 \(1/N\)-Lipschitz。

终端处，\(r\le0\) 时 \(L_0=\phi_N\)；\(r>0\) 时
\(L_0-\phi_N=\psi(r)\)，因而
\[
L_0+b=\phi_N+\psi_+(r)\mathbf1_{r>0}\ge\phi_N.\tag{1}
\]

在通常的核记号 \(K_cu(r)=\mathbb Eu(r-Z_c)\) 下，用同一 \(Y\sim\mathrm{Exp}(1)\) 耦合
\(Z_c=Y+cJ_c\ge Y\)。因为 \(B_{t-1}=G^{t-1}b\) 非减，
\[
K_cB_{t-1}\le GB_{t-1}=B_t.\tag{2}
\]
又 \(m(c)=\mathbb EZ_c\le\mu\)，所以
\[
K_cL_{t-1}=1-r/N+(t-1)\mu/N+m(c)/N\le L_t.\tag{3}
\]
相加即 \(K_cW_{t-1}\le W_t\)，对所有动作成立。结合 (1) 可按反向归纳支配全控制 Bellman 值。

## 闭式核验

令 \(x=r-\rho\)。对于 \(t\ge1\)，若 \(x\le0\)，\(B_t(r)=0\)。若 \(x>0\)，对 Gamma 密度在 \([0,x]\) 积分得到
\[
B_t(r)=e^{-r}\frac{x^t}{\Gamma(t+1)}
+(r/N-1)F_t(x)-\frac tN F_{t+1}(x),\tag{4}
\]
其中 \(F_j\) 是 Gamma\((j,1)\) CDF。关键恒等式是
\(\mathbb E[G_t\mathbf1_{G_t\le x}]=tF_{t+1}(x)\) 与
\(\rho/N-1=-e^{-\rho}\)。式 (4) 在 \(x\downarrow0\) 连续归零。\(t=0\) 应单独定义为 \(B_0=b\)；除非另行采用退化 Gamma 的 CDF 约定，不能把 (4) 原样当作 \(t=0\) 公式。

## 取最小值与拼接的限制

两个全局 supersolution \(W^{(1)},W^{(2)}\) 的逐点最小值仍是 supersolution：
\[
\sup_cK_c(W^{(1)}\wedge W^{(2)})
\le\min\{\sup_cK_cW^{(1)},\sup_cK_cW^{(2)}\}
\le W^{(1)}\wedge W^{(2)}.
\]
且两者都须在终端支配 \(\phi_N\)。因此可合法把本屏障与一个**已全局验证**的锐校正屏障取最小值。

但局部校正、数值拟合、只在一个状态窗成立的比较，或未修补终端的仿射函数，都不能直接参加该最小值拼接；否则 supersolution 不等式会在离开窗口或终端处失效。本屏障本身只把 \(L_0\) 在 \(r>\rho\) 的终端缺口补齐；它没有提供另一屏障，也不自动产生 \(g\log N\) 级的锐校正。任何声称的拼接收益仍须给出第二个全局屏障及其与 (4) 的交叉区域验证。

## 补充：自然管两侧的全控制边界

以下是新的粗全策略边界，而不是锐拼接结论。令
\[
r_-=t\mu-A\sqrt{t\log t},\qquad
r_+=t\mu+A\sqrt{t\log t},\qquad N\ge t,
\]
其中 \(A>0\) 固定、\(t\to\infty\)。设 \(U_t\) 为全控制 Bellman MAX 值，而 \(v_t\) 为固定 \(c=\mu\) 参考值。则一致于 \(N\ge t\)，
\[
0\le L_t(r_-)+B_t(r_-)-U_t(r_-)
\le t^{-A^2/(2\mu)+o(1)}+Ce^{-ct},\tag{5}
\]
并且
\[
U_t(r_+)\le t^{-A^2/(2\mu)+o(1)}.\tag{6}
\]

为核对 (5)，写终端为 \(\phi_N=L_0-a+b\)，其中
\(a=\psi_+\mathbf1_{r>0}\le1\)。因为 \(U_t\ge v_t\)，且参考和的 b 项非负，
\[
L_t+B_t-U_t\le L_t+B_t-v_t
\le\Pr(S_t^\mu<r_-)+B_t(r_-).\tag{7}
\]
固定 \(c=\mu\) 时 \(\operatorname{Var}Z_\mu=\mu\)，并且负参数矩母函数在零处的展开为
\(\log M_\mu(-\theta)=-\mu\theta+\mu\theta^2/2+O(\theta^3)\)。以
\(\theta=A\mu^{-1}\sqrt{\log t/t}\) 作 Chernoff，得
\[
\Pr(S_t^\mu<r_-)
\le t^{-A^2/(2\mu)+o(1)}.\tag{8}
\]
又 \(\rho>N-1\) 且 \(r_-\le t\mu\)，所以在 \(B_t\) 的非零情形
\[
x=r_--\rho\le t\mu-(N-1)\le t(\mu-1)+1.
\]
并且 \(b(u)\le u/N\le\mu\) 于积分支撑上。因此
\[
B_t(r_-)\le\mu\Pr(G_t<x)\le Ce^{-ct},\tag{9}
\]
因为 \(\mu-1<1\) 且 Gamma 下尾在任意固定小于一的线性位置指数小。这里常数独立于 \(N\ge t\)。

为核对 (6)，令
\[
M_*(\theta)=\sup_{c\in[0,\infty]}\mathbb Ee^{\theta Z_c}.
\]
在 \(\theta\downarrow0\) 时，最大点位于 \(c=\mu+O(\theta)\)，故
\[
\log M_*(\theta)=\mu\theta+\frac\mu2\theta^2+O(\theta^3).\tag{10}
\]
这是全控制条件矩母函数界，并非错误的 c-随机序。若
\(1/N\le\theta\le1\)，终端的两支分别满足
\(\phi_N(r-s)\le e^{\theta(s-r)}\)：在 \(s<r\) 用 \(\theta\le1\)，在 \(s\ge r\) 用
\(1+(s-r)/N\le e^{(s-r)/N}\le e^{\theta(s-r)}\)。逐步条件化于是给
\[
U_t(r)\le\exp\{-\theta r+t\log M_*(\theta)\}.\tag{11}
\]
在 \(r=r_+\) 取同一 \(\theta=A\mu^{-1}\sqrt{\log t/t}\)，它对充分大 \(t\) 自动满足 \(\theta\ge1/N\)（因 \(N\ge t\)），将 (10) 代入 (11) 即得 (6)。

式 (10) 需要作为独立的、在小正 \(\theta\) 一致的矩母函数引理引用或证明；一旦具备，它和上述估计都不提供 \(g\log N\) 级的尖锐系数。
