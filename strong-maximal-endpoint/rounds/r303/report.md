# 第303轮　任一删轴方向上的近门槛反向定位

2026-11-12。对有限正体积指标源 $E\subset\mathbb R^n$ 及固定中央紧区间 $p\in I$，写 $n=d+k,\ 1\le k\le\sqrt d$。定义实际超额 $r=\log(|\{M_S1_E>p\}|/|E|)-B_n(p)$。若 $r\ge0$，则对每一种保留 $d$ 轴的选择均存在 $q$，使
\[
 |q-p|\le C_I\frac{\sqrt{s\log d}+\log d}{\sqrt d},\quad
 \log\frac{\int|\{M_d1_{E_z}>q\}|dz}{|E|}-B_d(q)
 \ge r-C_I(k/s+k/\sqrt d)-d^{-1}.
\]
参数满足 $2\le s\le c_I\sqrt d$ 且 $k/s\le\log d$。这是同一个源在指定方向上的实际截面积分；门槛 $q$ 可随方向改变。特别是固定 $k$ 取 $s=\log d$ 时，损失 $o(1)$ 且门槛只移动 $O_I(\log d/\sqrt d)$。

**反用幂包络。** 沿用[第302轮](../r302/report.md)的 $L=-B_d'(p)$、$\delta=s/L$、$b=p-\delta$。设 $V(q)=\int|\{M_d1_{E_z}>q\}|dz$。弱矩传播给
\[
 |\{M_S1_E>p\}|\le L_s^{sk}\sup_{b<q<1}((q-b)/\delta)^s V(q),\quad s\log L_s\le\mu+C/s.
\]
选择右侧达到上确界的 $e^{-1/d}$ 倍以上的门槛。上确界有限，因为 $b$ 在固定紧区间、$a\le1$ 且 $V(b)$ 有限。令 $x=(q-p)/\delta$。在 $q\le p+\varepsilon$，第302轮 Taylor 比较保留负余项给
\[
 s\log(1+x)+B_d(q)-B_d(p)\le -c_I s\min(x^2,|x|).
\]
对 $-1<x\le1$ 用二次界，对 $x\ge1$ 用线性界，先定中央邻域和 $\varepsilon$，再取 $c_I$ 足够小。该段一般指标粗界 $V(q)\le |E|\exp(B_d(q)+C_I\log d)$ 由[第110轮的秩估计](../r110/central110.md)所引第91轮全局包络逐截面积分得到。对 $q>p+\varepsilon$ 使用 $V(p+\varepsilon)$；基线损失至少 $c'_I\sqrt d$，$s\log(1+C_I\sqrt d/s)$ 在 $s\le c_I\sqrt d$ 上至多为该损失的一半。另有辅助坐标费用 $O_I(k/s+k/\sqrt d)$。若 $k/s$ 很大，这一费用不能忽略。

由假设 $k/s\le\log d$，远端不可能达到所选值；近端必须满足 $s\min(x^2,|x|)\le C_I\log d$。代回 $q-p=\delta x$ 得所述窗口，保留非正的基线损失即得截面超额下界。更一般地，$k/s=o(\sqrt d)$ 时用 $\log d+k/s$ 替换窗口公式的两个 $\log d$，证明相同。

**按源质量定位。** 在上述门槛，令 $D_z=|\{M_d1_{E_z}>q\}|/(|E_z|e^{B_d(q)})$，忽略零截面，按概率 $|E_z|dz/|E|$ 抽取 $z$。若上面得到的对数下界为 $h$，则 $\mathbb E D_z\ge e^h$ 且 $D_z\le d^{C_I}$，从而 $\Pr(D_z\ge e^h/2)\ge e^h/(2d^{C_I})$。这一结论对每个删轴方向成立，但不保证不同方向的高超额截面经过同一源点。一般端点系数差仍为 $(\mu-1)/8$。下一步检验这些按质量计量的例外截面在不同方向上的交叠，并研究是否能在随机删轴时控制门槛漂移的符号。
