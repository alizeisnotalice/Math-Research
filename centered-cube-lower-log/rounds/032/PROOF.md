# 分层原子乘积及其任意正混合的强最大值预算

本稿在随机jitter与近重叠双峰的本机检查之后提出。它扩展第30轮的局部簇证明，但不调用该轮的小球或方差下界。因此无需有限域标签、Sidon资格、全窗口交叉引理或原安全尺度的平坦方差。

## 定理与适用范围

固定 \(p_*\in(0,1)\)、\(A<\infty\)。取整数 \(K=K_n\to\infty\)，令
\[
 v=\frac{9p_*^2n}{29K^2}\ge1,
 \qquad \frac{(1+\log n)^2}{K}\longrightarrow0.              \tag{1}
\]
在第 \(i\) 个坐标，任取 \(0\le p_i\le p_*\)、每格一个原子
\(a_{ij}\in[j/K,(j+1)/K)\)、\(0\le M_i<K\)、实系数 \(|c_{im}|\le A/m^2\)，以及任意有限 \(w_i\ge4\)。置
\[
 \delta_i=(w_iKv)^{-1},\quad
 D_i(x)=\sum_{m\le M_i}\pi m c_{im}\cos(2\pi mx),
\]
\[
 \rho_i(x)=1-p_i+\frac{p_i}{K}D_i(x)
                 +\frac{p_i}{K}\sum_{j=0}^{K-1}\beta_{\delta_i}(x-a_{ij}).
\]
沿固定的 \(K_n\) 序列，当 \(n\) 充分大时，这些密度非负、每周期质量1，门槛对所有 \(p_i\)、原子、系数、宽度与混合选择一致；门槛允许取决于(1)中的小量条件沿此序列成立的位置，不主张对任意延迟增长的序列共用一个数值 \(n_0\)。令 \(F=\prod_i\rho_i\)，并令 \(I=[10/9,11/9]\)。存在只依赖 \(p_*,A\) 的常数 \(C\)，使
\[
 \int_{\mathbb T^n}\sup_{L\in I}A_LF\,dP^{\otimes n}
       \le C(1+Kv).                                       \tag{2}
\]
对这些乘积密度的任意可测概率正混合，(2)仍成立，允许混合中的原子、系数、\(p_i\)、宽度、各坐标平移均不同；均须满足共同 \(K,v,p_*,A\) 的包络。特别地，\(K\asymp\sqrt{n/\log n}\) 时，弱尾系数至多 \(C\sqrt{n\log n}\)。

这是固定窗口、固定分层几何与系数包络的族内上界。它不是一般 \(C_n\) 上界，仍允许不超过半次对数的族内增长，也不证明这种增长可以达到。

## 1. 不需要标签的共同预算

每格一个原子使未平滑的周期原函数
\(P_i(x)=\#\{a_{ij}\le x\}-Kx\) 在一个周期内满足 \(|P_i|\le1\)。箱卷积保持此界，故
\(H_{i,L}=P_{i,\delta_i}(x+L/2)-P_{i,\delta_i}(x-L/2)\) 满足 \(|H|\le2\)。精确平均公式为
\[
 Y_{i,L}=A_L\rho_i=1+\frac{p_i}{KL}(H_{i,L}+B_{i,L}),
 \quad B_{i,L}=\sum_{m\le M_i}c_{im}\sin(\pi mL)\cos(2\pi mx).
\]
由系数包络，\(\|B\|_\infty\le A\sum m^{-2}\)，\(\|D_i\|_\infty\le\pi A(1+\log K)\)。所以共同背景至少
\[
 1-p_*-\frac{p_*\pi A}{K}(1+\log K)\ge(1-p_*)/2>0.
\]
并且在 \(J_2=[10/9,4/3]\) 上最终 \(1/2\le Y\le3/2\)、\(|\log Y|\le C p_*/K\)。所有平均准确积分为1。

若 \(t=|L-J|<2/K\)，未平滑长度增量是两个各长 \(t/2\) 的扫入计数之和减 \(Kt\)。每个短区间至多遇两个原子，故总计数 \(N\le4\)、\(EN=Kt\)，给
\(\|H_{i,L}-H_{i,J}\|_2^2\le4Kt\)。箱卷积是正概率 \(L^2\) 收缩，此界不变。低频部分由 \(\sum m^2c_{im}^2\le A^2\sum m^{-2}\) 给统一Lipschitz \(L^2\) 界。

无需箱峰互不重叠。只要 \(\delta_i\le1/K\)，任意一点附近宽度 \(\delta_i\) 的区间至多碰到三个格子，峰重数至多3。因此
\(\|P'_{i,\delta_i}\|_\infty\le3/\delta_i+K\)。它只把第30轮的导数常数放大固定倍。令 \(\varepsilon=p_*/K\)、\(\xi_{i,L}=\log Y_{i,L}\)。任意确定锚点 \(a\) 下的律 \(dQ_{i,a}=Y_{i,a}dP\) 满足
\[
 E_{Q_{i,a}}(\xi_{i,L}-\xi_{i,J})^2\le C\varepsilon^2Kt,
 \qquad \|\xi_{i,L}-\xi_{i,J}\|_\infty
        \le C\varepsilon\min(1,t/\delta_i).                 \tag{3}
\]
前一式由固定 \(Y_{i,a}\le3/2\) 改变积分测度，不需 \(Q_{i,L}\)。

## 2. 非同分布坐标的局部指数矩

先设全部 \(4\le w_i\le n\)。在长度至多 \(\Delta_0=(Kv)^{-1}\) 的cell中固定左锚点 \(a\)，令 \(\mathbf Q_a=\bigotimes_iQ_{i,a}\)、\(S=\sum_i\xi_{i,a}\)。每个坐标准确满足
\(E_{Q_{i,a}}e^{\xi_{i,L}-\xi_{i,a}}=1\)，所以其锚点相对均值不大于0。原始局部最大增量
\(M=\sup_{L\text{ in cell}}\sum_i(\xi_{i,L}-\xi_{i,a})\ge0\)
受相应中心化过程的绝对最大值控制。

第 \(j\) 层dyadic边的总方差由(3)不超过 \(C2^{-j}\)，因为 \(n\varepsilon^2/v=29/9\)；单项中心化幅度不超过
\(C\varepsilon\min(1,n2^{-j})\)。各坐标独立即可应用Bernstein，不要求同分布。逐层边数至多 \(2^{j+1}\)，取尾参数 \(u+3(j+1)\) 后并集可求和。第30轮 `../030/PROOF.md` 第5节的同一逐边推导遂给
\[
 \mathbf Q_a\{M>C[1+\sqrt u+\varepsilon(1+\log n)^2
                         +\varepsilon(1+\log n)u]\}\le2e^{-u}. \tag{4}
\]
所有无限细层均保留。有限 \(\delta_i>0\) 下路径连续，dyadic极限合法。(1)使尾积分中的常数项有界、线性项系数趋0，所以
\[
 E_{\mathbf Q_a}e^M\le C                                 \tag{5}
\]
对全部来源和cell一致。这里没有使用任何单坐标方差下界、小球或高频—低频交叉界。

由于 \(d\mathbf Q_a=e^S dP^{\otimes n}\)，准确有
\[
 E_P\sup_{L\text{ in cell}}\prod_iY_{i,L}
       =E_{\mathbf Q_a}e^{\sup_L\sum_i\xi_{i,L}-S}
       =E_{\mathbf Q_a}e^M\le C.
\]
不超过 \(1+|J_2|Kv\) 个cell覆盖 \(J_2\)，用 \(\max\le\sum\) 给 \(C(1+Kv)\) 的强积分上界。

## 3. 不同坐标的任意薄宽度使用共同扩长

取 \(\delta_*=(nKv)^{-1}\)，并令 \(\widetilde\delta_i=\max(\delta_i,\delta_*)\)。对应 \(4\le\widetilde w_i\le n\)。对每个坐标的正背景与原子来源，前轮逐原子夹逼给
\(\nu_{i,\delta}(W_L)\le\nu_{i,0}(W_{L+\delta})\)
及 \(\nu_{i,0}(W_L)\le\nu_{i,\delta_*}(W_{L+\delta_*})\)。

若 \(\delta_i\le\delta_*\)，因此
\[
 \nu_{i,\delta_i}(W_L)
 \le\nu_{i,0}(W_{L+\delta_i})
 \le\nu_{i,\delta_*}(W_{L+\delta_i+\delta_*})
 \le\nu_{i,\delta_*}(W_{L+2\delta_*}).
\]
最后一步只用非负测度的区间包含。若 \(\delta_i>\delta_*\)，同一扩长关系直接由正性成立。故所有坐标使用共同长度 \(L+2\delta_*\)，并且
\[
 \prod_iY_{i,\delta_i,L}
 \le(1+2\delta_*/L)^n
             \prod_iY_{i,\widetilde\delta_i,L+2\delta_*}.
\]
共同倍率至多 \(e^{2/(a_0Kv)}\)，\(a_0=10/9\)。由 \(Kv\to\infty\)，最终 \(I+[0,2\delta_*]\subset J_2\)，遂得到(2)。没有为不同坐标选择不同尺度，也没有交换原子极限与上确界。

## 4. 任意概率正混合与随机结构

设 \(F_\omega\) 为上述乘积源的可测族、\(\mu\)为概率律，令 \(F_\mu=\int F_\omega d\mu\)。Tonelli给每周期积分1，从而其密度几乎处处有限。对固定尺度的线性平均与正性，
\[
 \sup_{L\in I}A_LF_\mu
 \le\int\sup_{L\in I}A_LF_\omega\,d\mu(\omega).
\]
正长度紧窗口上积分路径连续，可取可数稠密尺度的sup后使用Tonelli。将一致的(2)积分，即得同一帽。有限混合、连续混合、各坐标不同来源、潜变量引入任意相关的jitter/宽度/系数，均在此范围。

若仅混合同一原来源的等权 \(R\) 个平移，还可结合第31轮的弱混合引理得到 \(C_p\sqrt n\min(1+\log R,\sqrt v)\)。这里的 \(\sqrt n\) 单分量弱界依赖第30轮的原来源条件；不能把它直接套给本定理内任意jitter分量。本定理的 \(C\sqrt{n\log n}\) 强帽自身不需那些标签条件。

这指出随机化研究的明确边界：保持每坐标等质量一格一原子、固定系数包络与 \(K\asymp\sqrt{n/\log n}\) 时，任意上述正混合仍不能达到高于半次对数的固定窗口下界。要取得更高增益，必须离开某个已列假设；要取得较小的正对数增益，当前强帽并未排除，但本轮没有构造或证明。
