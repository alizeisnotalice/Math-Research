# 低频投影补偿短独审

审查对象仅为 round036_compensation/ 的 PROOF.md、probe.py、positivity_certificate.py 与 positivity_certificate.json。未修改冻结文件，未重复模拟或重放整份证书。

**结论：通过其明确限定的审查范围。未发现源密度、Fourier因子、真实box平均、H/B缩放或正性证书范围上的明确错误。** 这不表示新的平坦Gram或概率链已经闭合；第一项剩余数学缺口仍是归一化后远lag自相关 \(C(L)\) 的统一小帽，PROOF已将它保留。

## 1. Fourier因子与真实平均核对

采用 \(\widehat f(m)=\int f(x)e^{-2\pi imx}dx\)，有
\[
 \widehat{P_\delta'}(m)=z_m,\qquad
 \widehat P_\delta(m)=z_m/(2\pi im).
\]
真实正负频投影的导数为 \(2\Re\sum z_me^{2\pi imx}\)，所以
\[
 Q'=K(\mu_\delta-1)-2\Re\sum z_me^{2\pi imx}.
\]
与 \(D_M=\sum\pi mc_m\cos(2\pi mx)\) 合并后，密度低频的实部系数确为 \(\pi mc_m-2z_m\)。这里没有遗漏或多出一个2。

对中心box平均，
\[
 A_LQ'=\frac{Q(x+L/2)-Q(x-L/2)}L,\qquad
 A_LD_M=B_L/L.
\]
因此 \(Y_L=1+p(H_{Q,L}+B_L)/(KL)\) 正确，来源不随L改选。

probe中原cross项 \(-4c_m\Im(\widehat P_\delta(m))\sin^2(\pi mL)\)，以及移除的HH能量 \(8|\widehat P_\delta(m)|^2\sin^2(\pi mL)\)，符号与实Fourier展开一致。高频Q与低频B严格正交，但投影同时移除了HH能量，不能只从原Gram中删除cross。

## 2. 中心化、目标与H/B归一化

程序先去primitive均值，再用
\[
 C_0=R_\delta(0)-(\int P)^2-2\sum|\widehat P_\delta(m)|^2.
\]
这与Q包含零频投影、\(\int Q=0\) 一致。令 \(\alpha=6C_0\)，则：

- H方案乘 \(s=\alpha^{-1/2}\) 于Q'，实际峰权重变为ps、背景变为1−ps；密度振荡系数为 \(\pi mc_m-2sz_m\)。目标仍为原 \(\gamma L^2/3\)，误差 \(-2C(L)/\alpha-T(L)\) 正确。probe明确保存peak_mass_H_rescaled，未继续冒称峰质量p。
- B方案保留峰p、背景1−p，把D乘 \(\sqrt\alpha\)，目标同时变为 \(\alpha\gamma L^2/3\)，对应参考方差为 \(\alpha v\)。误差 \(-2C(L)-\alpha T(L)\) 正确。probe保存target_scale与effective_v，未把较小绝对误差直接称为额外平坦收益。

准确地 \(G_B=\alpha G_H\)，两方案相对各自目标的相对误差相同。所有来源质量仍为1，因为修正项为周期零积分导数/非零频三角项。

缩放必须以 \(C_0>0\) 为前提；这不是对任意原子/任意宽度自动成立。PROOF后文已明确要求C0远离0，具体九组区间均严格为正，因此这里没有将该条件偷换成一般结论。

## 3. 正性证书：范围适当，未冒称Gram认证

逐项核对：

- 原子以保存的U为精确dyadic，再以有理式 \((j+U_j)/K\) 定义；证明没有把NumPy三角系数当解析精确值。
- z使用箱乘子 \(\sin(\pi m\delta)/(\pi m\delta)\)。b系数的trig参数5m/27、m/18分别对应原 \(\sin(10\pi m/27)\)、\(\sin(\pi m/9)\)，无参数/因子错误。
- R0、均值及平滑能量损失以Fraction计算；短gap的修正计入，未无条件假定箱峰不重叠。
- 基本区间运算使用70位定向舍入；negation使用copy_negate避免默认精度舍入；sqrt向两侧扩一ulp。cos48/sin49多项式配合给定余项帽是充分的。
- 四种来源分别使用与真实峰/背景一致的三角正性下端，特别是H方案背景为1−ps。

本次另实际运行了很小的只读算术检查：证书脚本SHA256与JSON记录一致；36个保存下端全部正，最小为
\[
 0.01680328388484000343929516594345568088\ldots.
\]
Fraction精确核对Machin正切恒等式，以及 \((22/7)^{50}/50!\)、\((22/7)^{51}/51!\) 均小于 \(10^{-35}\)，通过。没有运行新的来源或尾模拟。

证书及PROOF都清楚限定：只认证这些具体解析修正来源的正性/质量；不认证NumPy Gram曲线、连续L补点、目标方差平台、尾概率或任何 \(C_n\) 下界。这个范围没有被扩大。

## 4. 第一未闭合项

投影与常数归一化后仍有 \(-2C(L)\)（H方案另除alpha）。现有正性证书不能控制它，严格正交也不能替代远lag估计。若要将本补偿升级为更宽v窗口的平坦方差/下界方法，首先需要新的、对同一合法来源成立的剩余 \(C(L)\) 统一界，再支付完整概率链。当前文件对此没有作已完成断言。
