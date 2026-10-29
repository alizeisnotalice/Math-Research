# 第39轮：异表计数分类的固定尺度稳定性与后验熵的平方根对数阶

2026-10-29。本轮没有改善已接受的中心立方体√n下界。新的严格成果有两项：真实立方体来源的计数分类在固定宏观尺度下，具有不随维数、分类数增长的L¹稳定性；另一个独立Gaussian辅助模型的最佳后验倍率，其均值准确增长为√log R。后者尚未转成真实输入的弱型增益，不能称为C_n增加了对数。

## 1. 真正立方体核的无表数损失估计

取任意有限族非负一维周期密度ρ_b，周期质量均为1。共同源类别χ∈{0,1}满足∫ρ_bχ=η∈(0,1)，η不随表b改变。将源坐标类别计数S=Σχ(y_i)的取值集合分成任意固定分割E_b，允许空段。记真实π_b=Pr{Bin(n,η)∈E_b}，定义
\[
 F(y)=\sum_b\rho_b^{\otimes n}(y)\mathbf1_{S(y)\in E_b},
 \qquad F_{\rm prior}=\sum_b\pi_b\rho_b^{\otimes n}.
\]
二者周期质量准确为1。对每个固定L≥1，令m=⌊L⌋，本轮证明
\[
 \boxed{\|A_LF-A_LF_{\rm prior}\|_{L^1(\mathbb T^n)}
 \le\sqrt{\frac{\eta(1-\eta)}{m(m+1)}}
 \le\frac1{2\sqrt{m(m+1)}}.}
\]
特别在1≤L<2时可取1/(2√2)。整数L的实际差为0。结论允许各表不同、各表的query倾斜律不同；不需要Sidon或原平方根构造的匹配参数。它是指定来源族、固定尺度的定理，不是一般C_n上界。

证明先用6个Fraction小例核对count posterior的加权Gram恒等式，再用4组异表例及2个共同基准空间例核对不同query律下的推广。关键如下：单坐标binary channel的后验p有均值η，令其signal fraction为f=Var(p)/[η(1−η)]。标准化Bernoulli乘积的对称第k阶正交基，在条件期望下的平方范数准确为fᵏ。因此每个count bin的后验相对方差是f的非负系数多项式。

当第b表各自的signal为f_b≤f_*<1时，不同query律Q_b下仍有
\[
 \sum_b\pi_b E_{Q_b^{\otimes n}}
 \left(\frac{T_b}{\pi_b}-1\right)^2
 \le\sum_{k=1}^n f_*^k\le\frac{f_*}{1-f_*}.
\]
这里T_b是该表第b个bin的真实后验，空段不除π_b=0。跨表的Cauchy给L¹帽√(f_*/(1−f_*))。宏观窗口的完整周期部分使p落在显式区间，结合均值η直接得f_*≤η(1−η)/[m(m+1)+η(1−η)]，即上式。全部推导见 [channel/PROOF](round039_channel/PROOF.md)，[交叉独审](round039_posterior/CHANNEL_REVIEW.md)绑定其最终哈希。

第38轮实际source在L=7/6的四表bin-signal约为.015842–.015866。它低于连续signal，不能直接作连续帽。本轮以Fraction得到连续与bin之间的解析加量≤1.76846993×10⁻⁶；保存的bin值本身仍是普通FP，因此加上解析差额也不构成完整数值认证。严格宏观定理不依赖这些浮点数。

## 2. 可增长的后验熵：辅助量确有√log R阶

固定0<τ<1，s=√(1−τ²)，Y、Z独立标准正态，X=τY+sZ。按X的标准正态先验分成R个等概率分位段，记
\[
 H_R(y)=\max_b\frac{\Pr(X\in I_b\mid Y=y)}{\Pr(X\in I_b)}.
\]
先实际积分R2至4096；随后用单峰性将最大bin搜索缩减为包含x*=y/τ的bin及其两个邻居，算至R=2²⁴。最后完成自含上下界证明：对每个固定τ，全部整数R→∞时
\[
 \boxed{EH_R\sim
 \frac{2\tau}{\sqrt{1-\tau^2}\sqrt\pi}\sqrt{\log R}.}
\]
有限上界为3+2τa_R/(s√(2π))，a_R=Φ⁻¹(1−1/R)。上界来自有限参数区间的密度比包络积分及两个尾段；下界在中央query区间用分位bin宽一致趋零和精确二次曲率匹配。Mills界和quantile阶也在文内推导，不调用未读的渐近定理。见 [Gaussian阶数证明](round039_channel/gaussian_order_PROOF.md)。

τ=1/3时，EH_R在R2、16、256、4096、2²⁴约为1.21635、1.55329、1.83691、2.05445、2.54349。R=2²⁴时EH_R/√log R约.62361，理论极限系数约.39894；有限R仍有明显O(1)偏移，不能据这些小数声称已进入极限。三候选算法与全bin算法逐点最大差为0；半步网格细化的EH最大差约1.66×10⁻⁷。全部是普通FP数值控制，不是区间积分证书。

dyadic R的点态包络则准确为s⁻¹exp(Y²/2)，故任何q≥1都没有R一致的Lq矩，包络本身的weak-L¹仍有限。这说明“先取best posterior，再用一个不存在的一致强积分帽删去它”这条路线不成立。它没有证明这份√log R可以转移给立方体。

也有真实binary product channel的无dimension一致L¹(max)反例：任意非退化固定channel，signal f>0、p的essential range为d，则Chebyshev打包给
\[
 \liminf_n\frac{E\max_s(T_s/\pi_s)}{\sqrt n}
 \ge\frac{3f}{8d}>0.
\]
该channel可由真实周期cube位移构成，但这里最大化包括极小先验count；转成整数分位分类可需指数大的R。没有证明多项式R达到√n。尤其所有表同源时，尽管best posterior增长，gated与prior仍完全相同。

## 3. 真实权重必须保留：两项校准

**有限PB包络。** 另一条线先做12类log-DP探路，再证明自含有限n定理。p_i∈[.3,.7]、η∈[.49,.51]，令Δ=Σp_i−nη、D=Σ(p_i−η)²。D>0时，对任意整数集合E，
\[
 \Pr_p(S\in E)\le\frac{20}{3}
 \exp\!\left(\frac{5\Delta^2}{7D}\right)
 \Pr\{\operatorname{Bin}(n,\eta)\in E\}.
\]
D=0时比例严格为1。常数没有显式R，但D趋零、漂移大或跨表同时好事件仍需收费，不能直接积分这个指数。证明用精确倾斜、Fourier局部质量帽及Binomial整数mode的Chebyshev下帽，覆盖n=1和k=0,n端点，见 [posterior/PROOF](round039_posterior/PROOF.md)。两路复核未发现实质错误；根代理另以909个有理参数配置核对倾斜方差代数，并独立重算精确反例。

被否定的简单候选是PB(E)≤2exp(2d²)π(E)+2/√n。n=256、η=1/2、p_i=9/16、d=1、E=[147,256]给严格有理残差PB(E)−16π(E)−1/8=.0877721969…>0；因2e²<16，原候选失败。保留数值纠错记录：首次线性longdouble DP在n4096的远尾下溢，改log-DP后重算同一批输入，未删失败输入、未追加seed。

**耦合Gaussian位移。** 为检验能否把best posterior熵免费乘到尾上，主线另定义同一个真实概率律：X落入第b分位段时输出Y=t_b+τX+sZ，选择t_b=c−τx_b将各段代表中心送往c。固定表的正确回归均值是τ(y−t_b)，而非τy。随R→∞，gated输出密度在L¹中趋于N(c,1−τ²)，prior趋于N(c,1+τ²)。单尺度weak泛函对L¹收敛连续，故此特定策略在固定c、τ下的weak比值趋于有限常数，不能靠R增大获得无界倍率。

24组确定积分中，τ=1/3、c=8的经验比值从R2约1.0624到R64约1.1151，连续极限网格值约1.1162。dy从.002细到.001时，两个weak系数各改变约0.4%，报告明确没有按展示位数作精度承诺。这不是cube极限定理，也不排除其他位移或c随R变化。完整推导见 [coupled/PROOF](round039_coupled/PROOF.md)。

## 4. 主目标仍缺什么

本轮证明的是固定L的L¹约束，另精确找到了一个辅助后验量的√log R熵。二者均不能替代对同一个真实L¹输入及原Lebesgue查询律的共同尺度尾证明。目前没有c,a>0和合法F_n能使全部充分大n都有W(MF_n)/||F_n||₁≥c√n(log n)ᵃ。

把新定理对J个尺度直接求和，只给误差J/(2√2)：J=233时约82.38。更准确地，对0<ε<1，事件包含及Markov仅给
\[
 W(M_{\mathcal J}F)\le
 \frac{W(M_{\mathcal J}F_{prior})}{1-\epsilon}
 +\frac1\epsilon\sum_{L\in\mathcal J}
 \|A_LF-A_LF_{prior}\|_1.
\]
当J约为Kv时仍保留Kv费用；在v约log n的原预算下，这允许√n√log n，未排除较小的正log增益，也未将其实现。第一处未闭合的研究步骤是把分类熵与真实来源的高尾/局部尺度事件耦合，在共同尺度下保留可累积增益，同时支付prior和查询体积。没有将∫sup与sup∫交换，也没有声称平方根已是双边阶。

## 5. 实际核验、其他对话和复现

两名复用子代理分工推导并交叉审查；主线独审PB包络和Gaussian阶数，并运行有理反例及coupled模型。全部数值均在本机运行，无尾MC、无大维缓存。原agent环境Python3.15.0/NumPy2.5.4/SciPy1.18.1，主线另用Python3.12.14/NumPy2.3.5。Gaussian基础诊断约2.46秒，三候选大R约1.23秒，异表Fraction检查约.36秒，PB log-DP约.14秒，coupled24例及细化约2.77秒。具体运行不视为一般证明的替代品。

只读核对“上界”第643轮正式报告及第644轮正在归档的状态；目前没有可直接迁移的新完整上界结论。另读“验证下界算法成果评估”的已有记录，它同样未提供新的下界定理。未向其他对话发送消息、未恢复暂停任务。

运行 `python3 verify_saved.py` 可在临时副本中重放有理检查、PB扫描及coupled积分；加 `--full` 可重放SciPy通道及大R三候选数值。所有原件只读，输出可用 `--output /绝对路径/新文件.json` 保存。详见 [REPRODUCE](REPRODUCE.md)。本轮精简归档证明、必要源码、最小结果及纠错链；不包含大缓存，历史论文和无关未跟踪文件均未改。
