# 第39轮：分类通道的熵与固定尺度约束

本轮得到一个新的严格接口：共同class质量η、任意count分割、各表来源可不同且query倾斜可不同，固定宏观L≥1时，真实gated cube平均与实际prior平均满足

\[
\|A_LF-A_LF_{prior}\|_1\le
\sqrt{\eta(1-\eta)\over\lfloor L\rfloor(\lfloor L\rfloor+1)}.
\]

L∈[1,2)的universal帽为1/(2√2)≈.353553，无n或R损失。该结论不控制跨L的sup费用，不给weak maximal排序或Cn新阶。另一方面，best posterior likelihood确有R增长的熵，不能用并不存在的一致Lq(max)帽删去这一机会。完整自含推导见 `PROOF.md`。

## 先运行的计算

`probe.py` 在写证明前实际运行2.462秒，Python3.15、NumPy2.5.4及SciPy normal-CDF/quantile。没有尾MC。Gaussian通道τ=.2、1/3、.6、.8，R2至4096；Y在[-12,12]以.005步长做确定积分，普通FP未认证。代表结果：

| τ=1/3，R | E H_R | E H_R² |
|---:|---:|---:|
|2|1.21634643|1.50357475|
|16|1.55328700|2.78603954|
|256|1.83690528|5.03821437|
|4096|2.05444675|8.62376203|

随后证明：对任何0<τ<1，dyadic R上 `sup H_R(y)=exp(y²/2)/sqrt(1−τ²)`。因此所有q≥1的R一致矩都发散；0<q<1的fractional矩有精确有限帽。包络weak-L1仍为1/sqrt(1−τ²)。这是可用的后验熵来源，不是“gated无收益”的证明，也不是cube极限定理。

独立二元product channel p∈{1/3,2/3}用真实有限PB卷积计算，全部count singleton分类的结果如下。这里后验query law是严格定义的binary channel，不是Gaussian代理。

| n | E max_s(T_s/π_s) | 加权L2 trace | 描述性weak-L1(max) |
|---:|---:|---:|---:|
|31|2.57681963|1.12500000|1.07959207|
|63|3.24281960|1.12500000|1.07005728|
|127|4.18189511|1.12500000|1.06534033|
|255|5.50663268|1.12500000|1.06299552|

n255的max二阶矩约9.315×10^11，加权trace仍1.125，展示unweighted max与prior加权平方帽的区别。有限整数quantile R=1024的E max约1.94–1.97，并未数值显示随n的√n增长；√n证明使用全部count，或足够大的、可指数增长的R。

## 核心正接口

先实际运行 `trace_checks.py`：6个Fraction例（η1/2或2/5，n1/3/5）验证精确

`Σ_s π_s E_Q L_s² = Σ_{k=0}^n f^k`，

其中f=Var_Q(p)/[η(1−η)]。自含证明用Bernoulli-count的正交基：条件后的第k阶基像互相正交，其平方范数为f^k。由此每个count bin的diagonal variance是f的非负系数多项式。

父线提出多表扩展后，先运行 `multitable_checks.py`（.36秒）：η=.5/.4、n3/5、R4，四表不同query law、不同p分布、不同f。两个相同f但不同分布的完整count Gram逐项完全相等；选择各表自己的第b个bin后，prior加权diagonal仍被共同f*的coarsening trace控制。n3在共同四点base上取不同A_b，直接枚举的gated/prior L1也满足帽。空段保留，π_b=0不作除法。

证明进一步给出

`Σ_b π_b E_{Q_b^n}(T_b/π_b−1)² ≤ f*/(1−f*)`，

`∫|Σ_b A_b(T_b−π_b)| ≤ sqrt(f*/(1−f*))`。

共同η与共同count分割是必要接口；各表的Q_b不需要相同。宏观窗口的p范围及均值η给 `f*≤η(1−η)/[m(m+1)+η(1−η)]`，m=floor L，从而得到开头的真cube结论。没有将不同表的最大值当成同一个可选来源，也没有把Gaussian相关性直接塞入cube。

## 真实来源signal诊断

只读38轮 `Source.evaluate(7/6,b)` 四次，不建R/J大cache。K255/v8/R4/w4/M3、Q48960、实际η≈.498047800889。按query bin的真实倾斜质量 `(mean_plus+mean_minus)/Q` 计算pbar。四表bin-channel signal为

| table | f_bin（普通FP） |
|---:|---:|
|0|0.0158626020113|
|1|0.0158632267005|
|2|0.0158418121256|
|3|0.0158656379883|

normalizer观察为1，E_tilt pbar与实际η最大差在浮点末位。bin p范围约.426至.571。这是bin-channel signal；由全方差，连续query f不小于它，不能直接拿它当连续上帽。

`gap_checks.py/results`另以Fraction核验父线给出的连续-bin接口：class range≤B/L、total kernel≥1/L时，`0≤f_cont−f_bin≤B²/[η(1−η)]`。本参数B=29879/44945280，严格η∈[.49,.51]，增量上界1.76846993×10^-6。这个解析差额不认证保存的FP bin f本身；三角函数、全部kernel运算与后续矩统计尚未完整向外认证。无需这些数值，universal宏观L1结论已成立。

## 后验max的精确反例量词

对任意固定非退化binary product channel，f>0、d=essrange p，Chebyshev打包给

`E max_s L_s ≥ (3/4)(floor(n/D)+1)`，`D=floor(2d sqrt(n)/f)+1`，

故liminf Emax/√n≥3f/(8d)>0。它来自相隔O(√n)的条件count均值和互不相交的query统计量窗口。对η=.5可取整数quantile R=2^n获得所有非空singleton；一般R>min(η,1−η)^−n足够。未证明多项式R即可达到√n。

这一channel可由一个固定真实cube来源的“源点+均匀cube位移”Bayes coupling实现；正密度、单弧类别、非整数L保证f>0。然而同源所有表时gated与prior恰好相同，即使best posterior增长。改变各表以利用这份熵，需要同时支付真实source/query倾斜及共同L的几何。新多表L1界正是其中一个约束，没有排除跨尺度新结构。

## 只读互审及残留任务

父线 `round039_coupled/PROOF.md` 已审：联合律、回归τ(y−t_b)、gated项不再乘π、分位中点步函数L1、prior N(c,1+τ²)极限及weak连续性均通过；未重跑其确定积分。它仅控制固定c、τ下这一特定位移策略，无cube推广。

研究代理的 `round039_posterior/PROOF.md` 也已只读互审通过：7/20的log二次帽、49/119的odds比较、倾斜variance≥(2/5)V、Fourier质量帽与Bin mode/Chebyshev导出的20/3、k0/n及n1端点、D恒等式与两个反例量词均正确。该式是保留variance-defect退化的点态posterior/prior包络，未冒充可积或weak控制。

研究代理反向独审本稿§2–5、§7后记录通过，详见 `../round039_posterior/CHANNEL_REVIEW.md`，审查的PROOF SHA256为 `6fefe580ec45be62cc5e0b1d00623c50cabe35f8ba0a58fce71f5d0ae78b5ba6`，与当前文件一致。固定尺度trace/L1使用的是**完整、未另行条件于query safe事件的单位质量product倾斜**；条件或截断query后若失去E p=η，不可重新套用同一signal/basis公式。

本轮未改变他人文件、未Git、未开代理、未追加尾MC。真实共同尺度max的正对数增益仍是研究目标；上述固定L帽不允许交换∫与sup_L，也不免除J~Kv随n增长后的费用。
