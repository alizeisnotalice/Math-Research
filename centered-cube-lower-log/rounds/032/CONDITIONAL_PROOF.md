# 按查询区间计数的条件矩下界与有界贡献估计

本接口解决有限尺度数值中的一个计算问题：维数很大时，避免逐一生成全部查询坐标，同时保留到真实Lebesgue查询测度的下界方向。以下恒等式和不等式严格成立；本轮浮点执行没有进行区间舍入认证，故数值文件只是接口验证，不是确定性下界证书。

## 1. 条件计数与真实来源

令 \(Y_{jr}(x)>0\) 是一维周期概率密度 \(\rho_r\) 在第 \(j\) 个确定中心长度 \(L_j\) 的平均，\(1\le j\le J\)、\(1\le r\le R\)。共享来源是 \(F=\sum_r\pi_r\prod_i\rho_r(x_i)\)，\(\pi_r\ge0\)、\(\sum_r\pi_r=1\)。真实共同尺度平均为
\[
 G_j(\mathbf x)=\sum_r\pi_r\prod_{i=1}^nY_{jr}(x_i).
\]
各坐标查询在圆周上独立均匀。分圆周为有限个区间 \(I_q\)，\(p_q=|I_q|>0\)，令 \(N_q\) 为落在各区间内的查询坐标数，则 \(N\sim\operatorname{Multinomial}(n;(p_q))\)。定义准确的一维区间积分
\[
 m_{qjr}=p_q^{-1}\int_{I_q}Y_{jr}(x)dx,
 \quad c_{qjrs}=p_q^{-1}\int_{I_q}Y_{jr}(x)Y_{js}(x)dx.
\]
对称性使共同平均在坐标置换下不变。给定计数后，可固定任意含相应区间数目的坐标排列；各坐标在其区间内独立均匀。因此
\[
 M_j(N):=E(G_j\mid N)=\sum_r\pi_r\prod_qm_{qjr}^{N_q},
\]
\[
 T_j(N):=E(G_j^2\mid N)=\sum_{r,s}\pi_r\pi_s
                            \prod_qc_{qjrs}^{N_q},
 \qquad V_j=T_j-M_j^2\ge0.                                \tag{1}
\]
这里 \(M_j\) 是条件期望，并未被直接当作真实中心立方体平均。

## 2. 用条件方差支付离散化损失

固定阈值 \(\lambda>0\)，设
\[
 B_j(N;\lambda)=
 \begin{cases}
 (M_j-\lambda)^2/[V_j+(M_j-\lambda)^2],&M_j>\lambda,\\
 0,&M_j\le\lambda.
 \end{cases}
\]
若有经认证的方差上界 \(\overline V_j\ge V_j\)，以它替代分母中的 \(V_j\) 仍给下界。Cantelli给
\[
 P(G_j>\lambda\mid N)\ge B_j(N;\lambda).
\]
为展开所用不等式：对中心变量 \(Z\)、方差 \(\sigma^2\)、\(t>0\)，对任意 \(a>0\) 有
\(P(Z\le-t)\le E(Z-a)^2/(t+a)^2=(\sigma^2+a^2)/(t+a)^2\)。取 \(a=\sigma^2/t\) 得 \(\sigma^2/(\sigma^2+t^2)\)；\(\sigma=0\)直接处理。

取仅依赖 \(N\) 的最小序号极大选择 \(j_*(N)\in\arg\max_jB_j(N;\lambda)\)，条件于 \(N\) 后该尺度已经确定，故
\[
 P\{\max_{1\le j\le J}G_j>\lambda\}
 \ge E_{N\sim P_N}\max_jB_j(N;\lambda).                    \tag{2}
\]
真实连续窗口的极大函数还不小于此有限尺度极大函数。因此(2)是有效的下方向接口，无需给连续尺度恢复的上误差。它不声称有限网格已经解析了连续上确界。

## 3. 只用均值与区间振幅的可扩展方差上界

不必形成所有 \(R^2\) 交叉矩。若在 \(I_q\) 上 \(Y_{jr}\) 的振幅至多 \(d_{qjr}\)，Popoviciu不等式给
\(\operatorname{Var}_{I_q}Y_{jr}\le d_{qjr}^2/4\)。该式由 \((Y-a)(b-Y)\ge0\) 先得 \(\operatorname{Var}Y\le(EY-a)(b-EY)\le(b-a)^2/4\)。对分量积 \(H_{jr}=\prod_iY_{jr}(x_i)\)，由区间条件独立性
\[
 \frac{\operatorname{Var}(H_{jr}\mid N)}{[E(H_{jr}\mid N)]^2}
 \le \prod_q\left(1+\frac{d_{qjr}^2}{4m_{qjr}^2}\right)^{N_q}-1
       =:D_{jr}(N).
\]
对共享混合使用条件 \(L^2\) 三角不等式，得
\[
 V_j(N)\le\left[\sum_r\pi_r\left(\prod_qm_{qjr}^{N_q}\right)
                         \sqrt{D_{jr}(N)}\right]^2.        \tag{3}
\]
它保留所有相关性，未把不同相位当成独立变量。若已知区间内导数上界，可用 \(d\le |I_q|\sup|Y'|\)。箱平滑断点处先分片，再对每段处理；不能跨未控制的跳跃套导数。

## 4. 计数层的精确倾斜与阈值无关的贡献上界

给定 \((j,r)\)，单坐标按 \(Y_{jr}\) 倾斜时，区间概率为 \(p_qm_{qjr}\)，总和为1。故计数可以直接按 \(\operatorname{Multinomial}(n;(p_qm_{qjr}))\) 精确抽取。均匀选一个 \(j\)，再按 \(\pi_r\) 选分量，所得计数提议 \(Q_N\) 相对 \(P_N\) 的密度为
\[
 D(N)=\frac1J\sum_jM_j(N)>0.                              \tag{4}
\]
这来自多项式概率中相同组合因子的准确消去，不依赖Gaussian近似，也不需要生成区间内的全部 \(n\) 个坐标。

设 \(B(N;\lambda)=\max_jB_j(N;\lambda)\)。若 \(B>0\)，必有某 \(M_j>\lambda\)，于是 \(D\ge M_j/J>\lambda/J\)。因 \(B\le1\)，得到关键贡献界
\[
 0\le T_\lambda(N):=\frac{\lambda B(N;\lambda)}{D(N)}\le J,
 \quad E_{Q_N}T_\lambda=\lambda E_{P_N}B
 \le\lambda P\{\max_jG_j>\lambda\}.                       \tag{5}
\]
在 \(B=0\) 处贡献定义为0。若以 \(\alpha\) 比例混入均匀计数提议，则上界改为 \(J/(1-\alpha)\)。必须在提议中包含下界选择使用的全部目标尺度；只含三个锚点而对任意新尺度取max，不支持(5)的 \(J\) 上界。

独立抽 \(S\) 个提议计数、对事先确定的 \(\lambda\) 求贡献均值 \(\widehat T\)。由有界变量Hoeffding，概率至少 \(1-\alpha\) 时
\[
 \lambda P\{\max_jG_j>\lambda\}
 \ge\widehat T-J\sqrt{\frac{\log(1/\alpha)}{2S}}.            \tag{6}
\]
双侧区间相应使用 \(\log(2/\alpha)\)。多阈值同时陈述需另付并集费用，或用独立样本选阈值再验证。式(6)的确定数学输入必须是准确或经区间认证的 \(M_j\)、方差上界和提议密度；本轮普通浮点程序没有完成这层舍入认证，不能把其点估计直接登记为严格数值下界。

## 5. 本轮接口验证与剩余问题

本机NumPy实际运行第31轮的正箱峰加第一低频模型，在 \(n=8,32,64\)、\(Q=16,64,256\)、9个共同尺度、5个固定阈值上比较条件矩下界与直接查询，共45组。增加bin数后，本次观察到下界更接近真实有限网格尾；这不是对任意分割的单调性定理。另取一个固定区间排列、65536个区间内查询，解析条件均值和方差与直接计算吻合。

计数倾斜实验使用 \(n=64,Q=64\)、全部9个目标尺度和4相位、8192个提议计数。5个阈值的所有实测贡献均小于1.20，低于理论上界9。程序仍为浮点诊断，不能把观测最大值1.20替代Hoeffding中的确定上界9。

有限 \(J\) 的窗口弱系数本身至多 \(J\)，因为每个固定尺度的平均质量为1。要检测 \(\sqrt n\) 以上的增长，\(J\) 也必须足够大，或补充有效的全窗口概率聚合证明。免去逐坐标成本不等于已经解决所有尺度数、相位数和维数的联合复杂度；本轮没有新对数下界。
