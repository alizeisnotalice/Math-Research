# 第18轮：固定比例斜坡的乘法渐近与新的绝对下界项

日期：2026-11-06（Asia/Shanghai）。状态：固定比例指数斜坡族的严格乘法渐近；一般C_n只获得下界，尚无乘法渐近等价。

## 1. 主要结果与范围

令mu=log B、sigma=sqrt(mu)、L=log n。以唯一正根定义
\[
 \overline\Phi_N(z_n)=\sigma z_n/\sqrt n,\qquad d_n=\sigma\sqrt n\,z_n.
\]
Phi_N是标准正态CDF，与输入Orlicz模数Phi_n不同。沿用第17轮软根A_mu(q_n)=1及真实斜坡评分R_n^ramp(mu)=exp(q_n)/I_n。
记
\[
 a=(3-\mu)/6,\qquad b_0=(5\mu-6)/6.
\]
则
\[
\boxed{q_n(\mu)=n\mu-d_n+a z_n^2+b_0+o(1),}\tag{1}
\]
并且
\[
\boxed{R_n^{\rm ramp}(\mu)\sim
 \sqrt{2/\pi}\,e^{b_0}n^{-1/2}
 \exp\{n\mu-d_n+a z_n^2\}.}\tag{2}
\]
这里的乘法等价只适用于一个实际下界族的最优阈值评分，不能把左端替换成C_n、Q_n或自适应上界mathcalR_n。

由于z_n²~L，第17轮实际比例偏移收益给一个更强的一般输入下界
\[
\boxed{\log C_n\ge n\mu-d_n-\frac{\mu+3}{24}\log n+o(\log n).}\tag{3}
\]
这将第14轮下侧的O(log n)损失替换为显式系数。一般上侧仍为nmu−d_n+O(n^(3/8)L^(1/8))，尚未匹配(3)。

## 2. 第三累积量与CDF的一项绝对展开

固定Z=Z_mu=Y+mu J，p=P(J=1)=(mu−1)/mu。独立性直接给
\[
 E Z=\mu,\quad\operatorname{Var}Z=\mu,\quad
 \kappa_3=2+\mu^3p(1-p)(1-2p)=\mu(3-\mu).
\]
设X=(Z−mu)/sigma、gamma=κ3/sigma³，psi为X的特征函数。
本轮使用以下自包含估计：对n步标准化和的CDF H_n，全实轴上一致有
\[
 H_n(x)=\Phi_N(x)+\frac{\gamma}{6\sqrt n}(1-x^2)\phi(x)+O(n^{-1}),\tag{4}
\]
phi是标准正态密度。这里余项为绝对误差，不是相对尾误差。

### 直接Fourier证明

显式连续指数因子给
\[
 |\psi(u)|\le(1+u^2/\mu)^{-1/2}.
\]
所有固定阶矩有限，零附近的解析log特征函数满足
log psi(u)=−u²/2+gamma(iu)³/6+O(u⁴)。取固定充分小epsilon，n步特征函数chi_n(t)=psi(t/sqrt n)^n在|t|≤epsilon sqrt n满足
\[
 \left|\chi_n(t)-e^{-t^2/2}\left(1+\frac{\gamma(it)^3}{6\sqrt n}\right)\right|
 \le\frac Cn(t^4+t^6)e^{-c t^2}.\tag{5}
\]
证明细节：log余项n O((t/sqrt n)^4)=O(t⁴/n)，其真实部分可由C epsilon²t²吸收到负二次指数中；三次项纯虚，利用|e^{iv}−1−iv|≤v²/2给t⁶/n项。
在|t|>epsilon sqrt n，代换u=t/sqrt n并从(1+u²/mu)^(-n/2)中抽出一半幂，得到chi_n/|t|积分为O(e^-cn)；Gaussian修正项的相应尾积分同样指数小。
因此特征函数差除以|t|后的全线积分为O(1/n)，包括t=0附近，因为(5)除以|t|后可积。

近似CDF的导数为phi(x)[1+gamma(x³−3x)/(6sqrt n)]，其特征函数正是(5)的Gaussian修正项；该导数是总质量1的可积有符号密度，不必在所有x上非负。
对两个CDF之差作Fourier反演，绝对值由上述除以|t|的积分控制，即得(4)。也可先对差的密度反演，再以两端0的极限确定积分常数。
这给出所需完整余项，不借用普通Berry–Esseen或只有O(1/L)相对误差的鞍点公式。

## 3. 在隐式正态根附近保留常数项

写z=z_n、beta=d_n/n、f_G=phi(z)/(sigma sqrt n)。Mills关系和根方程给
\[
 z^2\sim L,\qquad f_G\sim z^2/n,\qquad n^{-1}=o(f_G).
\]
对任意固定实数b，测试
\[
 q_b=n\mu-d_n+a z^2+b,\qquad\Delta=a z^2+b.
\]
其标准化位置为−z+Delta/(sigma sqrt n)。把(4)代入并Taylor展开，得到
\[
 F_n(q_b)=\beta+(a+b)f_G+o(f_G).\tag{6}
\]
具体地，Gaussian位移贡献Delta f_G；三阶项贡献a(1−z²)f_G，两者相加只剩(a+b)f_G。
Gaussian二阶位移及三阶项位置变化相对f_G均为O(L^(5/2)/sqrt n)=o(1)；(4)的O(1/n)也为o(f_G)。所以不能删去三阶项中的“−1”：它影响软根的最终常数b_0。

## 4. 软平均的常数项

置d=nmu−q_b=d_n−Delta，U=E exp((S_n−q_b)/n)。单步正log mgf在1/n处的Taylor给
\[
 U=\exp\{d/n+\mu/(2n)+O(n^{-2})\}.
\]
展开并仅保留f_G尺度，得到
\[
 U-1=\beta+(\mu/2-a)f_G+o(f_G).\tag{7}
\]
其中−Delta/n贡献−a z²/n，d²/(2n²)贡献mu z²/(2n)；b/n、mu/(2n)及三次指数项都为o(f_G)。
不能因为最后两者较小就同时删去d²/(2n²)，后者在本尺度决定常数项。

令G=E[e^((S_n−q_b)/n)1_{S_n<q_b}]。负倾斜后一步方差不超过原方差mu：指数部分方差下降，Bernoulli参数≤1/2且随倾斜下降，p(1−p)也下降。
故log E exp(−lambda(Z−mu))≤mu lambda²/2，对所有lambda≥0成立。Chernoff后尾积分给
\[
 0\le F_n(q_b)-G\le\frac1n E(q_b-S_n)_+
 \le\frac\mu d e^{-d^2/(2n\mu)}.
\]
在本测试点d=d_n−O(L)，e^-d²/(2nmu)=O(z²/sqrt n)，所以右边为O(z/n)=o(f_G)。
结合(6)—(7)，精确软平均A_mu(q_b)=U−G满足
\[
 A_\mu(q_b)-1=(\mu/2-2a-b)f_G+o(f_G).\tag{8}
\]
取b=b_0±epsilon，b_0=mu/2−2a=(5mu−6)/6，(8)两侧分别有相反的严格符号。
A对q严格递减，故唯一根落在这两点之间。任意epsilon>0均成立，证明(1)。
此处不需要对第17轮粗软根做未经证明的高精度代换，而是直接在明确测试点上比较符号。

## 5. 模数的精确首因子与实际输入下界

第17轮分母积分可写为
\[
 I_n=\int_0^\infty\exp\{(n-1)[\log(1+u/n)-u/n]\}\,du.
\]
令u=sqrt n v，逐点极限为e^-v²/2。为合法使用支配收敛，对x≥0有
log(1+x)≤x−x²/(2(1+x))。
将v≤sqrt n及v≥sqrt n两区分开，n≥2时积分核分别被e^-v²/8及e^-v/8控制，因而全域可被e^(-min(v²,v)/8)支配。
所以
\[
 I_n/\sqrt n\longrightarrow\int_0^\infty e^{-v^2/2}dv=\sqrt{\pi/2}.\tag{9}
\]
(1)减去log I_n给(2)。第17轮已证明以q_n−1/n取得严格平均，再对每个固定n先截断轮廓、后取物理深层极限，仍实现该上确界评分；损失o(1)不改变乘法渐近。

由(1)、z_n²~L及(9)，固定mu评分满足
log R_n^ramp(mu)=nmu−d_n−mu L/6+o(L)。
再加第17轮明确偏移的(mu−1)L/8+o(L)收益，就得到(3)。
这是真正可积输入给的一般下界，未把受限家族的等价式提升成一般最优常数的等价式。

## 6. 本机软件与复现

MATLAB R2026a复用第17轮实际求得的软根，重新稳定求z_n并比较(1)—(2)，避免重复昂贵混合求和：

|n|q_n−nmu+d_n−a z_n²|固定mu评分/(2)右侧|
|--:|--:|--:|
|256|0.3835055290|1.327215669|
|1024|0.3364924610|1.289227838|
|4096|0.2989255668|1.252470893|
|1e4|0.2787370713|1.231196284|
|1e5|0.2376967849|1.186043705|
|1e6|0.2081907891|1.152889391|
|1e7|0.1863984853|1.128447861|

第一数值列的理论极限b_0=0.06538711897，收敛很慢；第二列极限为1。数据不承担证明，也未给出进入渐近区间的有限n阈值。
Wolfram15符号核验前三累积量与Edgeworth CDF修正的导数，并对n=256、1024用原始Gamma–二项全和独立求固定mu软根，60/100位精度加倍检查；它没有直接使用MATLAB算出的根作为精确输入。
复现：MATLAB `addpath('code'); sm_edgeworth_ramp_checks;`；Wolfram `wolframscript -file code/edgeworth_ramp_checks.wl`。
Wolfram两档精度的固定mu软根吻合至25位，MATLAB Code Analyzer无报告。浮点结果、二项截断和Gamma尾均未区间认证；严格结论来自第2—5节分析证明。算法确定性，无随机种子。

## 7. 剩余差距与下一轮两项动作

一般下界目前log C_n≥nmu−d_n−(mu+3)L/24+o(L)，一般上界仍为nmu−d_n+O(n^(3/8)L^(1/8))，不能据此宣称一般C_n的乘法渐近等价。
下一轮：① 检验反馈Bellman模型在log n尺度的额外收益，尝试构造可验证的离散屏障，避免只用形式扩散方程；② 检验固定斜率指数斜坡以外的截断平台及近阈值幅度优化能否再提高真实下界，并严格结算分母。
所有外部概率来源仍见sources.md；本轮额外需要的Fourier CDF误差及I_n首因子均已直接证明。无新增优先权主张，无冗余跟踪材料需清理。
