# 直接捕获质量带交叉审计：真实有限赢家否定单 pair 高斯衰减

2026-10-07，独立新稿。根 direct_capture_mass_band_square 稿在首次读取时尚未落盘，先按根提供的精确定义工作；终态前已全文核对落盘稿 §1–6，完整 μ、严格幅度门、下开上闭质量带、first moment 与 gate 范围一致。只新增本 prefix，不改原稿/总账，不重跑已有 cosine、spike 或 log-Gram 数值。

主结论：原 hard 幅度带 τ<M≤2τ、原尺度 {a,2a}、完整非负 L1 输入及真实有限 winner 并不足以保证
\[
 \int T_kT_l\,d\mu\le C
   e^{-c(k-l)^2/n}\sqrt{I1_kI1_l}
\tag{1}
\]
中的统一 C<∞、c>0。根提供的余弦机制可以严格完成：k=0,l=n 的左侧归一化比值≥35/384，独立于 n。因此 (1) 的所需 e^{-cn} 衰减不可能。此反例只有两个非空质量带，整体平方仍有常数预算；没有否定一般 aggregate √n·I1 或 actual FIRST 空间付款。

## 1. 原接口及不变的基础预算

对真实选中 Q(x,R(x))，令 m(x)=μ(Q(x,R(x)))、M(x)=m(x)/R(x)^n，E={τ<M≤2τ}。B>0 预固定，
\[
 D_k=\{x\in E:2^kB<m(x)\le2^{k+1}B\},\quad
 T_k(y)=\tau\int_{D_k}\frac{1_{Q(x,R(x))}(y)}{m(x)}dx,
 \quad I1_k=\tau|D_k|.
\tag{2}
\]
保留完整 μ、实际 winner 和全部 source endpoints；不按带重启来源。Tonelli 给 ∫T_kdμ=I1_k。捕获 y 的 receiver 且位于 D_k 时有
R^n=m/M≤2^{k+1}B/τ，故 receiver 在以 y 为中心、side (2^{k+1}B/τ)^{1/n} 的 cube 内；结合 m>2^kB 给 T_k≤2。因此
\[
 \sum_k\int T_k^2d\mu\le2I1,\quad
 \int T_kT_ld\mu\le2\min(I1_k,I1_l).
\tag{3}
\]
式 (3) 使用来源一次，但不含带差衰减。由于 R∈[a,2a]，总可用 massband 数约 n，不可从逐带 diagonal 免费得到 √n。

## 2. 完整原 L1 来源、硬幅度门、共同有限尺度

固定任意整数 n≥2，
\[
 \epsilon=1/8,\quad \tau=B=3/4,\quad L=8n,\quad
 f(y)=1_{[-L,L]^n}(y)\,[1+\epsilon\cos(\pi y_1)],
 \qquad \mu=f\,dy.
\tag{4}
\]
f≥7/8 于来源盒，且 f≤9/8；完整源质量包括近均匀背景。L 为偶整数，cos 的完整积分为0，因此
\[
 W=(2L)^n.
\tag{5}
\]
原尺度族是真正共同的有限集合 {1,2}，a=1；h_R 为原 n 维中心 cube 平均，真实 winner 由这同一输入的两个响应比较决定，按任意原 tie 规则保留零集。令 M=max(h_1f,h_2f)、E={M>τ}。全空间 M≤9/8<3/2=2τ，故完整 E 已是原幅度带，不能使用 τ=1/2（那会让正余弦部分违反上幅度带）。

若 x 属于 receiver core C=[-L+1,L-1]^n，两 query 都完整在来源盒，准确有
\[
 h_1*f(x)=1+\frac{2\epsilon}{\pi}\cos(\pi x_1),
 \qquad h_2*f(x)=1.
\tag{6}
\]
第一式只平均第一坐标，其他坐标背景不变；第二式第一坐标正好整周期，余弦平均0。前系数为严格正数，且 2ε/π≤ε（也可直接由平均 |cos|≤1 得幅度≤ε）。所以 core 上：

- cos(πx1)>0 时唯一 winner 为1，m=M∈[1,9/8]⊂(B,2B]，即 D_0；
- cos(πx1)<0 时唯一 winner 为2，M=1、m=2^n∈(2^nB,2^{n+1}B]，即 D_n。

零 cos 的超平面 Lebesgue 测度为0；没有为正体积输出选择 tie。core 全部位于完整 E，两个部分都正体积。

更强地，整个 E 上只有这两个质量带：若真实 R=1，m=M∈(3/4,9/8]；若 R=2，m=2^nM∈((3/4)2^n,(9/8)2^n]。因此 boundary、全部额外输出仍自动落 D0 或 Dn，没有被忽略的其他 massband。

完整 E 包含 C 且包含在 [-L-1,L+1]^n，因为 source 盒外距离∞>1 的 receiver 不可能捕获来源。于是
\[
 I1_k\le\tau(2L+2)^n \quad(k=0,n).
\tag{7}
\]
没有用周期 torus 代替真实 Lebesgue 空间，也没有把背景 W 免费化；boundary 真输入/真 winner 完整保留。

## 3. 同一 source 上的交叉下界

取内层 source 盒 F=[-L+2,L-2]^n。对 y∈F，所有接收 x∈Q_1(y) 或 Q_2(y) 都在 C 内。令
\[
 A(y_1)=\int_{y_1-1/2}^{y_1+1/2}
                      1_{\{\cos(\pi x_1)>0\}}dx_1.
\]
按真实正/负 winner 的 receiver 分区给
\[
 T_0(y)\ge\frac{\tau}{1+\epsilon}A(y_1),
 \qquad T_n(y)\ge\frac{\tau}{2}.
\tag{8}
\]
第二式中 Q2 的第一坐标长度2，恰一个整周期，负 cos 区长度为1；其它 n−1 坐标长度2。因此原 kernel 积分为
τ·2^{n-1}/2^n=τ/2。第一式使用 m1≤1+ε、其它坐标 side1，不删 source 权。

A 是2周期三角函数：在 [0,1] 等于1-y1，在 [1,2] 等于 y1-1；一个整周期积分为1，平均1/2。F 的第一坐标长度2L-4是偶整数，所以精确
\[
 \int_F A(y_1)dy=|F|/2.
\tag{9}
\]
又 μ(dy)≥(1-ε)dy 于 F，故
\[
\begin{aligned}
 \int T_0T_n\,d\mu
 &\ge\frac{(1-\epsilon)\tau^2}{4(1+\epsilon)}|F|\\
 &=\frac7{64}(2L-4)^n
 \ge\frac{21}{256}W.
\end{aligned}
\tag{10}
\]
最后一步是 (1-2/L)^n=(1-1/(4n))^n≥3/4，Bernoulli 不等式保证。所有 source profile 用同一完整 μ、同一真实 finite winner。

结合 (7)，得到
\[
\begin{aligned}
 \frac{\int T_0T_n\,d\mu}{\sqrt{I1_0I1_n}}
 &\ge\frac{7}{48}
        \left(\frac{L-2}{L+1}\right)^n\\
 &=\frac7{48}\left(1-\frac3{8n+1}\right)^n
 \ge\frac{35}{384}.
\end{aligned}
\tag{11}
\]
所以 (1) 对 |k-l|=n 要求 35/384≤C e^{-cn}，任何固定 c>0、C<∞ 最终都失败。上下幅度门、共同 B、两个完整来源 profile、严格 winner 和真正 L1 全域都已核验；不是对任意 selector 的抽象 incidence 反例。

## 4. 精确失败原因与仍有效的 aggregate 合同

质量差这里主要来自相同近均匀源在两个物理尺度的体积 2^n，而非来源群体的指数分离。正/负余弦 receiver 分区沿第一坐标互相交替。将各自输出通过原 h1、h2 拉回同一原 source 后，两者大范围同时正且常数阶。source-overlap 不会随着质量差 n 自动消失。

此输入全 E 仅 D0,Dn，(3) 给 S=T0+Tn≤4，所以
\[
 \int S^2d\mu\le4\int Sd\mu=4I1.
\tag{12}
\]
因此 (11) 不否定 aggregate ∫(ΣT_k)²dμ≤C√n I1；它只排除使用质量带索引差的统一 Gaussian/offdiagonal Schur 证书。

一般正确的弱式是
\[
 \int(\sum_kT_k)^2d\mu
 \le2\sum_{k,l}\min(I1_k,I1_l)
 =2\int_0^\infty N(t)^2dt,\quad
 N(t)=\#\{k:I1_k>t\}.
\tag{13}
\]
它准确保留 receiver mass 的分布，没有按來源重收 W，但最坏仍付非空带数约 n；没有证明 N 的加权平方只能 √n。带差 Gaussian 被否定后，需要 actual overlap 的整体结构或 weighted block 预算，而非把普遍 pair 衰减当免费几何性质。本文没有把 (13) 的未付强度冒充一般进展。

## 5. 其他源族及下界压力的查重范围

已读原 spike+background/scale-aware/真实 log-Gram 稿。单尖峰在广半径 shell 中，各 band 体积可指数不同，公共尖峰交叉可能被 √I1kI1l 归一化压小；这不能证明 (1) 对全部输入成立。tiny radial helper 的旧 proxy 放大亦不能用作本题的直接 T_k 反例，必须回到全 μ 的 exact m(x)。

余弦例使用完整近均匀背景，W 与 |E| 同阶，避免旧单峰背景指数质量压倒占用的漏洞；无需再为容易的多峰特殊族重复扫描。它已经否定 (1) 的一般量词。下界总表 A/B 已全文查读：RS/Sidon、非均匀正联合权重、增长 Cantor/XOR 等族都不允许暗限 source 标签数或误设 independence。本例只用单原正 Fourier 模式，且是过强 pair 合同的反例，不是一般弱下界改进。未运行或重跑下界模型；一条明确解析 counter 已足够终止强合同探路。

原 actual FIRST/fullfuture/LCA/history 门未在本输入上认证，可能被其它旧 paid 资格覆盖。因此本稿不作为 R_angle 反例，也不新增原主账费用。

## 6. 新有理证书登记

三轮 n=4,16,64，L=8n，ε1/8、τ=B3/4。解析 Fourier 恒等式 (6) 和 sign partition 周期自含；不需要 π 浮点值。新 Fraction 守卫将检查真正 h1/h2 的 sign receiver 占用长度（piecewise interval intersection）、A 的完整周期面积、全 source/core/query 包含、原幅度与 massband 残差、(10)/(11) 的几何常数及 (12) 两带 aggregate 范围。

无随机种子，无 Monte Carlo，无周期 torus 替代、无 n 维 receiver 网格，旧模型不重跑。数值是解析原核组件与常数守卫，不声称 sampling actual FIRST。

## 7. 终态收据

新守卫注册后一次完成，99 项 Fraction 检查 PASS：三个维数各23项，公共原 sign-interval/三角面积/注册组件30项；耗时约0.00037秒，无 live handle。数值检查不评价 π 浮点值，原余弦 Fourier 恒等式由 (6) 解析给出。

三个 n=4,16,64 的精确归一化交叉下界均≥35/384，输入质量为真实 W=(16n)^n。结果保存各轮有理实际常数，不拟合 c 或 C，也不借 finite probes 宣称 asymptotic：否定 (1) 的所有 n 推导由 (11) 自含证明。

- 脚本：[direct_mass_band_cross_audit_20261007_exact_guard.py](direct_mass_band_cross_audit_20261007_exact_guard.py)，SHA-256：d1f74900011f0af828880221da434f29c13f565798f05fd6106ffca03a81e62c。
- 预登记：[direct_mass_band_cross_audit_20261007_registration.json](direct_mass_band_cross_audit_20261007_registration.json)，SHA-256：929a95a1af827895831939de0cfd5169de698af0dddb52aa7babf62ccc7b2eae。
- 终态：[direct_mass_band_cross_audit_20261007_results.json](direct_mass_band_cross_audit_20261007_results.json)，SHA-256：9b398a90e8201bdd1dcd76ab9de465129651de6a978e87b17888d7dc53dfd100。

一般 aggregate √n 交叉预算仍需研究；本次可靠产出是一个保留全部真实 hardwinner/source/window/amplitude 合同的单 pair 反例和准确较弱 (13)，没有原 actual 付款或弱下界改进。
