# upper_singular057：正推前降维、完整 Green 的 power/log 包络与交叠可积性

本批对五维、16 条未单位化 body 方向、共同 Cauchy 时间的模型证明新的全时间上包络。记
\[
V=\{(1,\epsilon_2,\ldots,\epsilon_5):\epsilon_j\in\{-1,1\}\},
\quad\psi(\xi)=\sum_{v\in V}|v\cdot\xi|,
\quad p_t=\mathcal F^{-1}e^{-t\psi},\quad G=\int_0^\infty p_t\,dt.
\]
Fourier 逆变换含 \((2\pi)^{-5}\)。若 \(x\notin S_2:=\bigcup_{v\ne w}\operatorname{span}(v,w)\)，则
\[
\boxed{
G(x)\le {16\over\pi\sqrt5\,|x|^2}
\sum_{v\in V}{1\over r_v^2}
\left[A_0+A_1\sum_{w\ne v}{1\over|P_{v^\perp}w|}
\log{2r_v\over\rho_{vw}}\right],}
\tag{1}
\]
其中 \(r_v=\operatorname{dist}(x,\mathbb Rv)\)、\(\rho_{vw}=\operatorname{dist}(x,\operatorname{span}(v,w))\) 都是全空间 Euclidean 真距离。严格常数为
\[
A_0=288P_4+{8P_3S\over\pi}< {27\over10000},\qquad
A_1={16P_3\over\pi}< {9\over25000},\qquad
{16\over\pi\sqrt5}< {23\over10},
\tag{2}
\]
\[
P_2={1\over256\pi},\quad P_3={1\over1024\sqrt2\pi^2},
\quad P_4={3\over65536\pi^2},\qquad
S={5\sqrt5\over4}+5\sqrt{5/6}<{15\over2}.
\]
在奇异点按右端为无穷的扩展值解释，不以零乘无穷的未定义形式代入。

该式给出完整正单位面上的新范围
\[
\boxed{G(1,\cdot)\in L^p([0,1]^4)\quad(1\le p<2).}\tag{3}
\]
更精确地，任一 S1 点附近是可积的 \(r^{-2}\) 径向幂乘角向对数；远离 S1 时仅有有限和的 S2 对数。多个 pair-plane 在非 S1 点相交不会产生对数的乘积。本文不声称 \(p=2\) 可积，也不把 (3) 转成未经计算的数值范数。

本次恢复没有重新运行已有几何或 Arb 常数。已有 `results.json` 和 `check.py` 保留原字节；另作只读证据审计。没有完整管积分的数值证书，没有新的 I、全域 g 或 C5 上界，正式 \(C_{\rm cert}\le C_5\le32\) 不变。后续管预算必须另行登记。

## 1. 通用正测度降维：不取 Cauchy 二阶矩

设 E 是任意非零 Euclidean 子空间，P 是正交投影，\(a_v=Pv\)。令独立 \(Z_v\) 的密度为
\[
s_t(u)={t\over\pi(t^2+u^2)},\qquad X_E=\sum_va_vZ_v,
\]
其 E 上的密度记为 \(q_t^E\)。零向量 \(a_v=0\) 的变量先积分掉，不能在后面写 \(|a_v|^{-1}\)。符号方向的对称性给
\[
\sum_{v\in V}vv^T=16I_5,\qquad
\sum_{a_v\ne0}a_va_v^T=16I_E.
\]
因此线性映射 \(Au=\sum_va_vu_v\) 的算子范数为 4，对每个有限 u 有
\[
|Au|^2\le16\sum_{a_v\ne0}u_v^2.\tag{4}
\]

取非负紧支撑连续测试函数 f。把 (4) 乘以 \(f(Au)\prod_vs_t(u_v)\) 后，在非负积分中使用
\(u_v^2s_t(u_v)\le t/\pi\)。对固定非零 \(a_v\)，所替换的变量变成 Lebesgue 测度；其余变量仍为概率测度。沿单位向量 \(e=a_v/|a_v|\) 的坐标变换 \(r=|a_v|u_v\) 的 Jacobian 是 \(|a_v|^{-1}\)。正交于 e 的剩余推前恰是 \(q_t^{E\cap v^\perp}\)，故
\[
\int_E f(z)|z|^2q_t^E(z)\,dz
\le {16t\over\pi}\sum_{a_v\ne0}{1\over|a_v|}
\int_E f(z)q_t^{E\cap v^\perp}(P_{E\cap v^\perp}z)\,dz.
\]
Tonelli 适用于这些正测度。右边局部有限，且下节 Fourier 界给每个正时间密度连续，故测试函数不等式提升为处处成立的
\[
\boxed{|z|^2q_t^E(z)\le {16t\over\pi}
\sum_{Pv\ne0}|Pv|^{-1}
q_t^{E\cap v^\perp}(P_{E\cap v^\perp}z).}\tag{5}
\]
这里没有断言 \(\mathbb E Z_v^2<\infty\)，没有计算发散矩之差，也没有将带权推前错误归一化为概率。第三次使用 (5) 时只删去真正投影为零的方向。

## 2. 统一商峰值与第三步的几何常数

017 的已证单位球 \(\{\psi\le1\}\) 是 280 个真顶点的凸包。保存表的顶点平方范数分布为
\[
1/128:40,\qquad1/144:80,\qquad1/162:160.
\]
凸包中任一点的范数不超过最大顶点范数，所以
\(\psi(\xi)\ge8\sqrt2|\xi|\)。任意 E 上的商指数就是 \(\psi|_E\)，没有因投影而另加尺度。若 \(d=\dim E\)，Fourier 绝对值及径向积分给
\[
0\le q_t^E(z)\le(2\pi)^{-d}\int_E e^{-8\sqrt2t|\xi|}\,d\xi
=P_dt^{-d},\quad
P_d={(2\pi)^{-d}|S^{d-1}|\Gamma(d)\over(8\sqrt2)^d}.
\tag{6}
\]
d=2,3,4 即为 (2) 的 P 常数。

固定 v。其余十五个 \(|P_{v^\perp}w|^2\) 的分布是 \(16/5\) 五个、\(24/5\) 十个，所以倒数长度之和恰为 S。首坐标均为 1 的 body 向量中，没有三个不同向量线性相关：否则三个不同四维 cube 顶点落在同一仿射直线上；一条直线一旦改变某个 ±1 坐标，就只能在该坐标的两个允许值处命中 cube 顶点。

因此对不同 v,w，第三次降维恰有十四个非零方向。三列 ±1 向量的任一非零 3×3 子式都被 4 整除，故至少有一个子式的平方不小于 16。Cauchy–Binet 给三向量 Gram 行列式至少为 16；两向量 Gram 行列式 \(25-(v\cdot w)^2\le24\)。Gram–Schmidt 的体积比于是给
\[
|P_{\operatorname{span}(v,w)^\perp}u|^2
={\det\operatorname{Gram}(v,w,u)\over\det\operatorname{Gram}(v,w)}
\ge {2\over3}.\tag{7}
\]
这是统一的解析下界，不需要枚举旧 560 个 triple。

## 3. 三次降维、两个时间端点与对数常数

第一次在 E=R5 使用 (5)，\(|v|=\sqrt5\)，得
\[
G(x)\le {16\over\pi\sqrt5|x|^2}\sum_v H_v(x),\qquad
H_v(x):=\int_0^\infty t q_t^{v^\perp}(P_{v^\perp}x)\,dt.\tag{8}
\]
所有项非负，交换时间积分不要求事先知道 G 有限。以下固定 \(r=r_v>0\)、\(\rho=\rho_{vw}>0\)。在三维商 \(E=\operatorname{span}(v,w)^\perp\)，(5)–(7) 给
\[
q_t^E(P_Ex)\le\min\left\{{P_3\over t^3},\ {B\over\rho^2t}\right\},
\quad B={16\over\pi}\,14\sqrt{3/2}\,P_2
={7\sqrt{3/2}\over8\pi^2}.\tag{9}
\]
这里第二项是第三次降维后二维峰值的结果，未丢弃第三步；单独使用三维峰值会在时间零点产生不可积的 \(t^{-1}\)。

取预登记 \(T=r/24\)。在 \([T,\infty)\)，四维峰值直接给
\[
\int_T^\infty tq_t^{v^\perp}\,dt\le {P_4\over2T^2}={288P_4\over r^2}.\tag{10}
\]
在 \((0,T)\)，第二次使用 (5)，再代入 (9)。令
\(\tau=\rho\sqrt{P_3/B}\)，则所需的标量积分准确为
\[
J(T,\rho)=\int_0^T\min\{P_3/t,Bt/\rho^2\}\,dt
=\begin{cases}
BT^2/(2\rho^2),&T\le\tau,\\
P_3[1/2+\log(T/\tau)],&T>\tau.
\end{cases}\tag{11}
\]
小时间端点由 \(Bt/\rho^2\) 可积，大时间端点已由 (10) 覆盖。没有把某一个分支延伸穿过它不适用的交点。

因为 \(0<\rho\le r\)，且
\[
{B\over576P_3}={14\sqrt3\over9}<4,
\]
两种情形统一给
\[
J(T,\rho)\le P_3\left[\frac12+
\log_+\left(\sqrt{B/(576P_3)}\,{r\over\rho}\right)\right]
\le P_3\left[\frac12+\log{2r\over\rho}\right].\tag{12}
\]
把 (12) 乘 \(16/(\pi r^2)\)，按十五个 w 求和，再加 (10)，得到 (1) 与 A0、A1，所有系数与原 `check.py` 一致。

## 4. S1、S2 的实际无穷值与局部上下阶

记 \(S_1=\bigcup_v\mathbb Rv\)。非零 S1 点可用一个非零系数表示，S2\S1 点可用两个非零系数表示。移去一或两方向后，其余方向仍张成 R5，因为
\[
\sum_{u\notin\{v,w\}}uu^T\succeq(16-10)I_5=6I_5.
\]
对应余项密度满足稳定缩放 \(h_t(z)=t^{-5}h_1(z/t)\)，且 \(h_1\) 连续、\(h_1(0)>0\)。后者由可积的严格正 Fourier 被积函数在零点的积分得到；连续性还给原点小球上的正下界。

在固定最小表示 \(x_0=\sum_{j=1}^ma_jv_j\)、\(m=1,2\)、\(a_j\ne0\) 上，把每个被移走变量限制在 \(a_j\) 附近长度为常数乘 t 的区间。其 Cauchy 密度在那里至少为常数乘 t，余项密度至少为常数乘 \(t^{-5}\)。因而当 \(|x-x_0|\le c t\)、t 足够小时，
\[
p_t(x)\ge c' t^{2m-5}.\tag{13}
\]
这给 S1 上 \(t^{-3}\)、S2\S1 上 \(t^{-1}\) 的不可积下界，所以 canonical G 在全部 S1/S2 上保留实际无穷值。对附近点积分 \(t\gtrsim|x-x_0|\) 还分别给 \(|x-x_0|^{-2}\) 与 \(\log(1/|x-x_0|)\) 的下增长率。常数只在固定表示及其局部范围内使用。

在远离 S1 的紧集，全部 r_v 一致正，(1) 给有限和的 pair-plane 距离对数。靠近某个 S1 点时，选相应 v，以 \(z=P_{v^\perp}x=r\omega\) 记
\[
\rho_{vw}=r\,\operatorname{dist}(\omega,\mathbb R P_{v^\perp}w),
\qquad\log(2r/\rho_{vw})=\log(2/\sin\theta_w).\tag{14}
\]
所以不能写成在整个 S1 cap 内统一的纯 \(C r^{-2}\) 上界：角向接近 S2 时对数无界。由于每个 p_t 连续且非负，Fatou 给 G 下半连续；在任意靠近 S1、但不同于 S1 的 S2 点处 G=∞，从 S2 外逼近也必无界。任何仅依赖到 S1 距离、并在该非 S1 点附近有限的纯幂上界都因此失败。删除 S2 本身这个零测集合也不能修复该失败。

## 5. 完整面上的可积性与所有交叠

仿射面 H={x1=1} 与每条 Rv 只交于 v。令 h=x−v∈H0={h1=0}。从 \(|v_{2:5}|=2\) 得
\[
|h|/\sqrt5\le |P_{v^\perp}h|\le|h|.
\]
映射 H0→v⊥ 的四维体积 Jacobian 为 \(1/\sqrt5\)，所以面测度等于 \(\sqrt5\,dz\)。S1 项的 p 次幂局部积分因而由
\[
C\int_0^\varepsilon r^{3-2p}\,dr
\int_{S^3}\left[1+\sum_{w\ne v}\log(2/\sin\theta_w)\right]^p d\sigma(\omega)
\tag{15}
\]
控制。每条线在 S3 上只有两个点；点附近表面积有 \(\sin^2\theta\,d\theta\) 因子，故任意有限 p 的 \(\log^p(2/\sin\theta)\) 可积。有限和的 p 次幂可用凸性按有限项拆开；没有要求不同方向的对数奇性互不相交。径向积分则恰在 p<2 时有限。

在不含相应 v 的有界面区域，r_v 一致正。每个 pair-plane 与 H 相交为一条仿射线，限制到 H 的真距离与到该仿射线的距离相当；三维横截面上的 \(\int_0^\varepsilon s^2|\log s|^pds\) 有限。在多个 pair-plane 的相交点，同样拆分有限和的 p 次幂即可覆盖全部交叠。因此 (15) 与这些区域组成完整有限覆盖，证明 (3)，尤其 q=5/4 的完整 I 可积。这是解析证明，不是空间采样。

有限首批几何与此一致：完整 H 上有 16 个 S1 点，正面仅 v+=(1,1,1,1,1)。正面 S2 为十五条
\[
L_D:\quad x_j=t\ (j\in D),\quad x_j=1\ (j\notin D),\quad0\le t\le1,
\quad\varnothing\ne D\subseteq\{2,3,4,5\}.
\]
这十五条线的 105 对互异线仅在 t=1 的共同 S1 角点相交。但 t=0 的端点还被 \(2^{|D|-1}\) 个不同 pair-plane 命中：D 外两向量均为 +1，D 内两向量符号相反，取无序对即得该数；这些端点自身不在 S1。若在端点表示 \(x=av+(1-a)w\)，D 内的零强制 a=1/2 与相反符号，D 外的一强制同为 +1，这也证明没有遗漏其它命中平面。全部端点 incidence 共 40 项。由有限和可积性，它们全部已被支付。

正面十五线的完备性消费 034 的 120 pair-plane 分类；本批没有重枚举该分类。原首批只核验十五线、105 线对及十五端点命中式，不以这些新核验代替旧完整分类。

## 6. 下一批管预算的确定接口

034 的 \(\ell^\infty\) 管并集记为 Tδ，\(0<\delta\le1/3\)，其完整体积为
\(V(\delta)=108\delta^3-243\delta^4\)。不能只用 (3) 的有限性宣称某一数值 δ 已通过预算。一个可完全正积分实现的接口如下；本批只推导，不数值评估它。

令 \(a=27/10000\)、\(b=9/25000\)、\(s=15/2\)、\(f=23/10\)、\(q=5/4\)，并定义三个一维正积分
\[
J_0={4\pi\over4-2q}\int_{\mathbb R}(1+u^2)^{-q}du,
\quad
J_1={4\pi\over4-2q}\int_{\mathbb R}(1+u^2)^{-q}
 [\log(2\sqrt{1+u^2})]^qdu,
\tag{16}
\]
\[
M_5=8\pi\int_0^4 t^2[\log(10/t)]^5dt,
\quad
K_0=5^{2/5}15^{4/5}2^{6/5}
\left[aJ_0^{4/5}+bsJ_1^{4/5}\right].\tag{17}
\]
这些积分有限。对 v+ 项，上述投影将 Tδ 包含到十五个过原点直线的半径 2δ 的 Euclidean 圆柱并。沿圆柱轴、三维横截面先积分，分别给 (16) 乘 \((2\delta)^{4-2q}\)。对数轴与圆柱轴不同时，在每个径向球面用层蛋糕：对数的超水平集是同轴的两球帽；任意同面积集合与超水平集的交集不超过两者面积的最小值。因此其积分不超过对齐两球帽的积分。先按十五圆柱求和，再用 Minkowski，得到中央项范数至多 \(K_0\delta^{6/5}\)。这一步支付全部角向对数，不把管内的最小 ρ 错置为正数。

对另外十五个 v，在正面有 \(|x-v|\ge1\)，故 \(r_v^{-2}\le5\)。pair-plane 的面内距离 d 与真距离满足 \(d\le\sqrt5\rho\)，而 \(r_v\le\sqrt5\)，从而 \(\log(2r_v/\rho)\le\log(10/d)\)。在四维 face cube 内，沿该仿射线的投影区间长至多 2，横向距离至多 4，所以完整面上 \(\int\log^5(10/d)\le M_5\)。Hölder 仅用于这个显式对数函数，得到
\[
\boxed{
\int_{T_\delta}G^{5/4}
\le f^{5/4}\left[
K_0\delta^{6/5}+75aV(\delta)^{4/5}
+75bsM_5^{1/5}V(\delta)^{3/5}\right]^{5/4}.}\tag{18}
\]
这里 75=15×5，所有方向都计入。Minkowski 前的每个项可先截断，再用单调收敛，不假设待证的管可积性。由于 V≤108δ³，(18) 给显式正积分常数形式的 \(O(\delta^{3/2})\)，另两项的范数分别为 \(O(\delta^{12/5})\)、\(O(\delta^{9/5})\)。没有进行 J1 的数值积分、半径选择或 1e−6 比较；(18) 本身不表示 1/1024、1/2048 或其它指定半径已可用。下一批可以加强粗常数、认证这三个标量并报告实际通过范围，不能把本批意图写作已运行结果。

## 7. 恢复证据、实际运行与归档边界

原登记在计算前完成。原 `check.py` 的 SHA256 为
`c6f18a8524ed18e0c0986fbd20f030acdd4d663bcb3fc35136067dce8adbdfb6`；登记 SHA256 为
`f21af7692026df9630d72e0a2fdc95aa60e7dd392c77251dc622e7aaf9811d24`；已有结果 SHA256 为
`e751cf4f92c3d0b0695c7ebe3ae4739c27210149bb69df5361d4b1c93fbb9faa`。

保存结果记载原 PID 55532，两个完成阶段分别耗时约 0.00691675 秒、0.000795417 秒；原 checkpoint 与结果字节相同。记录包含完成状态与全部输出，但没有保存原 shell exit code，因此本报告不补写“原进程退出 0”。几何与标量结果的数学含义限于脚本实际核验的有限对象；(1)、(3)、(18) 的证明是本报告的解析推导，不声称脚本进行了形式化证明。

保存的 256-bit Arb 球给
\[
A_0=0.00264698100754130319\ldots,\quad
A_1=0.000356332479709168660\ldots.
\]
完整球与 dyadic 外端仍保存在原结果中；小数只供显示。另登记 `recovery_registration.json`，实际运行 `verify_saved.py`：只用标准库核对本地和七个外部输入 hash、checkpoint 字节、八个常数已存 dyadic 外端的正性顺序、四个预选 rational cap。没有导入生产器，没有新的 Arb、旧几何、Green 或管标量调用。新审计的运行回执见 `recovery_audit.json`。

029/034/017 是本报告相应完整模型、管分类和真顶点表的只读输入；登记中其它历史 α=3/2 报告仅为原路线背景，不把其核、常数或奇异集结论移植到本 Cauchy 证明。来源路径和 hash 全列于原登记及审计。

再次核验已存证据的入口为
```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-5-bounds/upper_singular057/verify_saved.py
```
`check.py` 是原有限批的复现源码，不是恢复阶段应自动重跑的命令。所有新文件仅在本目录；未修改旧锁定输入，未提交或推送，未使用 MATLAB。白名单、各文件字节数和 hash 见 manifest。
