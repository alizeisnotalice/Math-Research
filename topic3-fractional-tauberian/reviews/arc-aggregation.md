# 第 22 轮：一般四角的独立尺度共同不可行证书

## 结论

设 \(\alpha=1-t\)、\(c=2\alpha^2-1\)。四角
\((\pm r,\pm s)\) 的对角 pair 容量多项式为

\[
 G_{r,s}(z)=(z-2r)(z-2s)-2\alpha z+1.
\]

在边界

\[
 (\alpha-2r)(\alpha-2s)=c,
\tag{1}
\]

上有 \(G_{r,s}(\alpha)=0\)。对角线附近的一个真正二维参数窗可以延拓第 21
轮的四 patch 证书：取

\[
\theta=\frac\pi4+\zeta t,\quad |\zeta|\le\frac18,
 \qquad
 r_b=\frac{\alpha-\sqrt{c\tan\theta}}2,
 \quad s_b=\frac{\alpha-\sqrt{c\cot\theta}}2,
\]

并沿真实法向内移

\[
r=r_b-d_0t^3\cos\theta,\qquad
s=s_b-d_0t^3\sin\theta.
\tag{2}
\]

正式引理取中心内移 \(d_0=1/64\)，法向半厚与切向半宽常数
\(h_0=k_0=1/4096\)，并取 \(q=\alpha+et^2\)、\(e=1/8\)。存在绝对常数
\(t_0>0\)，使得对 \(0<t<t_0\)，四个以真实法向厚 \(2h_0t^3\)、真实
切向宽 \(2k_0t^{3/2}\) 的对称 patch 同时满足：

1. 每块完全位于相应的 pair body \(T\) 内；
2. 对任意从四块中各取一点的四元组，任何共同输入的四个独立见证尺度均小于
   \(M=2q-\alpha\)；
3. 四个边长 \(\alpha\) 的核心方体及尺度窗的质量下界严格大于一。

因而没有质量一输入能同时严格见证四块 patch。每块面积为正的
\(\asymp t^{9/2}\)。这里的结论是一个可用于沿对角边界连续布置的局部证书；它并不
自动覆盖所有仅满足 \(r/t,s/t\in[1/3,2/3]\) 的固定偏斜四角。

## 对角 pair 的精确机制

若两个相反中心的方体边长为 \(\ell,m\)，令 \(z=(\ell+m)/2\)。每个严格见证
捕获量大于 \(\alpha\ell\)、\(\alpha m\)，而交叠面积不超过
\((z-2r)(z-2s)\)。因此共同质量一输入必满足

\[
 G_{r,s}(z)>0.
\tag{3}
\]

在 (1) 的边界点 \((r_b,s_b)\) 上，\(G_{r_b,s_b}\) 的另一个根为
\(\alpha+2(r_b+s_b)\)，故

\[
 G_{r_b,s_b}(z)=(z-\alpha)\bigl(z-\alpha-2r_b-2s_b\bigr).
\tag{4}
\]

对 (2) 的 WolframScript 1.14.0 符号极限为

\[
 \begin{aligned}
 \frac{G(\alpha)}{t^3}&\longrightarrow2\sqrt2d_0,\\
 \frac{-G(q)}{t^3}&\longrightarrow2e-2\sqrt2d_0,\\
 \frac{-G(\alpha^{-1})}{t^3}&\longrightarrow
 2(1-\zeta^2-\sqrt2d_0),\\
 \frac{(\alpha+2r)(\alpha+2s)-4M(M-\alpha)-1}{t^2}
 &\longrightarrow2(1-4e-\zeta^2).
 \end{aligned}
\tag{5}
\]

在 \(|\zeta|\le1/8\) 时，所选 \(d_0,e\) 使四个极限都一致正。特别地，
\(r/t\to(1-\zeta)/2\)、\(s/t\to(1+\zeta)/2\)，它们落在
\([7/16,9/16]\subset[1/3,2/3]\)。

由 \(G(q)<0\) 与 \(G(\alpha^{-1})<0\)，并使用任一见证边长均位于
\((\alpha,\alpha^{-1})\)，(3) 强迫其平均尺度小于 \(q\)。故每一尺度均小于
\(M=2q-\alpha\)。最后一项正性正是四个核心方体的并面积减去四个尺度窗预算后
仍大于一的量。

## 真正法向和真正切向

记

\[
 u_\theta=(\cos\theta,\sin\theta),
 \qquad \tau_\theta=(-\sin\theta,\cos\theta).
\]

东北角 patch 可取为

\[
 P_{++}=\{(r_b,s_b)-d_0t^3 u_\theta+n u_\theta+v\tau_\theta:
 |n|\le h_0t^3,\ |v|\le k_0t^{3/2}\},
\tag{6}
\]

其余三块由坐标反射得到。对相反 patch 的扰动写为
\(h_0Nt^3u_\theta+k_0Vt^{3/2}\tau_\theta\)，其中 \(N,V\in[-2,2]\)。在 \(z=\alpha\) 处，\(G\) 的线性项
沿 \(\tau_\theta\) 完全消失，而法向项为

\[
 2\sqrt2d_0t^3-\sqrt2h_0Nt^3+O(t^4).
\]

切向二次项为负，大小为
\(-k_0^2V^2t^3/2+o(t^3)\)。代入所列显式常数后，
\(G(\alpha)>0\)，同时 patch 位于 \(T\) 的内侧。这个配比正是
法向内移与切向平方同为 \(t^3\) 量级，而不是把法向固定成 \((1,1)\) 后再作
\(O(t\eta)\) 的粗估。

在 \(z=q\) 时，切向一阶项还带因子 \(q-\alpha=O(t^2)\)；在
\(z=\alpha^{-1}\) 时，它是

\[
 -2(\alpha^{-1}-\alpha)(\cos\theta-\sin\theta)\sigma
 =O(t^2\sigma)=O(t^{7/2}),
\]

因为 \(\theta-\pi/4=O(t)\)。这严格小于 (5) 中 \(t^3\) 的负余量。法向扰动
由中心内移 \(d_0\) 和邻域常数 \(h_0,k_0\) 的显式余量吸收。因此 (5) 的
严格不等式稳定地延拓到整个四 patch。

## 并面积与为何不需要伪 \((1,-1)\) 切向

对先把真实切向替换为四十五度切向的模型，记横纵跨度为 \(X,Y\)，则

\[
 w=X+Y-2(r+s)\ge0.
\]

四角 deficit 的总和是 \(2w\)。按外包矩形的四角归属并估计，可得

\[
 XY\ge4rs+(r+s)w-4|r-s|\eta-4\eta^2.
\]

在 \(\alpha\ge3\max(r,s)-\min(r,s)+4\eta\) 时，和第 21 轮相同的角区覆盖计算给出

\[
 \left|\bigcup_{j=1}^4Q_\alpha(x_j)\right|
 \ge(\alpha+2r)(\alpha+2s)-4|r-s|\eta-4\eta^2.
\tag{7}
\]

真实切向和四十五度切向相差 \(O(t)\)，所以 (7) 还只损失
\(O(t\eta)\)；法向厚度另损失 \(O(t^3)\)。对
\(\eta=kt^{3/2}\)，这些损失是 \(O(t^{5/2})+O(t^3)=o(t^2)\)，被 (5) 的一致
\(t^2\) 质量余量吸收。故真正法向、真正切向既避免了 pair 端点的伪线性误差，也
不破坏面积矛盾。

## 适用范围与尚未展开的统一化

本报告严格记录的常数保证 \(|\zeta|\le1/8\) 的窄带。此前把较大比例偏斜说成会使
法向和 \(\pi/4\) 相距常数量，是不正确的：只要 \(r,s=O(t)\) 且 (1) 成立，
\(\alpha-2r=1+O(t)\)、\(\alpha-2s=1+O(t)\)，从而
\(\theta-\pi/4=O(t)\)。更精确地，若 \(r/t,s/t\) 有有限极限，则
\[
 \frac rt=\frac12-\frac\zeta2+O(t),\qquad
 \frac st=\frac12+\frac\zeta2+O(t),
\]
其中 \(\zeta\) 可以是固定常数；这正对应固定的偏斜比例。

因此在 \(z=\alpha^{-1}\)，真切向导数为
\[
 -2(\alpha^{-1}-\alpha)(\cos\theta-\sin\theta)=O(t^2),
\]
而不是 \(\Theta(t)\)。对 \(t^{3/2}\) 宽的 patch，其一阶扰动为
\(O(t^{7/2})=o(t^3)\)。本文只使用上述窄带证书；较宽参数区间不在本轮结论之内。

## 最终复核：扰动端点与量词

正式邻域以 \(p_\theta-d_0t^3u_\theta\) 为中心，法向、切向半宽分别为
\(h_0t^3,k_0t^{3/2}\)。从一对相反邻域各选一个点后，令
\(N,V\in[-2,2]\) 为两个扰动的归一化和。实际以 WolframScript 1.14.0 计算
\(G(z)=(z-2r)(z-2s)-2\alpha z+1\) 的 \(t^{-3}\) 极限，得到
\[
\begin{aligned}
 G(\alpha)/t^3&\to2\sqrt2d_0-\sqrt2h_0N-k_0^2V^2/2,\\
 -G(q)/t^3&\to2e-2\sqrt2d_0+\sqrt2h_0N+k_0^2V^2/2,\\
 -G(\alpha^{-1})/t^3&\to
 2-2\sqrt2d_0+\sqrt2h_0N+k_0^2V^2/2-2\zeta^2.
\end{aligned}
\]
与正式引理的三个式子一致。代入 \(d_0=1/64\)、\(h_0=k_0=1/4096\)、
\(e=1/8\)、\(|\zeta|\le1/8\)、\(|N|,|V|\le2\)，其最坏情形仍严格正；
故紧盒一致性给出某个 \(t_0>0\)，并不预设数值如 \(.01\)。
四个可直接取到的下界（前三个端点式及质量余量）数值分别为
\(0.0435035206\)、\(0.2051152922\)、\(1.9238652922\)、\(0.96875\)；
它们只是对上述精确符号不等式的交叉核验，并非证明的替代。

量词须按如下顺序使用：先任取四个邻域中各一点，并反设同一质量一输入对这四点都有严格见证；
再把相反点配对。每一对的端点符号强迫该对两尺度的平均小于 \(q\)，而每个尺度本来严格大于
\(\alpha\)，故**各个**尺度均小于 \(M=2q-\alpha\)。两对合并后才应用四方体并面积和尺度窗
预算，得到质量严格大于一的矛盾。另一方面 \(G(\alpha)>0\) 只负责所选的**每个点**落在
\(T\) 内；它不应被误写成对见证尺度的结论。

若在正文表述“任意共同平移”，该不变性只适用于“同一输入不能同时严格命中四邻域”的结论；
平移后的邻域当然不必仍包含于固定的 \(T\)。


# 第 22 轮：连续四角 circuit 的聚合几何复审

本记录只处理聚合、覆盖和随参数平移的几何逻辑。它不把尚未证明的局部四角证书当作定理。

## 结论

在一段长度 \(\asymp t\) 的对角侧弧上，若每个参数都有统一的四 patch 共同不可行证书（法深 \(\asymp t^3\)、切宽 \(\asymp t^{3/2}\)、单 patch 面积 \(\asymp t^{9/2}\)），则已有凸几何可严格给出聚合传递：\(\delta=\mathcal U-|H|<\kappa t^4\) 时，正比例的 circuits 的四 patch 都落在适当平移的凸包内，且这些平移后的 patch families 可取互不相交。孔洞计数遂迫使至少一个 circuit 被同一 \(H\) 的四点同时命中，矛盾。

因此，四阶路线的剩余实质工作是一般侧弧四角 certificate 的统一严格证明；并非聚合或平移依赖本身。下述“聚合引理”应在该证书完成后才能进入正文。

## 1. 参数、网格和孔洞计数

令 \(I\) 是避开侧弧端点的一段法向区间，\(|I|\asymp t\)，取网格尺度
\(\ell=t^{3/2}\)，相隔 \(L\ell\) 取点，数目
\[
 M\asymp t/\ell\asymp t^{-1/2}.
\]
每个 circuit 的四个证书 patch 面积至少 \(a_0t^{9/2}\)。若有
\(m\ge c_0t^{-1/2}\) 个好参数，且它们移回原坐标后的 patch 两两不交，则若每个
circuit 都没有被 \(H\) 四点同时命中，至少每个 circuit 有一个整 patch 落在
\(A\setminus H\)。从而
\[
 |A\setminus H|\ge ma_0t^{9/2}\ge c_0a_0t^4.
\tag{A1}
\]
与 \(|A\setminus H|\le\delta<\kappa t^4\)（\(\kappa<c_0a_0\)）矛盾。

## 2. 对称化 deficit 的局部 envelope

令 \(S=(A-A)/2\subset T\)、\(d=h_T-h_S\ge0\)、\(g=|T|-|S|\)。混合面积恒等式是
\[
 g=\frac12\int d\,(d\mu_T+d\mu_S)
 \ge c\int_I d(\theta)\,d\theta,
\tag{A2}
\]
因为对角弧上 \(d\mu_T=\rho_Td\theta\)、\(\rho_T\asymp1\)。

比单纯 \(L^1\) 更有用的是分布意义的半凹性：
\[
 d''+d=\mu_T-\mu_S\le\mu_T,
 \qquad d''\le\rho_T-d\le1
\tag{A3}
\]
于弧内。若 \(d\ll\operatorname{dist}(\theta,\partial I)^2\)，非负性和
\(d''\le1\) 蕴含每个 supergradient \(p\) 满足 \(|p|\le\sqrt{2d(\theta)}\)：
将半凹上二次支撑函数沿 \(-p\) 方向代入，若右端变负便与 \(d\ge0\) 矛盾。
因此 \(\sqrt d\) 是常数 \(1/\sqrt2\) Lipschitz 的（在有限值处，分布逼近即可）。

设 \(M_I d\) 是区间 Hardy--Littlewood 最大函数。于是若
\(M_Id(\theta)\le\lambda t^3\)，则对 \(|s|\le ct\)
\[
 d(\theta+s)\le a s^2+C_a\lambda t^3.
\tag{A4}
\]
证明如下：若左端大于 \(as^2\)，记其为 \(h\)，则上面的 Lipschitz 性给出长度
\(\asymp\sqrt h\) 的区间上 \(d\ge h/2\)。将该区间同以 \(\theta\) 为中心、半径
\(O_a\sqrt h\) 的区间比较，得到 \(h\le C_aM_Id(\theta)\)；否则直接由 \(as^2\)
控制。全局帽界 \(d\lesssim g^{2/3}=o(t^2)\) 保证上述二次支撑点不碰到 \(I\) 的端点。

弱型最大不等式和 (A2) 给
\[
 |\{M_Id>\lambda t^3\}|\le C\frac g{\lambda t^3}
 \le C\lambda^{-1}\kappa t.
\tag{A5}
\]
故 \(\kappa\) 足够小后，正比例的网格点满足 (A4)。这补足了先前“deficit-shadow”问题：
不能只写平均缺深，必须使用 (A3) 的半凹性。

## 3. Odd 部分的局部控制与参数平移

取 Steiner 平移后的 odd 支撑函数 \(f\)。其 \(H^1\) 能量满足
\(\int_I(|f'|^2+|f|^2)\le C\Delta\le C\delta\)。同时要求
\[
 M_I(|f'|^2+|f|^2)(\theta)\le\lambda t^3.
\tag{A6}
\]
其坏集长度同样至多 \(C\lambda^{-1}\kappa t\)。好点上
\[
 |f(\theta+s)-f(\theta)|\le C\sqrt\lambda,t^{3/2}|s|.
\tag{A7}
\]
在该点及其互补法向作线性平移，使 odd 部分归零。由法向横截的一致性和 (A7)，
支撑函数奇部在 \(|s|\le ct\) 上被
\(C\sqrt\lambda(t^{3/2}|s|)\) 控制；结合 (A4)，选定小 \(a,\lambda\) 即被
局部 \(s^2+t^3\) 支撑余量吸收。远法向仍使用已证的全局 cap-distance 与
\(\|f\|_\infty\lesssim\sqrt\Delta\) 控制。

令这种局部归零平移为 \(b_\theta\)。两参数的差满足
\[
 |b_\theta-b_\phi|
 \le C\bigl(\sqrt\Delta\,|\theta-\phi|^{1/2}
              +\sqrt\Delta\,|\theta-\phi|\bigr).
\tag{A8}
\]
首项来自 \(H^1\) 增量，第二项来自法向基底旋转。故
\[
 \operatorname{osc}_I b_\theta
 \le C(\sqrt{t\Delta}+t\sqrt\Delta)
 =O(\sqrt\kappa,t^{5/2})=o(t^{3/2}).
\tag{A9}
\]
而侧弧物理参数化的曲率半径 \(\asymp1\)，所以相隔 \(L\ell\) 的 circuits 在各自
平移后仍可用固定大 \(L\) 保持分离。不同 \(b_\theta\) 不造成逻辑问题：把每个
包含关系移回原坐标后，仍由同一个 \(H\) 命中这些彼此不交的证书 patch，而共同不可行性
对整体平移不变。

## 4. 条件聚合引理

假定一般侧弧版本的四角 certificate 对所有 \(\theta\in I\) 成立，并且其常数与
上述 patch 尺度一致。取同时满足 (A4)、(A6) 的好参数；(A5) 和其 odd 对应式表明
其个数至少 \(c_0t^{-1/2}\)。第 3 节把每个好参数的四 patch 放进相应平移的 \(A\)，
并保持跨参数分离。第 1 节的孔洞计数因此给出 (A1)。

这就是严格的聚合框架。其需要进一步人工核验的前提只有：一般 \((r,s)\) 侧弧上的四角
共同不可行证书、其法深 \(t^3\)/切宽 \(t^{3/2}\) 的统一常数，以及与邻近 circuits 的
物理分离常数。数值负优化不能替代这些前提的证明。

## 5. 特定连续四 patch 松弛的四阶方法边界

定义 \(\mathcal G_t\) 为所有有界开集 \(S\) 的面积上确界，条件为
\(S-S\subset K_{1-t}\)，并且对每个侧弧参数 \(\theta\) 及任意共同平移 \(b\)，
\(S\) 不会同时命中该参数的四个 certificate patch \(P_{\theta,j}+b\)。这不是 halo
常数：它刻意忘掉“同一输入实现所有点”的其余结构。

在一般侧弧 certificate 成立后，上述聚合证明只使用这些松弛条件。确实，令
\(A=\overline{\operatorname{conv}}S\)，则 \(A-A\subset\overline K\)、\(|A|\le|T|=\mathcal U\)，
并且 \(|A\setminus S|\le\mathcal U-|S|\)。故 (A1) 直接给
\[
 \mathcal U(1-t)-\mathcal G_t\ge\kappa t^4.
\tag{A10}
\]

反向候选为
\[
 S_0=(1-t^2)\operatorname{int}T.
\]
它是有界开集，且
\(S_0-S_0=(1-t^2)\operatorname{int}(2T)\subset K_{1-t}\)。对于任一 \(\theta\)，
相反两个 patch 在 \(u_\theta\) 方向所要求的最小投影差为
\[
 2h_T(\theta)-2(d_0+h)t^3,
\]
而 \(S_0\) 在该方向的宽度为
\[
 2(1-t^2)h_T(\theta)=2h_T(\theta)-2t^2h_T(\theta).
\]
在所取侧弧上 \(h_T(\theta)/t\) 有正的统一下界（中心值趋于 \(1/\sqrt2\)），
且 \(d_0=1/64,h=1/4096\) 时
\(2t^2h_T(\theta)>2(d_0+h)t^3\) 对充分小 \(t\) 成立。因此即使只看一对相反
patch，\(S_0\) 也不能同时命中；此结论对共同平移不变。故 \(S_0\) 是该松弛问题的可行集。

其面积亏损为
\[
 \mathcal U-|S_0|=(2t^2-t^4)\mathcal U(1-t)=4t^4+O(t^5).
\tag{A11}
\]
于是该**特定连续四 patch 松弛**的亏损阶恰为 \(t^4\)。这仅说明用这一类差体和
certificate 约束的聚合法至多达到四阶；\(S_0\) 不被声称为任何最大函数的 halo，
故它既不构造真实反例，也不收窄真实常数的三阶系数 \(10/3\) 上包络。

## 正式稿终审（`pb:arc-maximal`、`pb:quartic-gap`、`pb:quartic-relaxation`）

通过，未发现聚合或松弛范围的阻塞性漏洞。

1. `pb:arc-maximal` 的关键分布不等式
\(D''=\rho_T-\mu_S-D\le1\) 正确。结合 \(D\ge0\) 与内部 \(D=o(t^2)\)，
将半凹上二次支撑式向 \(-p\) 方向代入可得 \(|p|\le\sqrt{2D}\)；从而
\(\sqrt D\) 的局部 Lipschitz 控制给出
\(D(\phi+s)\le as^2+C_a\mathcal MF_D(\phi)\)。端点距离与
\(\|D\|_\infty=o(t^2)\) 的处理足以避免该二次比较越出 \(J_j\)。

2. 混合面积恒等式给 \(\int F_D\le Cg\)，而 Fourier 能量给
\(\int F_E\le C\Delta\)。两次一维 Hardy--Littlewood 弱型估计的坏参数总长皆为
\(O(\delta/(\lambda t^3))\)，所以先定小 \(\lambda\)、再定小 \(\kappa\) 后，
好集确有长度 \(\asymp t\)。正文正确使用从正测度好集的贪心分离选择；没有把测度信息
误换成预定网格上的好点数。

3. 两法向归零平移的局部增量、\(D\) envelope 与近法向二次 gap 的尺度相配。远法向的
\(ct^2\) gap 吸收 \(O(\sqrt\delta)=O(\sqrt\kappa t^2)\) 损失。平移振荡
\(O(\sqrt\kappa t^{5/2})\) 小于 patch 切宽 \(t^{3/2}\)，故固定大间隔贪心选择后，
同象限 patch 不交；不同象限距离 \(\asymp t\)。于是孔洞预算确实给出 \(t^4\) 矛盾。

4. 松弛命题的上界确实只复用差体包含、共同平移禁配和凸包孔洞预算。候选
\(B=(1-t^2)\operatorname{int}T\) 满足 \(B-B\subset K\)；其宽度严格小于相反 patch
所需投影跨度，因而连该对 patch 都不能同时命中。面积差
\((2t^2-t^4)\mathcal U=4t^4+O(t^5)\) 正确。正文明确把它限制为几何松弛，未误作 halo
或真实三阶系数的结论。
