# 029 U10两basis联合修正：全域条件接口独立审计

2026-10-25。结论：锁定文本 `upper/u10_joint_certificate.md` SHA256 `7d3ce12b4cafc3db7d07f3c8c9359f165a777d35f462a48d048f4b46eefae62e` 的三维充分接口通过；未发现需要修补的数学缺口。该结论是条件式：必须另外给定 w0,w1,D，并证明全域非负、全cube支配及剩余有界源域的联合非负。这里没有认证任何可行权重、区间覆盖、成本下降或新上界；有界域归约也不保证某种有限区间算法必会终止。

适用生成元严格限于原四个非单位体对角方向 v=(1,±1,±1)、各权1、L=−Σ|v·D|，以及原b0单位球和六个半径1/2、中心±3e_i/4的小球之和b1。本轮新加入面对角方向的frame不继承这里的具体常数；四维也没有本审计的这些具体常数证书。cb10点态μ区间、U10之外的候选/LP/网格不在本批。

## 来源和实际执行

已先登记 `029_registration.md`，五个来源锁在 `029_source_snapshot.json`，包括主011登记SHA `81656f9ab18ca9c08e1baaa074fa3fe972a317ea9887ba3cca1eb68111acf8c0`、011输入snapshot SHA `6376f1dd7d4de04ffaed7795e7adefce27013d7e2ebd2c8bfcc665ad851c0b35`。继承026已经独立通过的真实frame/正部势/概率源/全L1与共同尺度比较，以及027的真实Jensen下界；这些旧认证未重跑。没有把本轮可变ledger当算术输入。

独立Python第一次运行即通过，内部约0.125秒，无实现修补。代码不导作者程序，实际用SymPy1.14.0重建全部四个basis、精确有理多项式和常数，再以python-flint/Arb192bit检查固定实常数余量。完整准确结果在 `029_results.json`。来源收尾锁定以manifest为准；本审计不修改作者稿。

## 支持、C1,1、质量及源测度

六小球中心的15个距离平方中，12个为9/8，3个为9/4，均严格大于球直径平方1；闭球确实互不相交。b1≤1，b0≤1，因此B=w0b0+w1b1≤W=w0+w1，不可误用6W。b0与小球可以重叠，这不妨碍该上界或积分的线性。

支持在半径a=5/4的闭球内。继承真实G≥5/(16π³|x|²)以及log2<1，半径3/2闭球上G−c>1/(72π³)，故支持至Ω^c距离至少1/4；支持上G−c>3/(40π³)。这些是全域下界，不是取样。若点在frame线上G=∞，同样满足所需不等式。

单位球内Hessian的径向、切向本征值为−4+12r²和−4+4r²，r²∈[0,1]，其绝对值不超过8。球面上值与gradient均连续接零，所以全空间gradient Lipschitz≤8。半径1/2缩放后≤32，六球闭支持分离，沿任何线段逐片积分Hessian界即可给b1的同一全局常数32，而不是6×32。

更关键的是，C_c^{1,1}直接给方向二差分的sup界

\[
 |b(x+sv)+b(x-sv)-2b(x)|\le\min(C|v|^2s^2,4\|b\|_\infty),
\]

并给相应L1界，使用紧支撑弱Hessian的L1范数。除以s²后全s可积。支配收敛因此给Lb为连续有界普通函数，L1的对应支配/Fubini再给Lb∈L1且∫Lb=0。它与分布定义相同，没有遗漏球面atom或球外跳跃尾。这里的“无atom”来自全非局部积分估计，不只是局部Laplacian直觉。

独立径向积分得到M0=32π/105，M1=6(1/2)³M0=8π/35=(3/4)M0。平移和跳跃变量s=Rr给order1缩放L[b((x−z)/R)]=R⁻¹Lb((x−z)/R)。因此F=H−B为even、L1，且

\[
 LF=(\mu-LB\,dx)-\delta_0,\quad
 \nu=\mu-LB\,dx,\quad \nu(\mathbb R^3)=1,\quad\nu(\{0\})=0.
\]

正概率μ和无原点atom来自026，LB为L1密度；所以此处确实是有限signed测度。所有正负贡献须先联合求LB，不可逐基重复付μ。

## Ω内连续密度与Ω外自动结算

026给G在线外有限连续、在线上为∞且lower semicontinuous。因此h=min(G,c)全空间有界连续，在Ω={G>c}的每个点附近恒c（包括奇线与原点）。Q=(c−G)+连续，并在每个Ω内紧集附近为零。由μ=−Lh=LQ可得Ω内的连续普通密度

\[
 m(x)=\pi^{-1}\sum_v\int_{\mathbb R}Q(x+sv)/s^2\,ds.
\]

在距Ω^c至少δ的紧集，小跳|s|<δ/√3统一为零，剩余用c/s²作DCT支配。故q=m−LB连续，点态条件和a.e.条件在Ω等价。四方向完整|s|>S尾位于[0,8c/(πS)]；系数包含两条尾，未漏负s。

Ω外的结算成立，但必须按文本保留全域F≥0这个先决条件。令Z_-={G<c}，它开放且与支持有正buffer，F在该开集恒零。对φ≥0、suppφ紧含Z_-，统一小跳不贡献；剩余大跳的<LF,φ>等于非负F与非负φ的有限率卷积。用F∈L1和Fubini可严格控制全部积分，得到LF≥0于Z_-。这里δ0不在Z_-，故ν在Z_-非负。并未假称Lb在B支持外消失。

其余Z0={G=c}由G的正齐次degree−2处理。除有限条奇线外，每个射线在该level至多有一个半径，极坐标Fubini给Z0的三维Lebesgue测度零；奇线本身也零。于是LB dx不赋质量给Z0，ν|Z0=μ|Z0≥0。即使μ在此level存在奇异部分，也已是非负，毋须额外正则性。结合Ω内q≥0，得到全空间ν≥0。不能把“外部自动”理解成未证(A)时也成立。

## 近奇线log屏障：重建正常数

四选三的所有basis独立准确核det=±4，逆矩阵每行norm²恰1/2≤1；余下第四方向的basis坐标全为±1。写x=t v+z、z⊥v、|t|≤A、ρ=|z|∈(0,1)，以Cauchy时间u∈[ρ,1]条件于余下变量|C_u|≤u，其概率恰1/2。

两条transverse坐标绝对值≤ρ+u≤2u，密度各≥1/(5πu)。纵向坐标绝对值≤A+ρ+u≤A+2，密度≥u/[π(1+(A+2)²)]。乘上概率1/2和det因子1/4，得到

\[
 p_u(x)\ge\frac1{200\pi^3(1+(A+2)^2)u},\quad
 G(x)\ge C_A\log(1/\rho).
\]

所以给定M>0时，ρ≤exp(−M/C_A)确实保证G≥M；ρ=0继承G=∞。取M≥c+W+D/8可以支付近line的值支配。它只提供值和small-jump零区屏障，不能代替q的联合非负证明，也不能据此删除nearline源cells。

## 远尾全域归约：逐个常数复核

设S为四条frame线之并。若dist(y,S)≥d，且四个独立Cauchy参数中至多一个大于δ=d/(8√3)，余下三参数导致距某一frame线至多3√3δ=3d/8<d，矛盾。因此必有两参数大于δ。六pair分别补成det4 basis：两条被限制的Cauchy密度各≤t/(πδ²)，第三条不限制的密度≤1/(πt)，余下一条积分概率至多1，故每pair密度≤t/(4π³δ⁴)。求和并与无条件上界p_t≤(4π³t³)⁻¹取min，精确积分得

\[
 G(y)\le48\sqrt6/(\pi^3\operatorname{dist}(y,S)^2).
\]

pair-law密度分解可先按测度/a.e.表述；真实p_t连续，对距离严格内点成立后以d递增极限得到显示的闭边界版本，无额外常数。这里没有把时间概率误当空间密度。

对x=t v+z属于半径a=5/4的line tube，利用原frame和B的even/cube对称性取t>0、v=(1,1,1)，并取w=(1,1,−1)。四条完整frame线对应的dist²(v+τw,span方向)是

\[
\tfrac83(\tau^2+\tau+1),\quad\tfrac83(\tau^2+\tau+1),\quad\tfrac83,\quad\tfrac83\tau^2.
\]

τ∈[1,2]的最小值确为8/3；不是只检查两条线。t≥5/2时，扰动z/t使distance仍至少√(8/3)−1/2>1。因此s∈[t,2t]时G(x+sw)≤48√6/(π³t²)。

log2>2/3来自1/x在[1,2]的严格Jensen，√6<5/2。所以t≥60时，上界乘π³严格小于1/30，而(c/2)π³严格大于5/144，余量大于1/720。于Ω中可在单一横向direction上积分(c−G)≥c/2，得到m(x)≥c/(4πt)。

远点在B支持外，LB为非负jump-in。每direction与半径a球的参数截段长度≤2a/√3；命中时|s|≥|x|/√3−a/√3≥t−a/√3≥t/2。四方向完整求和给

\[
0\le LB(x)\le\frac{32aW}{\pi\sqrt3t^2}=\frac{40W}{\pi\sqrt3t^2}.
\]

此为包括全部球、两个s方向的外尾界，六小球没有再重复乘6。取T=max(60,160W/(√3c))，则m−LB≥0。对|x|≥R*=√3T+a、|z|≤a，三角不等式|x|≤√3|t|+a给|t|≥T。因此所有Ω内远尾tube均已解析支付；tube外各line不碰支持，LB=0；Ω外另由上一节支付。剩余源条件确实只需在Ω∩T_a∩B(0,R*)中证明。

该剩余域有界，并非自动有正slack、紧离Ω边界或拥有已生成的有限cell cover。未证Ω归属的cells不能套m公式；文本保留此限制是必要的。近线屏障、积分尾界也未替代这项条件。

## 全L1/尺度结论及终止范围

若(A)(B)(C)确已认证，K=F/D非负、even、L1，支配1_cube/8，原点外源为ν/D≥0。026已审finite-frame比较的假设全部满足，得到同一坏集/同一u、σ下所有共同r及一般f∈L1的弱型常数

\[
 C(K)=(\int H-w0M0-w1M1)/D.
\]

该逻辑不把D当成已经认证的8(Gedge−c)，不假设g=edge，也不预设成本改善。恒定Lebesgue、原方向权重与普通密度/Fubini/球质量/level零体积均保留；不能只换背景或方向而沿用数值常数。

最小复核命令（主导入本审计后，stdout-only；可选--save）:

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-finite-dim-bounds/audit/029_check.py
```

审计clone运行加 `--project /Users/zhengzhihao/Desktop/T/Math-Research-finite-dim/cube-finite-dim-bounds`。MATLAB脚本只交主串行作有限算术交叉，实际状态见 `029_matlab_receipt.json`；未用于替代这些分布/连续性证明。共享SymPy/Arb数学内核及026/027既有证明明确披露，无作者算术程序导入、无LP/网格/候选重新求解。最终manifest完成后本有限批停止。


主串行MATLAB已实际exit0（session3132），以eval(fileread)执行同SHA a79d155171219b25a2a3deac94838ece9c7e2683167fc95461b2bddfa0ebb7d4脚本，未clear/path。审计已读log，SHA13efd998d63475112708d4fc0620c064c3ba5301b46788c0b1382bc7c65b7795。实际四basis/质量/距离/粗余量assert均通过；展示supportgap=.0040682774611725068、buffergap=.0020973503569214264、tailgap_t60=.00013768400895848749。只作算术cross，不将float冒称落在192bit窄包围内。收尾再核五来源未变。
