# 033 全域minimum：136闭patch/59976节点独立fresh认证

2026-10-25。结论通过：固定三维body权1、face权τ=1/3的真实Green满足整个闭cube上的逐点下界

\[
G(x)\ge131/25000=0.00524\qquad(x\in[-1,1]^3).
\]

这比原目标13/2500强。独立fresh384bit计算的最小**认证下界**约为0.00524264797871796797512876059812575296875438805，来自patch[15,7]；并非证明真实minimum恰取该数或唯一取在该patch。本批只认证g接口，I及成本未在此闭合，不宣称U3<8。

## 锁源、实现和真实运行

来源先锁 `033_source_snapshot.json` 的14项：首manifest `320bdba75d2273568d55e06c2cc888bd310daa38aea60cff79486e0c8cd0c7cf`、Cheb manifest `3e36efe767eb8f87ca78bbee9ab5c640113eb5b76a2da28f8202b528155ecf20`、作者checker `f3f57f72959b6350337ed95a702d137b8e485d9f640fb8b2dff67168bc6a5db8`、report `f6be7132f37aebecb2c3c872cf7bd0af7a56a9c85de7a638d04bd5eee1b52bdf`，以及真实cb11模块 `ffbacd4f1a466f8f20d11dc59c458be2980b45be964751b8ea1a980cbbf1d849` 等。收尾再核全部未变。

未运行会写回作者目录的main，也未读取或写入其节点缓存。本审计新实现完整59976个rational-point real Hermite DD、DCT和Bernstein转换。只读取032已审锁模块的48对ray/rate/weight准备数据；没有调用作者point evaluator、Cheb checker或selfcheck。新点值先将ε²+b²等内层式保持准确Fraction，再送Arb；Bernstein幂系数由SymPy的精确Chebyshev多项式展开构造，独立验证21个完整多项式恒等式，而非只比较其21个点。共享的公式、准备数据、SymPy1.14.0、Fraction和python-flint/Arb内核明确披露，不称无共享软件内核。

首轮在441个新节点及所有代数预检后，因Python3.15的Fraction不直接接受flint.fmpq，在下端类型转换处exit1。该实现错误完整记录于 `033_execution_history.json`，未当作数学失败或成功证据。仅修为准确有理字符串转换；未改任何数学参数。第二轮从零完整重算59976节点，实际exit0，耗时62.0211秒；包括失败首轮总共60417次新点值计算。完整成功轮没有缓存读写、没有分支、没有参数搜索。

固定384bit、ε=1/50、degree20、rho3、128bit dyadic节点替代、16网格136闭patch均未改变。成功轮节点ball流SHA为 `4f45e2d83f49096acd58b66ea763498f9a0efdcb3c76580da71d8b8d4cfd71d7`；各patch准确下端及强阈值余量见 `033_results.json`。hash是运行记录，不替代下面的误差证明。

## 截断核下界与real DD evaluator

记本批Hε=∫ε∞ptdt，区别于031的正部势(G−c)+。真实pt≥0直接给Hε≤G，故minimum下界不需要支付G−Hε的短time误差。不能把这一步误当成I上界也无需短time余量。

032已完整审定96锥、48反锥pair及截断time实DD恒等式。每pair三节点b_j=r_j·(1,y,z)/k_j，W=2/Πk_j，且

\[
H_\epsilon(1,y,z)=\frac1{8\pi^3}\sum W\,f_\epsilon[b_0,b_1,b_2],
\quad f_\epsilon(b)=\tfrac12\log(\epsilon^2+b^2).
\]

新实现按准确Fraction判断等节点；重复两节点用[f(b)−f(a)−f'(a)(b−a)]/(b−a)²，三重用f''(a)/2。f'=a/(ε²+a²)，f''/2=(ε²−a²)/[2(ε²+a²)²]；零节点是正常情况。所有signed项保留，不从ball重叠猜“准确重复”，不裁剪负项。ε>0时全部式有限，原G在body/face线上的∞不产生覆盖缺口。

## 全纯tube与节点位移：正常数及max而非sum

逐ray准确核k_r≥3(|r_y|+|r_z|)，最小比值确为3。借032锥上线性及三角不等式，推出全频率

\[
\psi(\xi)\ge3(|\xi_y|+|\xi_z|).
\]

pt(0)=15/(544π³t³)。令F(y,z)=Hε(1,y,z)。若max(|Im y|,|Im z|)≤s<3ε，Fourier绝对值中额外指数至多s(|ξ_y|+|ξ_z|)≤sψ/3，故

\[
|F|\le\int_\epsilon^\infty p_{t-s/3}(0)dt
=\frac{15}{1088\pi^3(\epsilon-s/3)^2}.
\]

该统一绝对支配在严格tube的紧子集成立，给联合全纯性。real parts任意，不要求ellipse或Cauchy圆盘的实部留在patch内。两个虚部由同一个max s控制，不能额外错误付2s，也没有漏付一坐标。

在任意实点取单坐标Cauchy半径R=1/50，有ε−R/3=1/75>0，故每个实偏导由

\[
D=\frac{15}{1088\pi^3(1/75)^2R}<125.056215
\]

控制。每patch半宽h=1/32。21个true roots cos((2j+1)π/42)中，正10个以384bit Arb确认2^128倍的唯一floor；q_j为该准确dyadic。负根用准确反对称，中央根准确0。因此每根误差<2⁻¹²⁸、每个物理坐标位移<h2⁻¹²⁸。连接两节点的实坐标折线路径给函数值误差≤2Dh2⁻¹²⁸，无须假设无理点间独立。

这项小于2.296921e−38的误差以向上半径加在每个rational点的Arb值上；其后经整个DCT/Bernstein矩阵传播。它不是在变换后仅扣一次未经放大的节点误差。精确cos的DCT矩阵和所需node-value inclusion都同时成立，区间依赖遗失只会扩大包围。

## Chebyshev系数、aliasing及tensor误差

rho=3的每个物理Bernstein ellipse虚部高度s=h(rho−rho⁻¹)/2=1/24，ε−s/3=11/1800>0，得到统一M<11.906178395825。全纯性在其闭包邻域成立，故一维Cheb级数系数|a_k|≤2M rho⁻k，k≥1。

在N=21个roots上，每个T_k的值alias为某个±T_l、l≤20，或零（由cos(kθ_j)对2N的反射/符号周期）。因此插值后的sup≤1，高频项误差≤2。绝对收敛级数逐项相加给

\[
E_1\le4M\sum_{k=21}^\infty\rho^{-k}
=4M\rho^{-20}/(\rho-1).
\]

DCT cardinal式为ell_j(t)=N⁻¹[1+2Σ_{k=1}^{20}cos(kθ_j)T_k(t)]，所以Λ≤41。对tensor插值I_yI_z，恒等式(I−I_yI_z)=(I−I_y)+I_y(I−I_z)给全闭patch误差

\[
E\le(1+41)E_1=84M/3^{20}
<2.868313237155\times10^{-7}.
\]

这里已经包含两个坐标的插值误差；不再漏付第二维，也不重复扣两次E。该推导只需各方向ellipse上的统一M，实际tube对其乘积域也成立。

## DCT及准确Bernstein转换

采用a0=N⁻¹ΣV_j、ak=2N⁻¹Σcos(kθ_j)V_j的实际系数约定。有限roots的cos正交直接给逆关系；程序另外实际核了441个Arb矩阵entries包含对应单位阵值且半径<2⁻³⁰⁰。此有限数值检验辅助标准正交身份，不代替其证明。

设T_k(2u−1)=Σ_l p_kl u^l。独立SymPy提取准确有理p后，使用

\[
B_{jk}=\sum_{l\le\min(j,k)}p_{kl}\binom jl/\binom{20}l.
\]

其依据是u^l=Σ_{j≥l}[binom(j,l)/binom(20,l)] b_j^{20}(u)，由二项式恒等式直接成立。审计逐个k对完整多项式反向展开，21个恒等式全部准确为零。A是DCT，B是此转换，故二元系数矩阵为(B A)V(B A)^T，行列方向与真实node网格一致。

闭单位方形上Bernstein基非负、总和1，因此多项式下界为441个系数的最小下端。本审计先把每系数下端转为准确Fraction取min，减**完整Arb E**后再向下取端点；这个下端不大于“最小系数下端−E上端”。节点及矩阵四则误差都已在球内，另扣E恰一次。全部136个root直接严格超过131/25000；未执行作者的depth8 de Casteljau分支，本结论不依赖其未使用分支实现。

成功轮最大node半径约2.296919896e−38，最大Bernstein半径约2.042668306e−26。它们远小于最后强阈值余量，但接受依据是完整准确下端比较而非这些展示量。

## 强下界、闭覆盖与全cube

每patch为准确闭矩形[i/16,(i+1)/16]×[j/16,(j+1)/16]，0≤j≤i≤15。给任意0≤z≤y≤1，取i=min(15,floor16y)、j=min(15,floor16z)即有j≤i，故落入对应闭矩形；所有边界、diagonal与顶点均包含。diagonal patch核整个方形，额外区域不会造成漏洞。独立index集合完整且无重复，共136。

对非零cube点，取m=max|x_i|≤1，再用031已审cube群对称将x/m送到(1,y,z)、0≤z≤y≤1。**使用的是G的**正齐次G(x)=m⁻²G(1,y,z)，不是给截断Hε误用齐次。由Gface≥Hεface≥strong，得整个cube同一下界；原点G∞也满足。

最小认证下端的完整准确有理数在 `033_results.json/minimum_H_lower_exact`，并且

\[
\min_{136}\underline H-131/25000>2.6479787179\times10^{-6}>0.
\]

其patch索引[15,7]对应y∈[15/16,1]、z∈[7/16,1/2]。这里比较的是可复核dyadic下端与有理数，不是sampled minimum，不将node值最小冒充连续下界。

## 复现、MATLAB、共享与停止

主导入后，完整fresh复核命令（默认不写结果，--save才写本审计结果）：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-finite-dim-bounds/audit/033_check.py --save
```

审计clone运行加 `--project /Users/zhengzhihao/Desktop/T/Math-Research-finite-dim/cube-finite-dim-bounds`。每次必重新计算全部59976节点，不开放读取旧cache选项，不写作者文件；开始/结束都核来源hash。完整point/rate理论依赖032，真实过程与全L1/all-r框架依赖031；这些已审证明不在033重跑，也不因本g接口单独更新U。

MATLAB脚本仅由主串行核degree0/1/7/20完整Bernstein恒等式及M/D/E/位移常数，实际状态见 `033_matlab_receipt.json`。它不宣称第二次完整59976节点重算。共享prepared rays/rates及Arb内核如实披露，新point evaluator和多项式构造不是作者checker的重复调用。

本有限批只有一次有记录的实现类型修复，随后完整fresh成功；没有精度/degree/网格/ε/ρ/阈值放宽、分支或搜索。严格结论只为本τ1/3的全cube g≥131/25000。I/J/短time上界及新U3仍须独立完成，后续034另登记。本批最终manifest完成即停。


主串行MATLAB已实际exit0（session11147），同SHA20b23e380911a413a9ba9eb17b544de710bd5fd5b9dae261146237c9ad593e44。审计已读stdout，degree0/1/7/20准确转换assert通过，M=11.906178395824435、D=125.05621405336429、E=2.8683132371546148e−7、nodeerror=2.2969198930461325e−38；log SHA61b34a4d93cac482dafef44170560218dc2292b92d5c3317d3e4547cc350d3ab。float展示不冒称处于384bit窄包围内。收尾14来源及成功checker/result关联hash均未变，033有限批已停止。
