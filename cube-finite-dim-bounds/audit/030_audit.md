# 030 有限方向Cauchy Green的face奇性帽：独立审计

2026-10-25。结论：锁定 `proofs/round011_face_caps.md` SHA256 `b1564bf4c0cb4a6e2f72786ea5ffae30c922ad45c0cb9874aae6bfa55d4fa8de` 的分析引理、五组补基常数和十个帽预算准确区间通过。没有发现数学缺口。它给出face上G幂积分的遗漏预算，不是g下界、完整J值、正部源迁移证明或新上界；本轮新weighted迁移稿和cb11谱锥/Green表示均未在030中独审。

本批固定 n3 body+face、τ=1/10,1/3,1,3，以及n4纯body。所有方向不单位化，body权1、face权τ。零权方向从过程及cmin定义中删除，n4的τ=0是删去face的纯body情形。

## 来源、独立实现与实际执行

先读主011c登记，来源锁在 `030_source_snapshot.json`；主输入manifest SHA `d34ea143f8f399f7becfa70788e3761d57f6584284f842b0f62f24e6d4b63b74`。主Python SHA `6880de4799407f862b8d662ea816cc4ef9bdb126d640140888f7122d5ee6a1aa`、主结果 SHA `7aa6b0cfe3bd91f967cef08380ad90182972af24465705b00c579eb4eefbadd1`。完整七项来源在快照中；审计未导主程序。

独立Python首次运行通过，内部约0.329秒，未修补/重复。以SymPy精确矩阵det重建各unordered pair的所有body补基，再将pair和乘2；与主ordered循环不同。另逐一验证主保存的每个补基确实达到对应最大det，核对方向集、权、重数和全部face投影PSD矩阵。十个确定解析表达式用256bit Arb的正sqrt/cuberoot路线求包围，与主exp(log×power)实现分开，均严格落在主192bit端点之间。共享SymPy1.14.0、Fraction及python-flint/Arb算术内核，如实披露，不声称不同区间库。

## 空间事件密度与权重

令Zi为独立尺度ci t的Cauchy变量，Xt=Σvi Zi。权已放进一维law中，因此变量变换用未加权列vi，不能再把ci乘进行列式。全部|vi|≤√n。设a=1/(2m√n)。在|x|>1/2处，不可能所有|Zi|≤a，因为那会给|Σvi Zi|≤1/2。因此至少有一个大参数i；这也覆盖所有cube face点及其邻域。

记S_lines=∪i span(vi)、d=dist(x,S_lines)>0，δ=d/(2m√n)。对该i，如果其他|Zj|≤δ，则距span(vi)最多(m−1)√nδ<d/2，矛盾。于是存在ordered pair i≠j，|Zi|≥a、|Zj|>δ。事件可重叠，以求和作上界；这里没有假设事件互斥或把事件概率相乘。

每pair用n−2条body列补成基，最大绝对det记dij。补列权均1。条件于其余变量后，两条被限制的一维密度分别由

\[
\frac{c_it}{\pi a^2},\qquad
\frac{c_jt}{\pi(\delta^2+c_j^2t^2)}
\]

控制；n−2条补列各由1/(πt)控制，最后除dij。剩余独立变量是总质量1的概率卷积，故不增上述L∞密度界。这是在空间密度层面做变量变换，而不是将路径事件的概率当成pt。乘积中的时间次数是t^(4−n)，权是ci cj；换成cmin的分母后得到

\[
 p_t(x)\le\frac{C\,t^{4-n}}{\delta^2+c_{min}^2t^2},\quad
 C=4m^2n\pi^{-n}\sum_{i\ne j}c_ic_j/d_{ij}.
\]

大参数i的尺度项被舍掉只使上界变大。n3 τ=3时cmin仍是1，不能误取3；所有本批cmin均严格正。

密度不等式可先按a.e.意义由测度分解推出。无条件body基product密度卷积给pt连续；|x|>1/2与d(x)>0都是开条件，且显示右侧随d连续，因此从邻域a.e.延拓至每个所述face点。不能只在零体积face上使用一个未经延拓的a.e.密度版本。

body基的最大det独立枚举为4(n3)、16(n4)，其余变量卷积不增sup，故pt≤1/(d*π^n t^n)。两种pt界分别积分t∈(0,1)、[1,∞)，得到正文的n3 log式和n4 atan式，尾常数D=1/[d*π^n(n−1)]。两个时间积分原函数均由独立SymPy实际验证。

## face投影与Voronoi分片

对face xk=σ、σ=±1，不平行于该face的方向线交点q=σv/vk，其坐标是0,±1且qk=σ。对face切向z=y−q，zk=0，精确有

\[
\operatorname{dist}(y,\operatorname{span}q)^2
=|z|^2-(q\cdot z)^2/|q|^2\ge |z|^2/|q|^2.
\]

最后一步由矩阵(|q|²−1)I−qT qTᵀ半正定；独立程序核了全部face、所有q的本征值。n3 mixed对应本征值(0,1)或(0,2)，n4 body对应(0,3,3)。平行face的线在xk=0中，距face至少1。因此在r=dist(y,Q)≤1时，所有方向线的distance均至少r/√n，δ≥r/(2mn)。

n3每face的Q有8点：4个body角点和4个face方向边中点；共6faces、48个帽。n4 purebody每face8点，共8faces、64个帽。交点集合在各face内先去重，不能因方向符号表示把同一条线重复算作不同q。

预算所需的r必须是最近交点distance。先作Voronoi分片，使每片中r=|y−q|，才可对各片用相应径向包络，最后把它扩大到以该q为中心的完整(n−1)维球。角点/边点的完整球越出face只扩大预算。不能直接在重叠帽中给任意中心套r=|y−q|；正文已经保留必要的Voronoi步骤。并列最近点可确定性分配，有限超平面ties对face测度为零。

## 幂界、帽积分及严格区间的含义

n3由log(1+K/r²)≤log(1+K)+2log(1/r)，得到G≤A+B log(1/r)。A,B均正。ε≤e⁻¹保证log(1/r)≥1，于是

\[
G^{3/2}\le\sqrt2\{A^{3/2}+B^{3/2}\log^2(1/r)\}.
\]

在二维完整帽球中积分2πrdr，原函数

\[
\int_0^\epsilon r\log^2(1/r)dr
=\frac{\epsilon^2}{2}\{L^2+L+1/2\},\quad L=\log(1/\epsilon),
\]

端点r=0的极限为0。得正文E3。n4由atan≤π/2得到G≤D+B/r，再用p=4/3的凸性幂界2^(1/3)，三维完整球积分4πr²dr，奇项积分∫r^(2/3)dr=(3/5)ε^(5/3)，得正文E4。所有指数、球面因子和1/2均已逐式核对。

这些上界趋零。方向线与face的交点只是有限零测集合，可在J的Lebesgue积分中剔除；这不意味着在生成元的全空间分布源中删去奇线。本审计没有因此断言正部构核合法。

独立常数汇总：

| n,τ | S=Σordered ci cj/dij | Cπ^n | ordered det计数 |
|---|---:|---:|---|
| 3,1/10 | 491/100 | 5892 | 4:36,2:30,3:24 |
| 3,1/3 | 92/9 | 36800/3 | 同上 |
| 3,1 | 32 | 38400 | 同上 |
| 3,3 | 156 | 187200 | 同上 |
| 4,purebody | 11/2 | 5632 | 8:32,16:24 |

十个单帽表达式的独立展示值如下；使用预算时须取JSON中准确有理上端，不能取展示小数。

| n,τ | E(10⁻⁶) | E(10⁻⁸) |
|---|---:|---:|
| 3,1/10 | .002416197224702365 | 4.197515342145587e−7 |
| 3,1/3 | .0001985785386497470 | 3.430562947295845e−8 |
| 3,1 | .00004133538952610277 | 7.097307019014634e−9 |
| 3,3 | .0004449220695322016 | 7.639334264709341e−8 |
| 4,purebody | .00009927848049515193 | 4.608098864843413e−8 |

`030_results.json`另给Σfaces #Q×E的准确区间：n3乘48，n4乘64。它是遗漏积分的上界表达式，不是实际帽内积分的双侧包围；E的下端不能当作实际积分下界。外帽区域的G幂积分以及g的全face下界仍须另证，不能凭本表给成本<8/16。

## 复现与停止

主导入后最小命令（默认stdout-only，可选--save）：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-finite-dim-bounds/audit/030_check.py
```

审计clone运行时加 `--project /Users/zhengzhihao/Desktop/T/Math-Research-finite-dim/cube-finite-dim-bounds`。入口先核全部七项来源SHA，任何漂移立即停止。MATLAB仅交主串行作十个帽径向包络直接积分和投影PSD的浮点cross；具体实际状态见 `030_matlab_receipt.json`，不替代本文证明或256bit包围。原主MATLAB补基回执作为锁源输入，未冒称本审计的独立执行。

未运行谱锥/Green值、成本优化、LP或任何新grid；本有限批以证明和常数通过、manifest归档后停止。后续031候选按新登记单独锁源审计。


主串行MATLAB已实际exit0（session31145），同锁源ff6b9af7bd5db58e00f3e1a78c939a2e2139c8cf8949d82fde19c7e9d13f103b。审计已读完整十个径向积分stdout，max relative error=3e−15，投影PSD n3/n4 assert通过；log SHA db044af659638ce2772fc56c47fa198db121fcd7d5f2127501cf7841dd6fd220。浮点值不冒称位于256bit窄包围内。收尾七项主来源哈希均未变，本批已停止。
