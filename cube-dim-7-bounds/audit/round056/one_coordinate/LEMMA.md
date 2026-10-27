# R056 一维范围扩展与inverse-face矩：独立局部引理

**局部结论ACCEPT，依以下guards使用。** 只审核稳定Psi/J闭式、baseline矩、a.e./AC匹配桥以及phi的inverse-face矩；不审核target计数、列系数改进或完整C_n结论。已锁正式源`upper/round056/refined_inverse/PROOF.md`实际SHA-256 `6986ed9ec48b46f3136a9a56722d0092db4f2942d52e35c9f4de37eb279f399a`，接入核对见第8节。

## 1. 必须区分的两个概率空间

* **baseline target law：** W~Unif[-r,r]，仅在R=r时有E_WJ(W;q)=1及下述logJ矩。
* **full source law：** Z~Unif[-ρ,ρ]，用于计算E[phi/ell]、E[phi²/ell]。随后乘capture/ell并除以J，才得到某固定captured face的free-coordinate law。

不能把baseline的W改成Unif[-ρ,ρ]，也不能把cap条件化、soft-weight条件化或跨face混合后的坐标继续当作原独立free law。

## 2. 扩大后的局部参数域与稳定闭式

取0<d≤1/3、0<h≤d/8、ε=ζt∈[0,d]，固定target r。baseline及扩展区间为

q_b=r−d+h，Q=[q_b−2h,q_b]。

整个Q上要求

q≥1/2，R=q+d−h≥1，ρ=tq≥4r，
a=(κ−1)ρ=(1−ε)q>0，
**g=R−a=d−h+ζρ，0<g≤min(d,R/2)。**

正号ζρ不可删除。g≤d需要例如ε·max(Q)≤h；固定d、h=2^(−24)/(Ld)、ζ=n^(−3)、t~n/L时最终成立。R∈[r−2h,r]。q_b可稍越source band；定义和guards须延拓到整个Q，而不是只检查原source q端点。

对候选d=3/10及target r∈[1.075,1.225]，扩展q总范围是[0.775−h,0.925+h]，而不是只到原source endpoints。h≤d/8时它仍有q>1/2且R≥1.075−2h≥1；ρ≥4r及ε·max(Q)≤h则由正式t、ζ参数在充分大L下另行满足。这些是局部引理的实际接入guards。

对|w|≤r，定义

ell_R(z)=min(1,(R+κρ−|z|)/(2R))，
J(w;q)=(2ρ)^(−1)∫_[-ρ,ρ]1_{|z−w|<κρ}/ell_R(z) dz。

令

F(R,s)=−2Rlog(1−s/(2R))−s=∫_0^s u/(2R−u)du，
m=(|w|−a)_+，s=(g−m)_+。

则精确公式为

**J=1+[F(R,g)+F(R,s)−m]/(2ρ)。**

与原A/B或Psi表达相同。a>0且m≤g+2h≤5d/4≪ρ保证只丢一侧source尾；m>g时F(R,s)=0并自动计入中央段删减，没有漏掉|w|>R的匹配端点情形。

稳定实现仍可用u=s/(2R)≤1/6时的正级数

F(R,s)=2RΣ_(k≥2)u^k/k，tail after K≤2R u^(K+1)/[(K+1)(1−u)]。

## 3. d≤1/3时的a.e.导数与全匹配shift

在0≤s≤g≤d、R≥1下，直接微分给

F_s=s/(2R−s)≤d/(2−d)≤(3/5)d≤d，
|F_R|≤2u²/(1−u)≤d²/(2−d)≤(3/5)d²≤d²，
0≤F≤s²/[2(2R−s)]≤(3/10)d²。

R'=1、g'=ε≤d，故F(R,g)的q导数绝对值≤2d²。

若|w|≤r−2d，则a≥r−2h−g≥r−5d/4，整个Q上m=0。于是J=1+F(R,g)/ρ≥1，

|∂_q logJ|≤[2d²+F/q]/ρ≤(13/5)d²/ρ<3d²/ρ。

边带中令B=F(R,g)+F(R,s)−m，J=1+B/(2ρ)。m>0时s=(R−|w|)_+，s'为0或1，−m'=1−ε，因此a.e.

|B'|≤1+d+3d²≤5/3<2。

m=0时|B'|≤4d²<2。另有

|B|/q≤(6/5)d²+(5/2)d≤29/30<1，
J≥1−m/(2ρ)≥1−5d/32≥91/96>1/2。

所以边带|∂_q logJ|≤3/ρ，足以保留原常数

**|∂_q logJ|≤(64/ρ)[d²+1_{|w|>r−2d}]。**

在|w|=a(q)处仍存在capture-onset kink，一般没有classical q导数；函数局部Lipschitz且在闭Q上绝对连续。|w|=R处因F_s(R,0)=0是一阶光滑拼接。故导数界必须按a.e./AC理解再积分，不能写成处处可微。

如果0<x=n/t≤L且q≥1/2，便有n/ρ≤2L。对K_w≤4nd的任一target（最多n个坐标），

Σ|∂_qlogJ_i|≤(64/ρ)(nd²+4nd)≤640Ld。

积分长度≤2h得到

**|ΣlogJ_i(q)−ΣlogJ_i(q_b)|≤1280Ldh=5/65536**

当h=2^(−24)/(Ld)。候选d=3/10时K_w≤n≤4nd自动成立，不需另证binomial tail。这个简化只适用于d≥1/4，不能写成所有d≤1/3都自动成立。x也必须保持正且满足所述上界；旧R045的x=(1−32d²)L在d=0.3时为负，绝不能原样移植。

## 4. baseline的3/6/12常数保留

仅此节取R=r、W~Unif[-R,R]。Fubini由“给定source z时接受的target W比例=ell_R(z)”严格给EJ=1。

记T=(|W|−(R−g))_+、B_0=F(R,g)，J−1=(B_0−A(T)/2)/ρ。条件g≤R/2给

B_0≤g²/(3R)，A(T)≤(4/3)T，ET²=g³/(3R)。

因此

E(J−1)²≤2g⁴/(9R²ρ²)+8g³/(27Rρ²)
≤(11/(27R))g³/ρ²≤3g³/ρ²。

同时|J−1|≤g/ρ≤1/8。于是|logJ|≤2|J−1|及logJ≥(J−1)−2(J−1)²，结合E(J−1)=0、Jensen，得到

**−6g³/ρ²≤ElogJ≤0，Var(logJ)≤12g³/ρ²。**

这段推导只需R≥1、g≤R/2、ρ≥4R，常数不依赖旧d≤1/16的小域。pinned target坐标W=±R仍有−g/ρ≤logJ≤0，为确定项。自由target坐标的独立性仅在固定full target face上使用；此处不把source/target两个law交换。

全匹配还有逐点界

**|J−1|≤5d/(8ρ)**，

因为正excess≤2F≤(3/5)d²而删减m≤5d/4。故任何固定0<epsJ<1，足够大ρ后J≥1−epsJ成立，可用于下一节；不需对n坐标作union bound。

## 5. 新phi矩：fullsource及captured free coordinate

令ell_pin=1−g/(2R)、v_pin=−logell_pin、phi(z)=−logell_R(z)≥0。两个source尾合计密度为du/ρ，其中ell=1−u/(2R)。直接换元ell得精确式

**E_full(phi/ell)=R v_pin²/ρ，**
**E_full(phi²/ell)=2R v_pin³/(3ρ)。**

中央ell=1的atom没有贡献。没有漏掉双尾因子2，也没有把第一式再乘一个2。

对target坐标w_i，captured inverse law的density相对于fullsource为1_capture/(ell J_i)。若J_i≥1−epsJ>0，非负性允许去掉indicator，因此

μ_i:=E_inv,i phi≤B_1:=R v_pin²/[ρ(1−epsJ)]，
E_inv,i phi²≤B_2:=2R v_pin³/[3ρ(1−epsJ)]。

这是**每个free coordinate**的条件law。可选的粗上界为B_1≤(9/25)d²/[Rρ(1−epsJ)]、B_2≤(18/125)d³/[R²ρ(1−epsJ)]，由v_pin≤u/(1−u)、u≤1/6得到；保留精确v_pin通常更好。

## 6. 固定face、混合face必须分开

在一个给定captured signed source face j上，pinned phi恰为v_pin，其他n−1个coordinate在上述inverse law下独立。故V=−logG有

E(V|face j)=v_pin+Σ_(i≠j)μ_i，
Var(V|face j)=Σ_(i≠j)Var_inv,i(phi_i)≤(n−1)B_2。

**跨face混合不能直接把右式称为Var_mix(V)。** 正确的total-variance还包含Var_j(E(V|j))。

这里可给紧的额外项：所有face使用同一组μ_i，且pinned v_pin与j及符号无关，所以face j均值=v_pin+Σ_i μ_i−μ_j。因为0≤μ_j≤B_1，任意inverse-row face混合权下

**Var_mix(V)≤(n−1)B_2+B_1²/4。**

该项一般不为零。更简单而且往往更合适的是逐face应用Chebyshev：令M=v_pin+(n−1)B_1，则每个face均值≤M。任意s>0时

P(V>M+s|face)≤(n−1)B_2/s²。

再对face作混合，**同一概率上界保持成立**，无需声称混合variance就是独立variance之和。此论证不要求不同face的均值相同。

以上均针对fullsource、capture及inverse weight下的law。进一步条件化原corner cap或soft b之后，free-coordinate独立性一般消失；应先在上述law估计再另付删除/retention，或重建新law，不能直接复用独立方差。

## 7. 实际新核验与范围

本机SymPy 1.14.0一次新检查验证两个phi矩的原函数及零gap边界；Fraction核对d=1/3端点的5/3、29/30、91/96导数常数、11/27<3、shift=5/65536。另符号核对face均值的total-variance correction等于被排除μ_j的variance，并给出该修正可正的简单代数例。未重算旧J baseline积分，未复跑旧生产、未用共享MATLAB。

```sh
PYTHONDONTWRITEBYTECODE=1 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-dim-7-bounds/audit/round056/one_coordinate/verify.py
```

结果见`verification.json`。这些检查支持公式及常数实现；全域正确性来自以上解析证明，不把有限算例当全域证据。本文件不提升为完整R056定理审计。

## 8. 正式源局部接入与身份锁

正式源为geometry独立clone的`upper/round056/refined_inverse/PROOF.md`，最终本次受审身份为`6986ed9ec48b46f3136a9a56722d0092db4f2942d52e35c9f4de37eb279f399a`。先前草稿`5d639f31cb7c62d05a09612d63fe22fc1af408cc5f311eed72d803c7e0c12eb3`至`623b7b15967e1437e45c641d171c84a50b5847558ccade39d7913f6be78028b8`只在本职责外修正target per-face cap/gate口径；最终版再补全band延拓及actual increment说明。实际对三身份逐字/子段哈希比较，确认§2局部J段和§5至mixing段完全不变。对应subsection hashes分别为`dcac1ea7b5276151b7380866a2b5e2b9a4267f5d347fa43d040f602f8d9832ae`、`c9b672ed708fdd675b09188f0a4c8f65027b640768da5229a281c305699ecc85`。本局部锁定不表示额外承担这些新全链段落的审核责任。

源实际选择x=93L/100>0，故n/ρ≤(6/5)L<2L；未移植旧的负x表达式。target interval已trim为[43/40+2^(−20),49/40−2^(−20)]，h≤2^(−20)使[q_b−2h,q_b]全部位于[31/40,37/40]。在该源中baseline q_b本身也在source band内；本引理允许的小外延不再必要，但仍覆盖其应用。

源的局部finite gates给h≤min(d/8,1/10000)、ζt≤h/2、q≥31/40、ρ≥4r。它们推出R≥10749/10000>1、g=d−h+ζtq≤d、κ>1，满足本引理全部guards。J逐点界方面源使用较宽的|J−1|≤2d/ρ，本引理已给5d/(8ρ)，故其统一J≥1−e（e=1/1024）接入安全。这些逻辑蕴含的核对不等于重审全部29个finite gates或整定理。

源§5使用本引理的**逐face Chebyshev后混合**路线，没有错误地把unconditional mixture variance写成独立方差和。其具体常数也新增以Fraction独立核对：R/q≤43/31，g≤3/10及R≥10749/10000给ell_pin≥3083/3583；令u=500/3083、Q=u−u²/2+u³/3，则

v_pin≤log(3583/3083)≤Q<1/6，
[43/(31(1−e))]Q²<4/125，
[(2/3)43/(31(1−e))]Q³<1/300。

于是每face均值≤Q+4x/125、方差≤x/300。χ=7/200，且χ−4/125−1/2000=1/400，给cutoff与均值间隔≥x/400，Chebyshev为1600/(3x)。在源L_0=2^20、x=93L/100处，Q/x≤1/2000及1600/(3x)≤e两式均以准确Fraction核验，随后随L增大只会改善。

锁源阶段只新增上述局部Fraction检查、源哈希及subsection字节比较；原SymPy运行未重跑，也未调用作者生产/guard脚本。完整源码的其他证明责任仍交由整链审计，当前ACCEPT严格限于此局部接口。
