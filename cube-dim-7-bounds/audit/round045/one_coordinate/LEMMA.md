# R045 一维闭式及uniform-q乘积桥：独立局部引理

**局部结论ACCEPT。** 本文件只审核一维J、其baseline矩和q导数到typical-target乘积的桥，不审核R045完整新C_n定理、soft-cap payment或后续有限a/全算子链条。

对照源：`upper/round045/soft_inverse_column/PROOF.md`，实际SHA-256 `add435fdf2f55993eba9064bd50c493e4c42769c1c48463329c4015be2e66ab9`。源§4及§5截至乘积式(11)与以下更弱假设下的独立证明相容；特别已正确使用locally Lipschitz、a.e.导数和绝对连续性穿过kink。

## 1. 参数与精确闭式

固定0<d≤1/16、0<h≤d/8、0≤ε≤d、t>0及target半径r≥1/2。令

q_0=r−d+h，Q=[q_0−2h,q_0]，
R=q+d−h，ρ=tq，a=(1−ε)q，g=R−a。

假设整个Q上q,R≥1/2、ρ≥4R、0<g≤d。于是R∈[r−2h,r]、a>0。这里ε对应ζt，不随q变化。实际匹配条件|r−q−d|<h给q∈(q_0−2h,q_0)，并不要求q_0本身属于原source band；只需上述闭式及微分定义在这个小扩展区间成立。

对|w|≤r，设capture半径为ρ+a、

ell_R(z)=min(1,(R+ρ+a−|z|)/(2R))，
J(w;q)=(2ρ)^(−1)∫_[-ρ,ρ] 1_{|z−w|<ρ+a}/ell_R(z) dz。

定义非负稳定excess

F(R,s)=−2R log(1−s/(2R))−s=∫_0^s u/(2R−u)du，0≤s≤g。

令m_±=(±w−a)_+，E(m)=F(R,(g−m)_+)。则

**J=1+[E(m_+)+E(m_−)−m_+−m_−]/(2ρ)。**

证明：ell在两端长度g的source尾上逐线性减小，中央为1；a>0使capture至多从一个source端点删去长度m=(|w|−a)_+。因m≤r−a≤g+2h≤5d/4≪ρ，该删减不会穿到另一source尾。保留的普通长度为2ρ−m，两个尾的剩余inverse excess为F(R,g)及F(R,(g−m)_+)，即上述公式。边界开闭只影响零测集。

这与作者的A/B表达完全等价：B=F(R,g)，A(R,g,T)=T+F(R,g)−F(R,(g−T)_+)，故2B−A=B+F(R,(g−T)_+)−T。

数值稳定实现可不用相近量相减：令u=s/(2R)≤d≤1/16，则

F=2R Σ_(k≥2) u^k/k，
0≤F−2RΣ_(k=2)^K u^k/k≤2R u^(K+1)/[(K+1)(1−u)]。

## 2. q导数：明确的保守常数及kink

以下偏导固定另一个参数：

F_s=s/(2R−s)，
F_R=−2log(1−u)−2u/(1−u)。

由非负幂级数可得，在0≤s≤g≤d、R≥1/2时

0≤F≤d²，0≤F_s≤2d，|F_R|≤4d²。

又R'=1、g'=ε，故F(R,g)的q导数绝对值≤6d²。

作者A/B版本的偏导常数也由同一计算覆盖。对T≤g，A=∫_0^T f(y)dy，其中f=2R/(2R−g+y)，f_R=2(y−g)/(2R−g+y)²、f_g=2R/(2R−g+y)²。故|f_R|≤8d、f_g≤4、f≤2，积分给|A_R|≤8d²、|A_g|≤4d、A_T≤2。对T≥g，A=F(R,g)+T，直接有|A_R|≤4d²、|A_g|≤2d、A_T=1；两支一阶偏导在T=g吻合。B=F(R,g)满足|B_R|≤4d²、|B_g|≤2d。因此作者采用的16/8/2及8/4等更宽常数安全。

若|w|≤r−2d，则a≥r−2h−g≥r−5d/4>r−2d，整个Q上m=0。所以J=1+F(R,g)/ρ≥1，且

|∂_q log J|≤[6d²+F/q]/ρ≤8d²/ρ。

一般边带|w|>r−2d中令B=F(R,g)+F(R,(g−m)_+)−m，J=1+B/(2ρ)。在m>0时s=(R−|w|)_+，其q导数为0或1；m'=−(1−ε)。因此a.e.

|B'|≤6d²+4d²+2d+1<2。

m=0分支则|B'|≤12d²<2。同时

|B|≤2d²+5d/4≤11d/8，|B|/q≤11d/4<1，
J≥1−m/(2ρ)≥1/2。

于是边带a.e.有|∂_q logJ|≤3/ρ。综合给更宽松的作者目标

**|∂_q logJ(w;q)|≤(64/ρ)[d²+1_{|w|>r−2d}]。**

不能声称处处classically可微。对固定w>0，在w=a(q)处左、右J导数之差恰为

(1−ε)/[2ρ(1−g/(2R))]>0。

该点是真capture-onset kink。相反|w|=R处F_s(R,0)=0，两支一阶导吻合。J在有限多个分支上光滑且所有one-sided slopes有统一界，所以在闭Q上绝对连续；logJ因J≥1/2也绝对连续。上述a.e.界可合法积分，不能把“有kink”当成桥接失败，也不能忽略它而宣称处处导数。

## 3. Baseline R=r的均值、log均值及方差

在q=q_0时R=r。令W~Unif[-r,r]。对固定z，capture在target W中的比例正是ell_r(z)，所以Fubini严格给

**E_W J(W;q_0)=1。**

这不是按source取均值，也不适用于R≠r时仍按W~Unif[-r,r]取平均。

闭式还可直接校验这一归一化。若F_0=F(r,g)，则

∫_0^g F(r,s)ds=g²/2−(2r−g)F_0。

两端在g=0均为0、对g求导均为F_0。把它代入target中央段和尾段的J公式，恰得E(J−1)=0。

中央|W|≤r−g时|J−1|≤g²/ρ；尾部概率g/r≤2g，且尾部|J−1|≤g/ρ。因此

E(J−1)²≤(g⁴+2g³)/ρ²≤3g³/ρ²。

全域J≥1/2使|logJ|≤2|J−1|，且logJ≥(J−1)−2(J−1)²。结合Jensen，得到

**−6g³/ρ²≤E logJ≤0，Var(logJ)≤12g³/ρ²。**

作者只在R≥1使用这些常数，强于此处R≥1/2假设，无须修改常数。这里的log均值控制不可省略：仅有小variance并不足以证明乘积接近1。

target固定face坐标W=±r时，J=1+[F_0−g]/(2ρ)。因0≤F_0≤g/2，明确有

**−g/ρ≤logJ(±r;q_0)≤0。**

所以pinned坐标可作为趋零的确定项保留，不能错误地赋予它另一条Uniform law。

## 4. Typical target及uniform匹配桥

设N=n−1个free target coordinates独立Uniform[-r,r]，另一个坐标固定为±r。记ρ_0=tq_0，b_0=N d³/ρ_0²。全坐标log和S_0满足

|ES_0|≤6b_0+d/ρ_0，Var S_0≤12b_0。

对任意v>6b_0+d/ρ_0，Chebyshev给

P(|S_0|>v)≤12b_0/[v−6b_0−d/ρ_0]²。

实际target r≥1时，K=1+Bin(N,2d/r)。若nd≥1，则EK≤3nd、Var K≤2nd，从而P(K>4nd)≤2/(nd)。两个good事件可以取交，完全不要求它们独立。固定d且n→∞，q_0≥1/2、t=n/x、x≤L给b_0≤4d³L²/n→0及d/ρ_0≤2dL/n→0。因此对每个r，|S_0|≤1/8且K≤4nd的target比例为1−o(1)，误差常数在该r区间内一致。这里不声称一个事件同时对所有r成立；原row需要的只是给定target r的统一估计。

现在固定任意good target w。对所有q∈Q，a.e.

|∂_q Σ_i logJ(w_i;q)|
≤(64/ρ)(nd²+4nd)
≤320nd/ρ≤640Ld，

使用d≤1及n/ρ=x/q≤2L。由绝对连续性，任意matched q满足

**|ΣlogJ(w_i;q)−ΣlogJ(w_i;q_0)|≤1280Ldh。**

此为全部连续q的确定性结论，不是有限q点或网格抽样外推。若h=2^(−20)/(Ld)，右端恰为5/4096。与|S_0|≤1/8合并后，所有matched q都有|ΣlogJ|≤517/4096<1/4，故作者所需

**1/2≤∏_iJ(w_i;q)≤2**

成立。q_0可稍出source band，但整条bridge区间仍满足引理假设；这正是需要单独延拓定义的用途。

## 5. 与正式源及工具检查的对接范围

作者参数为ε=ζt、r=1+u、u∈[d/4,3d/4]，实际d≤2^(−24)（主拟选2^(−30)），h=2^(−20)/(Ld)。固定d后，充分大L使h≤d/8、ε≤d及εq_0≤h，且所有扩展q仍有q,R≥1/2、ρ≥4R、g≤d。因此本引理覆盖源§4–5所用参数，包括q_0越过source-band端点的小延拓。未给出整个R045定理的共同显式维数阈值。

本机SymPy 1.14.0实际核验F偏导、baseline积分primitive、onset kink的精确jump、|w|=R处C1拼接。512bit Arb在一组新的有理局部参数d=1/32、h=1/512、ε=1/2048、r=65/64、q=63/64、ρ=126，对两个deep位置和两个edge位置核验正性与64导数界，并以非负级数和严格有理tail核验F闭式。这些点处于本局部引理范围，用于检查实现；不模拟主定理的高维参数，不替代全域证明。

```sh
PYTHONDONTWRITEBYTECODE=1 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-dim-7-bounds/audit/round045/one_coordinate/verify.py cube-dim-7-bounds/upper/round045/soft_inverse_column/PROOF.md
```

结果见`verification.json`，其中冻结了§4和§5截至乘积桥的独立subsection hashes。未复跑Terra的4365次积分、未用共享MATLAB、未改他人文件或运行生产计算。结论不超出此一维/匹配桥职责。
