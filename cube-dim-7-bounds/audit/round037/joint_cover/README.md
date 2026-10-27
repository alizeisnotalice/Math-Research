# R037 四立方体负例与独立 splitting 接口审计

**ACCEPT最终稿** `upper/round037/joint_cover/SPLITTING_PROOF.md`，SHA-256：`5c74827b95956d881a7ef7f5202c11812c1bafa5fab9f96fdea64a317e79dfef`。

精确范围为：有限正体积开轴平行盒、λ>0、同一个原L1输入、固定正预算0<h_i≤M_i/λ、非负系数及守恒支付；所有需要选anchor或求平均的集合均非空。四cube负例对每个n≥2成立。它反驳对所有有限实例推广的half-gap断言，**不凭自身反驳注册有限batch的H或提供该batch统计**。splitting是合法的局部支付机制，尚不是新C_n界。

初稿SHA为`da365936198f2067e8820ee99787bde529ed433e9602e34b58e5f580d1aed17b`。独审指出的预算正性和“registered H”范围两处文字已修复；本报告接受的是上述最终身份。未重跑48case、旧toy或生产LP。

## 1. 四cube几何与真实partial-flow

中心±4e_1、±4e_2且原半径1时，任意两个Q的某坐标中心差至少4>2，故小cube两两严格分离。witness半径4时，opposite两对中心差8，开盒不相交；所有cross两点的各坐标差≤4<8，公共交集为正体积。图正是C4=K_(2,2)，无仅凭边界接触制造的边。额外n−2坐标区间完全相同，以下二维几何及体积比对所有n≥2原样成立。

所有被选cube的increment均为v=2^n，h_i=H=8^n，p=v/H=4^(−n)。tie顺序0,2,3,1，b=7、τ=3p/2时，四步eligibility load除以p依次为0,1,1,2。前三个可选，最后1不可选且单调性确保停止，故D=Q_1。这是真实legal selection，不是任意给定的q。

E_j=8Q_j：E_2和E_3包含整个D；E_0在第一坐标的下端点为−4，故恰覆盖D的一半。两类D atom的相对体积为1/2、1/2，成员分别为{0,2,3}、{2,3}。合法比例流在三成员atom各付1/3、两成员atom各付1/2；每份均≤p/τ=2/3，原first-rejection容量得到满足。因此c_0=v/6、c_2=c_3=5v/12，真实初始q/p=(7/6,0,17/12,17/12)，支付总量4v、实际负载31p/12，均正确。

## 2. 真最优mass与所有严格saving concentration停滞

C4的四条极大边各有正体积精确incidence cell，故每边放H/2是可实现的mass primal：每顶点收到H，总质量2H。dual z_2=z_3=1、其余零时每条极大边dual和为1，目标2H。原始与对偶等值严格证明M(V)=2H，无数值solver误差。

Lopt=|U|/M=2p也可直接在payment问题核对：四边load的平均值是(Σq_i)/2=2p，所以任何非负守恒q的max-edge load至少2p。把q_0均分给nonedge{0,1}后，q'/p=(7/12,7/12,17/12,17/12)，每边恰为2p，因而确实最优。原证书gap=7p/12，占L0的7/31>10^(−4)。

所有h_i相等，任意非空clique含正支付时，其anchor预算H，重分配后系数等于原clique系数和，s_C恒为零；全零clique也有s_C=0。故明确要求s_C>0的Algorithm D在**任何候选池**上都不能接受一步。这是严格saving规则的结构反例，不是在断言所有s_C=0移动皆无用。uniform输入f=(1+ε)1_(∪V_i)、λ=1让每个M_i=(1+ε)H，确保全部支付均对应同一个真实原输入。

## 3. Independent-module及一般boundary bill

F非空、独立且h_i>0时，t_F=Σh_iq_i/Σh_i是加权平均，因此非负且≤q_max；替换为公共t_F精确守恒。每个团最多含一个F顶点。如果每个外部正系数顶点对F全邻或全不邻，含任一F点的负载只能再付共同外邻团；任何共同外邻团又能与取得q_max的F点合并。故旧峰=max(L_rest,q_max+L_good)，新峰=max(L_rest,t_F+L_good)，严格改进条件与正文一致。

外部零系数在当前一步不变，可以忽略；F中的旧零系数可能被激活，所以module条件必须对**全F**检验。嵌入LP或后续步骤允许外部零变正时，需使用完整allowed-positive图、固定被删零系数，或每步重新核查。正文已包含这个必要范围。

任意独立F时，在一个F witness内部，把活跃外部项分成共同外邻与该顶点的B_i；它们各自是团，分别≤L_good和Δ_F。不在任何F内部时≤L_rest。再用q_max+L_good≤L(q)，即可得到命题3两层上界。这里无需假设两组能共同组成团；保守性只会加大上界。s_F≥Δ_F是充分非增条件，并非必要条件。尚无所有几何族的累计账单估计。

## 4. Nonedge pair精确最小化

对nonedge i,j，完整极大团列表可按含i、含j、两者皆不含三类划分。外部系数固定后，其最大负载恰为Φ(x)=max(A_rest,A_i+x,A_j+(P−h_ix)/h_j)，没有被遗漏的同时含i,j分支。每个顶点至少属于某个极大团，孤立点给相应A为零；最后一类为空可用零，因为q'≥0。

两条非恒定直线分别严格增、严格减，交点由h_jA_i+h_jx=h_jA_j+P−h_ix解得x=(P+h_j(A_j−A_i))/(h_i+h_j)。交点在区间内时最小化二者最大值；区间外取最近端点；加上常数A_rest不改变该点为最小点的事实，尽管可能使最优点不唯一。故命题4的clip公式及非增结论正确。P=0时区间退为{0}也成立。当A_i=A_j时退化为两点等系数splitting。

该结论只给固定外部系数的坐标最优，不给全局payment最优或重复迭代的有限终止结论。

## 5. 最小实际核验及复现

本机Python Fraction实际核对：四cube前两坐标严格交图、原Q分离、指定顺序eligibility数值、两atom内点的E成员、比例流容量、守恒q、C4四边load及显式mass primal/dual等值；并检验pair交点在内部、左侧、右侧三个代数情形。后两情形包含外部恒定branch主导的退化。只用有理运算，未重做高维atom化、参数搜索或LP。

```sh
PYTHONDONTWRITEBYTECODE=1 python3 cube-dim-7-bounds/audit/round037/joint_cover/verify.py cube-dim-7-bounds/upper/round037/joint_cover/SPLITTING_PROOF.md
```

精简结果见`verification.json`。二维有理计算只是实现/代数辅助；任意n结论来自上述被动坐标和解析体积比论证。没有生产batch统计、旧算例复跑、MATLAB、n=7特化或Git修改。
