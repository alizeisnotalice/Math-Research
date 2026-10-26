# B086 独立审计：rank≥4 aggregate 与统一物理 Hessian 矩阵

**结论：A086 完整解析引理及其精确有理接口通过本审计，建议主任务接受。** 未发现固定基、绝对值正和、次数、mixed multiplicity、PSD拼接、零/平行列或物理列权重方面的缺口。主任务保留最终接受权。本轮结果仅保存在本地，不进行远端操作。

通过范围是：对固定cutoff的rank≥4部分，以独立四列子集替代逐pattern常数求和，并给出统一的全空间C2及双侧Hessian矩阵界。该矩阵仅支付 `R_(≥4,T)`；**不是完整Green的Γ，不是已实例化的body32常数，更不是g、I或Cn结论。** 093的完整账目与最终成本S不属于本审计，未提前背书。

## 1. 封源和独立验证范围

锁定作者11文件（含manifest）与076两项上游来源，共13项，散列逐项吻合。

| 来源 | SHA256 |
|---|---|
| A086 manifest | `d99aa0c047803262155d37f4537d27bf01d70f95bc684028f3cf5afe57dea238` |
| theorem.md | `309c973bcbcbba206293d9ae47ae2a478c7dbb808c94e9e44510ac0a6d75d581` |
| aggregate.py | `bdc3772e5d8f28c0ace2928aeab3f74885fe992ea63ca847f837dc3974403120` |

完整路径见`source_bindings.json`。作者manifest明确排除的历史WIP不是最终证明或运行状态；本审计不将该旧零配置状态混入最终账目。全文核读theorem/report、静态审查aggregate/checker，并核查输入与执行来源。没有运行作者四个6×8 toy配置，没有复跑076或085。

两项全新轻量接口配置已在单次预登记程序中通过：一个新的正有理加权非轴六列frame，以及在该frame中插入零列。验证使用本地Python标准库Fraction，参照矩阵逆由独立的有限nilpotent级数给出，未使用作者Gaussian elimination生成参考答案。有限配置支持接口实现；全称数学结论仍由下列完整解析审查承担。

沿用主已授权的独立`.git`局部工作副本，完整checkout曾失败的既有状态不作虚假声明。本轮未再clone，仅写独占`audit/round086/rank4_aggregate_review/`。没有stage、commit、push、MATLAB或递归代理。

## 2. 固定选择和正和因式分解

**固定J与completion的论证成立。** 每个rank(A_L)≥4的L含某个实际独立四列子集，因此字典序首个j(L)存在且唯一。每个独立J又可在完整spanning frame中扩充成六列基；余下两列由固定字典序选定。N的列序明确为J在前、K在后。选择只依赖已固定的物理列阵和标签，不依赖x、t、方向或具体assigned pattern。

在同一个J类中，所有pattern使用同一仿射解 `w_N=N_J^−1(x−A_(N^c)z)`，因此空间微分只命中六个已解标量因子。未解因子虽然影响已解参数的值，但其自身仍是z的函数，对x的导数为0。其积分域不动。选中因子及前两导数在固定t下全局有界，未解因子由可积的Cauchy乘积支配，足以证明逐pattern的C2与固定积分微分。

**绝对值的顺序正确。** 作者先展开有序Leibniz命中索引并取每项绝对值，才扩大assigned集合为全部L⊃J。同J的方向系数与L无关，此时所有被扩大项均非负，故式(8)合法。J⊂L也保证扩大的pattern仍属于rank≥4；但同pattern可以被不同J majorant过计数，这不改变原始唯一分配。

扩大后才逐标量因式分解。四个J因子一直是large因子；两个completion保留 `S_r=|(cχ)^(r)|+|(cη)^(r)|`；未解因子未受微分，故其两种选择精确和为c_t。未解Cauchy乘积积分为1，包括m=6时的空乘积/零维积分。由此消去的是未解概率质量，而不是signed导数中的假抵消。

没有把rank条件变为依赖x的事件再微分，也没有在扩大后替换N或逐方向寻找更好的basis。零列不能进入独立J或可逆completion，只能留在未解变量；其χ/η和仍为c_t，积分为1。重复和平行列同样由实际rank测试正确排除非法基，保留其原始独立随机变量标签，不会破坏分解。

## 3. cutoff绝对导数、时间幂及k阶常数

**四条尾的导数支撑不丢失。** 对ℓ=cη，在|s|<R内函数与前两导数都为0，边界由C2匹配覆盖。在|s|≥R，Cauchy尾界和cutoff导数界给

    π(A0,A1,A2)=(R^−2,4R^−3,20R^−4),   |ℓ^(r)|≤Ar t.

常数4为2+2；常数20为6+2·2·2+6。每个被命中0、1、2次的tail因子仍保留一个t，而不是把tail条件在微分后丢掉。

completion必须同时支付χ和η两侧的绝对导数。作者式(5)的

    S1≤|c'|+2c|χ'|,
    S2≤|c''|+4|c'||χ'|+2c|χ''|

准确包含所有Leibniz成本。利用0<t≤1得到 `π(C0,C1,C2)=(1,1+4/R,2+8/R+12/R²)`。S0=c_t是真正恒等式，所以不应额外乘2；S1、S2只是不等式，没有假定它们等于|c^(r)|。

一次k阶项若q次命中四条J尾，其时间幂为 `t^4 t^(−2−(k−q))=t^(2−k+q)`。在t≤1时可粗化为t^(2−k)。这也覆盖二次命中同一tail，q=2而tail本身仍只贡献t；或二次命中同一completion，正确产生t^−3的那个因子。

以A0⁴C0²归一化得到 `q_J=1/(π⁶R⁸|detN_J|)`。二阶可能比值只有a2、b2、a1²、a1b1、b1²，作者h2全部覆盖。有序微分命中索引的方向系数绝对和，正好为各`||N_J^−1 u_r||1`的乘积；因此mixed二阶项的2和多重指标重数没有遗漏。式(12)–(13)分别给value、gradient、Hessian的t²、t、1支配，均可在0端积分。

支配常数对所有x统一、独立J有限，故时间积分不仅具有分布意义Hessian界，而是全空间C2，且可逐次积分其前两空间导数。该结论对固定rank≥4部分成立，不需要额外S2几何条件；在076应用中仍须使用同一个原模型和同一个cutoff。

## 4. 统一PSD矩阵及15项平方余量

对一个固定方向u，令z=N_J^−1u。保留实际时间幂，tail-tail、tail-completion、completion-completion项分别带t²、t、1。积分后各非对角对的系数为：

| 类型 | 对数 | 每对绝对乘积系数 |
|---|---:|---:|
| tail-tail | 6 | `2a1²T³/3` |
| tail-completion | 8 | `a1b1T²` |
| completion-completion | 1 | `2b1²T` |

用 `2|z_i z_j|≤z_i²+z_j²`，每个tail参与3个tail对、2个混合对，每个completion参与4个混合对、1个completion对。因此精确得到

    d_tail=(a2+3a1²)T³/3+a1b1T²,
    d_comp=(b2+b1²)T+2a1b1T².

作者式(18)的平方余量依次为6项系数a1²T³/3、8项系数a1b1T²/2和1项系数b1²T，恰好是全部15对，没有少付mixed项。

由此 `M4=Σ_J q_J N_J^−T diag(d_tail I4,d_comp I2)N_J^−1` 是在所有方向之前定义的一个PSD矩阵。对任意u的绝对二次型界都由同一M4给出，所以得到 `−M4≼D²R_(≥4,T)≼M4`。不能用不同方向各自选择的completion最小值拼一个矩阵；该证明没有这么做。若以后更换completion，必须整套固定并重新绑定数据。

trace支付额和物理拉回也正确：PSD矩阵的最大特征值不超过trace；有理A、R、T使π⁶M4、π⁶traceM4有理。face参数化E必须用EᵀM4E。这个trace仍只属于rank≥4部分，不包括rank3的cutoff-removal、Taylor、signed层或full-G大时尾。

## 5. 精确接口的代码边界

`aggregate.py` 的rank与determinant/inverse均为Fraction Gaussian elimination；row交换符号、pivot乘积和inverse同步行操作一致。`independent_four_bases`枚举所有四列tuple，rank=4才进入completion；completion顺序固定，输出N按J+K排列。最高层`certificate`调用完整枚举后生成矩阵，因而关闭J覆盖义务，不枚举scalar patterns。

`majorant_from_bases`的矩阵项使用 `inverse[k][i]*inverse[k][j]`，正是物理congruence中第k行的外积；前4行/后2行权重与实际N列序一致。共同π^−6显式留在返回值之外，R^−8与|det|^−1均保留。K012使用的Λ是inverse全部元素绝对值之和，虽粗但确为2→1范数上界。

**低层函数有明确调用前提。** 它不会自行证明传入records覆盖全部J，也不认证外部伪造的inverse/列序。有效调用必须使用`certificate`，或提供已独立认证、完整且无遗漏的records。这是作者已写出的接口义务，不是本审计替它补出的结论。绝不能把部分records产生的数字称为全aggregate。

正谱权必须先吸收入实际物理列。算法不擅自归一化列长，det与inverse都使用该实际列阵；因此正有理权重变化被完整保留。若调用者先把不同长度列统一成方向，将得到另一个模型的证书，不能转用于原模型。零列可被接受的范围仅为本文rank≥4引理，不据此擅自改写076其他部分的假设。

## 6. 两个新配置及实际计算证据

配置1使用一个新的upper-triangular shear H（详见登记），物理列阵 `A=H diag(1/2,2,3/2,5/3,7/4,3/5)`，detA=21/8，R=2/3、T=3/8。H=I+U且U⁶=0，参考inverse独立用 `H^−1=Σ_(k=0)^5(−U)^k` 构造；没有复用被审Gaussian elimination求参照。

因为这六列本身是一个基，全部15个四子集独立，各completion唯一。每个原始物理列在其中10次作为tail、5次作为completion。因此可不枚举scalar patterns，独立得到闭式

    π⁶M4 = [(10d_tail+5d_comp)/(R⁸|detA|)] A^−T A^−1,
    d_tail=4401/512,   d_comp=729/16.

接口的15条J/K、det、inverse的J+K行序、整个6×6矩阵、trace和K012全部与该参照严格一致。配置1覆盖m=6空未解积分的边界；它是新的非轴加权代数frame，不是旧axis核或body32实例。

配置2在该同一frame的零基index2处插入一个零列。仍只有15个独立四子集，零列不进入任何J/K，矩阵和K012与配置1完全一致。这核对了“零列留在未解概率质量中”的具体接口性质。两配置共生成并核对30个新toy basis records，未枚举任何scalar patterns，也未计算密度或Green积分。

程序一次完成两个配置，tool chunk `e737b7`，guard PID2217、child PID2218，exit0。全guard墙钟0.071527417秒，guard CPU0.045773秒、child OS CPU0.063425秒，总0.109198秒；无活进程。登记允许165秒工作、3秒kill/reap，179秒硬alarm包含最终收据，满足总墙钟180秒；child CPU27秒加guard CPU2秒低于CPU30秒上限。无科学失败、超时或自动重发。

执行只导入散列固定的`aggregate.py`，没有运行作者checker；`-B`防止在作者目录写bytecode。作者四项旧toy、076/085、body/native kernel均零调用。矩阵结果仅属于上述新小frame，不能称作真实body32 Γ或最终成本S。

## 7. 最终状态和剩余账目

本轮建议接受rank≥4 aggregate定理与其最高层完整证书接口。全称解析审计完整，有限新配置作为实现证据，不替代证明。未发现必须作者补证的缺口。

下一实际模型计算仍需完整物理frame、认证cutoff/time、完整J records和同一归一化；还须逐项支付rank3所有余项和原full-G上尾，匹配真实标量顶点，并比较完整误差/成本。此处没有将“存在有限常数”升级成实用Γ，也没有对093尚未审核的最终账目作结论。

最小本地交付为本报告、两个新配置checker、单次guard、登记、来源绑定、结果、执行/工具收据及manifest。single-use guard拒绝删除历史后隐式复跑。所有新文件均留在本地独占目录，无stage/commit/push。
