# 第39轮 channel 稿独审

结论：本次只读版本的§2–5及§7末段真实cube非恒定性论证通过，未发现实质错误。新异表结论是共同η、固定n、固定源count分割和同一固定L下的加权χ²/L1结论；不控制跨尺度max。未重跑MC，未改channel目录或冻结源码。

审稿对象为 `../round039_channel/PROOF.md`（2026-10-29落盘版本，SHA256 `6fefe580ec45be62cc5e0b1d00623c50cabe35f8ba0a58fce71f5d0ae78b5ba6`），并只读核对 `probe.py`、`trace_checks.py/results`、`multitable_checks.py/results` 的公式/证据口径。下面是独立检查链。

## 1. Count Gram及异表推广

标准化binary变量Z_i均值0、二阶矩1。不同子集乘积的内积为0，同一子集为1，故对称化后的e_k有正确的binom(n,k)^{-1/2}归一化；其只依赖S且n+1维，因此是全部count函数的基。

给定完整query后各latent C_i条件独立，v_k中的每坐标g_i最多出现两次。计算E(v_jv_k)时，子集不相同必有只出现一次的g_i，均值0；相同子集时恰得到f^k。没有隐含高阶g矩，也不要求query law本身仅由f决定。将singleton indicator展开时系数为π_s e_k(s)，除π_s后正是文中L_s展开，故完整weighted Gram归一化正确。

对集合E，variance为Σ_{k≥1} a_{E,k}² f^k/π_E²。这说明每个bin的diagonal随f单调；非对角Gram元素可能为负，稿子没有偷用其单调性。表b的a系数只依赖共同Bin(n,η)，因此各表Q_b不同也不妨碍逐项用f_b≤f*比较。

异表加权和中，对固定k，bin内Cauchy给Σ_b a_{E_b,k}²/π_b≤E e_k²=1，前提正是这些E_b是同一count空间的固定分割。于是加权χ²预算≤Σ_{k=1}^n(f*)^k。这里不需要将不同Q_b强行放成共同query律，亦无R因子。空段保留为零项，商式排除π_b=0，方向正确。

坐标通道不同但共同η的补充也正确：每个k的谱值替换为第k个f_i的归一化初等对称多项式，仍≤(max_i f_i)^k。若η_b或η_i改变，共同count基与同一π不再成立；稿子明确未作这种推广。

## 2. 固定尺度L1与真实来源

在共同基准μ下，各A_b是各自Q_b^n的单位质量密度。三角不等式后分别在各自Q_b^n下积分，再作两次scalar Cauchy，确实给

`∫|ΣA_b(T_b−π_b)|dμ ≤ Σπ_b sqrt(Var_b) ≤ sqrt(Σπ_b Var_b)`。

这个跨表步骤只用Σπ_b=1，不要求A_b相同或通道间独立。G质量为Σ E_b T_b=Σπ_b=1，prior也质量1，没有漏乘/多乘π。

真实cube channel的one-coordinate query密度是a_b=A_Lρ_b，posterior为A_L(ρ_bχ)/a_b。周期Fubini给∫a_b=1和E_{a_b}p_b=∫ρ_bχ=η。源坐标产品与各坐标独立cube位移构成独立复制通道；Bayes精确给产品cube均值A_b乘PB段后验。因而文中G、P分别是真gated和实际prior的同一个L平均，归一化正确。

第38轮共同D和每格峰质量相同，使所有表的cell law相同，χ为前127格时η共同且一般不等1/2。表内位置旋转不破坏这个性质。若另换成table-dependent背景或改变cell质量，必须重新查η；当前证明不自动覆盖该改动。

## 3. Macro signal帽与bin lift

L≥1、m=floor L的窗口可拆成m个完整周期加不足一周期的余弧。plus/minus的余弧质量上界分别η与1−η，得到posterior的两个端点。由(p−a)(b−p)≥0的期望导出的variance帽准确，代数简化为

`f≤u/[m(m+1)+u]`，`f/(1−f)≤u/[m(m+1)]`。

这给固定L的L1上帽≤1/[2sqrt(m(m+1))]，m1时1/sqrt8。整数L的p恒η、f0，确应更强地完全相同。该corollary不依赖coarsened数值f，也不要求原匹配下界参数资格。

bin-conditioned posterior确为pbar=E_Q(p|bin)，其中bin query概率须按a加权，而不是uniform bins。稿中m_±为uniform-bin核均值时，pbar=m_+/(m_++m_−)，与上述加权条件期望一致。给出的p−pbar代数恒等式正确。两个class振幅各≤B/L、a≥1/L时，绝对差≤B；全方差给f_cont−f_bin≤B²/u。B为明确解析range帽时该迁移成立。

本源N_transition≤2的range预算确给B=29879/44945280：常数项1/97920、峰项1/1530、D项110/(9·255·48960)相加。η界使用245/81作为Σ_{m≤3}20/(9m²)，系数1/1020也正确。稿子已明确，解析差额不能认证已保存float f_bin的值；约.01586仅为诊断，未偷换成连续signal上界。

## 4. §7末段的真实cube非恒定性

对严格正的周期ρ和proper非空单弧χ，h=(χ−η)ρ在弧内为正、弧外为负（均按a.e.解释）。若p恒η，A_L h恒0。L1函数的窗口积分绝对连续，其a.e.导数为两端h之差，故h对L mod1的平移不变。这迫使正号弧的indicator也对该平移不变。

proper单弧不可能被非零圆周平移保持到零测差；例如其第一个Fourier系数非零，平移不变即要求exp(2πiL)=1。因此L非整数时矛盾，p非恒定，f>0。严格正性和proper弧假设不可删，稿中均已保留。固定一维来源再让n变化的声明也没有冒充原匹配n依赖来源。

Chebyshev packing中，给定完整labels后U的条件均值仅依赖S，故再条件到S时不产生额外均值方差；间距和半径预算一致。它是posterior max的L1反例，与所有表相同时G=P的完全抵消不矛盾。

## 5. 结论边界

本审范围内不需要修正的公式或量词。固定L的∫|gated−prior|小，不意味着∫sup_L差小，也不意味着两者weak maximal排序。加权trace控制的是按真实π加权的diagonal预算，不是E max_b posterior/prior；因此与Gaussian/真实通道max的非一致矩现象相容。稿子末节保留了这些限制。
