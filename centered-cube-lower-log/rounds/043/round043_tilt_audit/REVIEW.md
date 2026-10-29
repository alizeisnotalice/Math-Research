# 第43轮脚本与结论独审

`checks.py`已实际运行，结果 `CHECKS.json` 为PASS。未修改主来源代码、PLAN、SUMMARY或其他轮次文件。

## 数学判断

主目标τ_ν≤1/[q(q+1)]成立，完整自含证明在 `PROOF.md`。正确归一化histogram prior是Mult(n,μ_ν)，q1的trace≤2。不能沿用原η：一般Jensen下界为(Σμ_c²/η_c)^n，只要μ≠η就指数增长。这个下界已用主代理新给出的真实正来源、t1 score倾斜有理例严格实现，审计另用全长度5/4物理中心窗独立积分得到相同分数。没有发现数学方向错误。

MGF倾斜恒等式精确成立，但必须保留m(t)^n；其后的trace付款只覆盖用新π_t归一化的type平方和。scalar event posterior和type向量平方和是不同对象，证明已分别列出。不同table的Q_t、m_t、μ_t仍不同，尚无共同query律或联合付款。

## 主9000数值脚本

读取并核对 `round043_tilt/probe.py`，其source helper为 `round042_channel/probe.py`；两个SHA与本轮PLAN一致。主posterior中r是各class的真实周期残余窗积分，a=(q+Σr)/L，p=(q/K+r)/(L a)，对应η=1/K。继承source每格质量1/K。full长度q+ell与canonical residual ell窗只相差固定query平移，积分/trace不变。

score倾斜使用Gauss权w×a^(1+t)，因为基准Q=a dx，指数1+t正确。随机极化和两点例直接定义节点上的离散query律，不必额外乘原Gauss权；它们是不同任意ν例。trace每次重新计算μ=w@p，并用ΣE p²/μ−1，使用了正确新prior。9000=360channel×25law，其中3240score倾斜、2880随机极化、2880两点，数量核对通过。

本审计只读调用其posterior输入生成，按完全相同RNG顺序重建9000个law，用独立centered variance `ΣE(p−μ)²/μ` 计算trace。各组极值与原SUMMARY相差<2e−12；逐law两种trace最大差2.10e−15。两点组原最小−2.22e−16是减1抵消的浮点残差，不能当成负方差。

Gauss–Legendre正节点权、source分段原函数、16阶积分及8/16/32/64加密依赖42helper。主SCRIPT的BASE备用路径适合归档到rounds/043后引用兄弟042；当前运行实际使用同一output根的42helper。复现必须同时保留该helper并核对SHA，不能宣称43probe完全无跨轮文件依赖。只读重放不调用main，因此未覆盖或重写原SUMMARY。

加密结果支持数值稳定性，不是向外区间或积分误差证书。数值扫描只覆盖列明K/q/ell/source/相位和25种节点律，任意ν、任意K的结论依据一般证明，不依据有限扫描。所有探索都是确定性积分/随机生成候选ν，不是MC估计C_n或盲测。

## 额外归一化实验

审读 `normalization_probe.py` 并重算7个Fraction二项展开及8个实际score-tilt μ/factor。原分母toy量恰(10/9)^n，新分母恰1。实际例的μ偏移与百万维log下界仍是浮点描述，不升级为strict physical证据。

审读 `exact_tilt_example.py`：source质量、残余窗积分、仿射乘积的精确积分式和独立80cell Simpson的网格对齐均正确。审计另实现真实全长度5/4的周期box交长，使用完全不同的窗端点±5/8分段，9段Fraction积分精确得到原JSON的m、num、μ和Jensen因子。这关闭了仅用Dirac/抽象posterior例的物理实现疑问。

## 交付

`PROOF.md` 给出一般严格证明与第一处未支付费用；`checks.py/CHECKS.json` 给出实际运行证据。无C_n增长阶提升，未提交Git，未改写来源文件。

## leading price补充独审

后续 `price_probe.py/PRICE.json` 已审读。主提出的谱定位成立：trH_η≤1+κ、γ≥||H_ηu||²=1+χ²(μ||η)、其余谱和≤κ；原先验hist trace在γ^n与γ^n/(1−κ)之间。u方向块分解的Schur公式进一步给γ≤1+χ²/(1−κ)。自含推导已加入 `PROOF.md` 第8节，没有把新prior的Markov收缩错误地用于H_η。

独立物理全长5/4积分另计算z_cz_d，得到与PRICE完全相同的精确Gram [[822002,721090],[721090,770000]]/1517091，并用其特征多项式精确认证γ∈[1+χ²,1+2χ²]。144个有限law的trace帽为Fraction检查，谱夹逼为实际NumPy检查，资格在JSON中分开记录。主百万维log price只是浮点辅助量；严格一般结论是γ^n的不可省代价，仍没有C_n新阶。
