# U08补充来源接口：Lin二维3.879具体核能迁移什么

2026-10-25，主代理明确追加的定点原文核验；不新增核搜索/优化。按academic-literature-search技能读取指定原文及其直链证书，未重跑Lean。GitHub实际HEAD读取72c022ba09032b522c878cca2f9b6c68217d77ea，以下源码链接固定到该版本；读取hash在u08_lin2_sources.json。

## 实际原文公式与证书

[作者PDF《c2<4》最后一页§3（PDFp9）](https://coolrmal.github.io/assets/papers/c2-less-than-4.pdf)确有显式核，不只是摘要。α=1，r=|u|+|v|，R=3797/2000，ε=159411/200000：

K(u,v)=(R/r−1)_+/(R−1)−ε(1−r)_+²，r>0。

原文声称K≥1于单位diamond、A K≥0原点外，A=−|D_u|−|D_v|。其半质量准确为R²/(R−1)−ε/6。本机Fraction单独重算得2787954611/718800000=3.8786235545353365<3.878624；这只是核对原文常数，不是新的3D/4D试验。

原始[Lean Certificate.lean](https://github.com/CoolRmal/centered-maximal-constant/blob/72c022ba09032b522c878cca2f9b6c68217d77ea/CenteredMaximal/Cauchy/Certificate.lean)列明γ=εR(R−1)=1087694569899/800000000000>4/3，r0=3457/4000。作者将所需符号化为γT(r)≤J(2r/R)，小r区间直接定号，大r区间用强凸性及r0处有理log包络；不是密集采样。J定义/尾密度的角向估计在[BaseDensity.lean](https://github.com/CoolRmal/centered-maximal-constant/blob/72c022ba09032b522c878cca2f9b6c68217d77ea/CenteredMaximal/Cauchy/BaseDensity.lean)，bump切片和角向Jensen在[BumpGenerator.lean](https://github.com/CoolRmal/centered-maximal-constant/blob/72c022ba09032b522c878cca2f9b6c68217d77ea/CenteredMaximal/Cauchy/BumpGenerator.lean)。具体代码只按所用段落读取，未在本机重建整个Lean项目。

[作者网页§4.3](https://coolrmal.github.io/articles/two-partial-balayage-principles/)直链[Palomar条目version3](https://palomar-registry.org/entry?id=PALOMAR-2026-09-19-000002&version=3)和[源仓库](https://github.com/CoolRmal/centered-maximal-constant)。PDF另列α6/5的spline候选并明确未获Lean核；本补充不把那个报告值当已审证书。

## 维数可迁移接口与不可直接迁移部分

1. **可迁移的是一般比较原则。** 对每个n，若真实坐标stable的障碍分解、K偶非负L1、AK≥0原点外及目标全域支配均成立，就能用∫K作常数。U08补足n3/4 α1的分析接口；不补出新的K。

2. **任意固定可逆线性换元保弱型常数，但必须同时搬运形状、归一和生成元。** 二维T(u,v)=(u+v,u−v)，|detT|=2，T(Diamond)=cube。原公式K的cube归一majorant为H(x)=K(T^-1x)/4，故∫H=(1/2)∫K。搬到square后的生成元符号为−(|ξ1+ξ2|+|ξ1−ξ2|)，不仍是原cube坐标轴Cauchy。一般n的线性推送也合法，但没有自动把n维crosspolytope变成cube的矩阵。

3. **diamond→square是二维特殊等价。** n维l1单位球有2n个顶点，cube有2^n个顶点；可逆线性映射保顶点数。n3为6与8、n4为8与16，故不存在上述形状等价。即使K≥1_l1ball，也不能把该支配直接宣布为cube支配，必须重新覆盖/归一并实际付成本。

4. **核心势形状和角向一维化不能照搬。** 二维真实Cauchy势恰为1/[2pi(|u|+|v|)]，这允许同一diamond半径下精确拆势、尾和bump。一般n的真实势应从G_n(x)=pi^-n∫_0∞t^n/Π_j(t²+x_j²)dt取出，齐次次数1−n，通常不是l1半径的单变量函数。特别n3/4在非零坐标轴上的该积分也发散（小t阶t^(2−n)），不能将二维1/r或直接r^(1−n)模板冒充真实势。该事实不妨碍G作为局部可积的分布势，但使源/截断的分析必须重做。原J/T的二维两切片角向不等式没有在本补充中证明n3/4版本。

5. **不能用两个二维常数直接相乘给四维弱(1,1)。** 尽管D2×D2能线性映成4D cube，若将两个有限正源核直接相乘，L=L_x+L_y下会产生−W[δ0(x)K(y)+K(x)δ0(y)]，负源支撑在两张非原点平面上，而不是仅总原点。准确地写L_2K=μ−Wδ0且μ({0})=0，两个正交项对这些平面的质量为0（K dx绝对连续），故负面源不能被忽略。这不是所有非乘积耦合核的障碍，也不是乘积算子弱型常数乘法定理。

## 具体后续建议（未运行）

真正可沿原模型推进的是：保留“真实势的原点负源+可准确控制的正尾−局部修正”结构，针对原cube或明确的线性预像形状，重建n3/4全部角向/轴源不等式与质量。选定修正时先给完整AK分布符号接口，再有限寻找成本。二维R/ε/J/T数值不能未经新推导直接用作n3/4参数。该建议未启动新试验，本批无新U3/U4。

## 主线程指定的下一轮方向（本轮不启动）

沿二维cube的真实推送生成元，改用cube体对角Cauchy frame：n3的v=(1,±1,±1)四方向，n4的v=(1,±1,±1,±1)八方向，L_V=−Σ_v|v·D|。这不是乘两个二维常数，也不是本轮新核。主另轮预登记并负责具体谱锥/真实势构造；本分支不重复。

U08分析链的可迁移理由是结构性的：有限spanning frame中任选n条线性独立方向，其独立一维Cauchy和经可逆矩阵有密度上界Ct^-n；其余方向的独立和对应保概率卷积，不增L∞。故同级transient势估计可用。有限spanning又给Σ|v·ξ|与|ξ|等价，生成元域H1、频率尾及有限方向截断仍可按U08付账；大尺度cutoff仍O(R^-1)。核比较的CND/无限Lévy源链仅用上述符号控制和每方向二差分，不要求原坐标乘积。这是另轮扩展接口说明，未独立审新frame全文，更未给其真实势公式、cube支配或低成本证书。
