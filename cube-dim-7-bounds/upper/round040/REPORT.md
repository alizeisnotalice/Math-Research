# 第040轮：一般覆盖精化、质量接口与三条范围反例

本轮严格完成一个old-x覆盖精化与三个共同mass/选择接口；都已分别独立审核。没有七维专项或Cn下界计算。原算子仍为Lebesgue、中心轴平行cube、坐标共同半径、一般L1与全部正尺度。R042后来改变scalar尺度及截断方式，取得更强阶，按独立轮次封存；不能把本轮模型成本误记作其算子改进。

## 正向精化：全反对角带的能量

对任意固定0<γ<sqrt2，存在Nγ，使Cn≤T−γ n sqrt(loglogn/logn)，T=nlogn−nloglogn。取L=logn,l=logL、x=L−3l,t=n/x、d=c sqrt(l/L)，c<sqrt2，h=1/(L d l)，u,v∈[εd,(1−ε)d]，连接|u+v−d|<h。原两侧goodset使用tol t/l、1/l和K统计，Chebyshev使坏target比例o1、source总坏质量o(k)，k=nhL^(−c²/4)。

deleted-product统一为L^(−c²/4)(1+o1)，捕获面总数2n(1−o1)，系数为n(1+o1)，不再只用n/2。准确反对角面积2hD−h²，D=(1−2ε)d，给总energy≥2knD(1+o1)，两侧degree≤2k(1+o1)，τ/k→0只需fixed c²<2。因此同一真实history缺额≥nD(1−o1)。若γ<c(1−2ε)，全finite-a实际捕获与原一般输入桥接给Cn≤F(A,B−γn sqrt(l/L))<W−该项。T−γ推论须先取更大中间γ支付W−T正误差，原proof已明确；不是把o1当有限零。

energy_refinement/PROOF.md及audit/round040/energy_refinement为全链ACCEPT。门槛只对固定strict参数共同充分大，未认证128/256或端点sqrt2，更未称global最优。R042更强阶会渐近盖过这个结论，直接比较完整cost另见其comparison。

## CF即使使用真实f预算仍会pointwise过支付

取R037的真实全selected共star高load历史及唯一fη=(1+η)λ1_unionstars。每个Mi=(1+η)λ|Si|，整个允许预算区间|Si|≤hi≤Mi/λ内都有wi/(1+η)≤pi≤wi。原all-source eligibility更强性质保证仍全部selected。于是拒绝集空，R032 cover-flow没有任何可调flow，exact证书L*=Σvi/hi，最优hi endpoint给Wold/(1+η)；取η趋0仍leading1饱和。

同一有限cover却可由共同top witness支付λU/||fη||≤c_a近1。这只界有限U或selected union比，不反向界整个levelset弱比。cheap高load probe可远离并按homothety缩到任意小质量而sup不变，说明仅凭‘存在高load点’的全局dichotomy不够；它不排除f加权热点质量法。截断输入时保持预算/选择的断言限fixed h或局部mass-covariant规则，不适用于任意全f oracle。

mass_sensitive/PROOF.md及独审audit/round040/mass_sensitive为ACCEPT；两scope修订的旧新身份在editorial_adaptation.json，原代码/计算receipt保持原真实字节，未重跑旧toy。新增weighted cluster primal以真实Mj作成本可合法支付，但dual只给此目标下界，不是Cn上界；尚无uniform估计。

## common-star index-cover的树障碍与体积修复

任意固定n构造深D二叉层状stars，R_l=8^-l、子centre偏父±R_l e1/2，原Qi半径R_l/a。a>8/3、1≤β≤2时所有小Qi严格不交、βstars兄弟分离，真实weight a^-n；Da^-n≤τ保证全部eligible。所有common-star簇都是祖先链。全列fractional index-cover最优严格K=1+D/2：非负chain dual与完整向下路径primal都明确，不是只对称猜测。

给每leaf等质量同一L1输入并足够总质量使root真见证，深k的Mj=F2^-k，则weighted Sf/F仍=1+D/2；即仅把cost改成真实f mass并未消除index过支付。即便common-star列按真实volume费c_J=Σvi/|Srep|，仍≥a^-n(1+D/2)。这些深度结论须保留Da^-n≤τ，不能fixed a/τ无限D。

明确2N列修复是containment H_j={i:ri≤rj且Qi⊂S_j}加singletons，按实际disjoint vi成本收费。原finite-instance map严格，保留旧W可行解；在树上root列一次覆盖全部，exact最优=a^-nΣ_l=0^D(2·8^-n)^l。不是任意无限簇/无结构LP替代，但尚无全几何Cn改进。star_lp/PROOF.md及audit/round040/star_lp独立ACCEPT；作者8例224列714断言与另一新near-boundary微例均实际运行，有限check不代替all-n解析证明。

## Fresh-mass新选择与仍待估的causal spill

按半径递减，只有真实fresh质量m_j=∫_{Vj去priorV}|f|≥θλ|Vj|才选择，θ∈(1/2,1)。冻结g_j=|f|1_{fresh region}，支持不交，且Σselected |Qi|≤F/(θλa^n)。被拒原witness由更早且R_j≥R的mass共同支付。以T_D,a为只累计R_j≥R、在Q(x,(1+1/a)R)取平均的causal maximal halo，归一level(1−θ)λ/(1+1/a)^n后得到单一spill常数Ψ，完整finite compact→generalL1/allr给Cn≤Ψ。粗Ψ≤(2+1/a)^n/θ不依赖Cn而非循环；它本身不改进旧界。用T≤Mg代回会产生递归损失，不能据此假闭合。

common-spike历史新选择只留首源，但一次拒绝可需要任意多源，各只贡献1/N captured mass；全n≥2的strict新strip几何反例与有理质量账本已给证，固定单源份额归属失败。selector/PROOF.md及audit/round040/selector独立PASS；另一新微检n4,N3,θ5/8等与作者三个旧例分开，未重跑旧baseline。uniform Ψ改进仍缺。

## 比较、证据与未完成部分

正向γ精化是完整算子上界，tree/CF/局部支付反例是指定certificate family障碍，fresh Ψ只是合法接口；三者不可混写。所有输入和原尺度保持，未借用四维线外连续/Green小时间证明。内部新旧比较用已证W和完整F函数；新正向γ无显式N，不能在任何展示维数填结果。公开范围沿已核原始文献，内部未称发表、首次或全球最强。

最小代码各theory目录verify*.py及独审microcheck，sourcehash和原实际receipt在MANIFEST；原receipt与scope文字补丁版本通过adaptation区分。新root-independent audit目录保留其来源路径身份，归档不称跨机器自动运行。没有生产大网格、旧失败scan或下界重启。R042已审的二阶矩新阶与显式门槛将独立整合；R041论文并行修订，goal继续active。
