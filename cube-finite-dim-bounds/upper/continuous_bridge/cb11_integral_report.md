# 11H：固定tau1/3积分认证的有限区间诊断（未达目标）

2026-10-25。原目标I=∫_[0,1]^2G(1,y,z)^(3/2)≤57/100000。本分支只cb11_integral_*，不攻g、不优化tau、不作点采样成本搜索。初始及唯一一次有原因表示修订均完成，各3952项，未追加网格/精度。**未认证I目标**；以下巨大的数是粗区间上界表达式，绝不是实际I或I下界。

## 1. 已执行范围与失败位置

空间从4×4方格出发，对(1,0),(0,1),(1,1)三个奇性角各沿触角子格二分26次，删三个边长2^−28的末级square。剩247个不重叠boxes，连同删格面积准确和1。每删格在半径10^−8帽中，因为sqrt2·2^−28<10^−8。

每个box付15个有限time区间和1个无限time尾项，247×16=3952≤4096。128bit Arb/Acb对整个空间/time矩形评价全部48反锥对，不删除负实部/zero knot。每个因子t−ib有严格正real part，先逐因子reciprocal后乘积，避开无害人工分母。

|机制|epsilon|time T|所得I上界表达式展示|结果|
|---|---|---|---|---|
|初始逐box短time|1/64|16|约83930.805|未达|
|唯一修订：全局Lq短time/Minkowski|2^−26|16|约1.0110877×10^36|未达|

初始H_epsilon幂积分粗上界约6.71；修订后约1.01×10^36。缩小epsilon后区间取消更差，不是实际H增到此量级。数据保留两次各247个remaining boxes、每box位置/面积/r_Q下界/H区间/误差贡献，见cb11_integral_results.json。

最坏box不是帽附近，而是[0,1/4]²，距最近face奇点平方≥9/16。其H区间上端初始约10.2，修订约5.07×10^24；真实密度的若干锥项在zero投影层相消，空间区间在相消前包围导致严重膨胀。帽细化或提高浮点精度不能被本结果当成已解决表示问题；需要整体confluent/合并评价。它也不是生成元或候选的障碍。

## 2. 帽与single quadrant的准确换算

030已审n3 tau1/3每face8个交点，六face共48帽。由坐标符号/置换对称，完整cube surface积分为24I。把三个删格反射到全部face/quadrants，内部不交（边界零测），且均包含于030的半径rho=10^−8帽。故遗漏I上界可取

    C_cap=(030 all-face48预算的准确上端)/24，

实际约6.86112589459×10^−8。未取E的下端作实际积分下界，未把48帽直接当一个quadrant或漏掉24。对H≤G，同一预算也包其删格部分。Voronoi最近交点的前提按030保留。

## 3. time大尾、短time与空间区间均已支付

固定tau1/3的精确peak是p_t(0)=15/(544pi³t³)。因此无限time尾

    ∫_16^∞p_t(x)dt≤15/(1088pi³·16²)，

约1.73689×10^−6，作为每box第16个项加入[0,tail]区间，未把time截断当0。

H_epsilon的前15段用Acb矩形range乘time宽度。空间y/z在整个box中同时包围，没有仅用中心点；p_t真值由完整cb11复锥和给出，配对只利用严格共轭。此representation的误差虽巨大，但没有未支付的time quadrature或空间cell误差。

030的face小time界给，r=最近交点distance、A=55200/pi³：

    0≤R_epsilon=G−H_epsilon≤A log(1+400epsilon²/r²)。

初始对每box用r下界，将此差加到H上端，再用area·upper(G)^(3/2)求和并加帽。epsilon约0.01即便r约1也付几十的Green误差，远大于目标。本界保守不等于真实R很大。

## 4. 一次修订：已闭合的全局Lq误差预算

改为先积分R的幂，再Minkowski。q=3/2，face按最近交点作Voronoi分片；全face8个中心各扩大为完整平面，再除4到positive quadrant。因此

    ∫_[0,1]^2 R_epsilon^q
      ≤4pi∫_0^∞ r[A log(1+(20epsilon/r)²)]^q dr
      =K epsilon²，
    K=1600pi A^(3/2) J，
    J=Γ(5/2)ζ(3/2)/2=3sqrt(pi)ζ(3/2)/8。

r>1没有被假定为真实face区域；它仅是把正径向包络扩到全平面的积分上界。原face最近距离≤1，030原假设足够。8/4、2pi、20²全数支付，与全48/24同义。

J身份：令v=log(1+u^−2)，积分分部得

    ∫_0^∞u log(1+u^−2)^qdu
     =q/2∫_0^∞v^(q−1)/(e^v−1)dv=Γ(q+1)ζ(q)/2。

端点因q>1消失，正指数级数允许Tonelli。epsilon=2^−26时实际Arb给K epsilon²约1.4557474945×10^−7。若B是整个quadrant的I_H=∫H_epsilon^q上界，则

    I≤[B^(2/3)+(K epsilon²)^(2/3)]^(3/2)。

实际准确条件检查：**若可证B≤9/16000**，则右侧严格小于0.000565930019122<57/100000。这个if前提目前未证；不能将数值H/G积分代入当I证明。

唯一修订使用此全局预算、同一247 boxes及15+1个time项，epsilon改2^−26。但H空间/time区间膨胀，B未获得足够小上界，故I仍未达；新增的内容是完整可复用的平均短time预算，而非一个成功I证书。

## 5. 主新face-marginal预算的分析核对（不新增实验）

本批计算后，主给proofs/round011_face_time_budget.md。本代理已逐步读核，不重复其det枚举/数值批：

- 在x1=1，active k个参数至少一个|Zi|≥1/k；restricted density与n−1条body补基、余概率卷积给p_t≤C_f t^(2−n)。用邻近plane的1−eta阈值再连续延拓可避免零概率slice误用。
- 真X1 marginal是scale At的Cauchy，L1(plane)≤At/pi；sup^(q−1)乘同一密度的marginal得到time指数1/n。
- Minkowski积分给n/(n+1)·epsilon^(1+1/n)；positive orthant的norm因子是2^(-(n−1)^2/n)，不是2^(-(n−1))。

未发现分析缺口。n3由主锁定k8,A16/3,C_f=(256/3)/pi³；epsilon1/2000可把允许I_H阈值提到约.000567072834，比本批距离法所需epsilon大得多。它是更有希望的后续预算，但尚未在此批另作time/grid实验，也没有把旧不同cutoff的H上界retroactive套入新epsilon。主新预算自身独审状态另记，不能因本读核称第三方独审通过。

## 6. 算术、来源、范围与下一步

cb11_integral_diagnostic.py实际只运行初始与一次修订，约0.2秒/次，Arb/Acb128bit、Fraction endpoints。保存的checked rational I upper包住对应上界表达；这些下端不是实际I下界。原cb11模块/旧证据未改，源hash与030准确端点在cb11_integral_sources.json。

最小budget checker只核覆盖、帽换算及上述条件式，不冒充3952项的独立重建或I目标通过。MATLAB只交解析预算cross-check给主串行，不重复空间/time诊断。共享Python/SymPy/Fraction/Arb及cb11锥表达明确披露，完整表达的独审状态不因本计算自动改变。

本批已按登记停止：固定tau，初始≤4096、一修订、所有247 remaining boxes和各预算已保存；没有点样本/Gauss趋势替代。主后续Cheb/Bernstein新机制将另登记，本分支未私自启动。当前尚缺连续H_epsilon幂积分上界，g由minimum独占处理，最终成本/全L1迁移由主及原独审承担。U3=8、U4=16本批均未改变。

主串行MATLAB实际exit0（session10679），执行script SHA b98241625e9598c5f42cf41a26b33750330e8e2c65d9fa3bf21c8dcdd4fdd802，J=1.7363680300455859、条件I=.00056593001912155804；仅解析budget校准PASS，实际I仍未认证。详见cb11_integral_matlab_receipt.json。
