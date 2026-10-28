# 第34轮独审：分块共享相位与多项式v预算

## 当前状态

已读取round034_block_plan.md及原来源/概率/主证明，并实际运行本目录两个小脚本。分块共同尺度、条件矩、全混合提议密度和贡献帽通过审查，后续已只读核对round034_block_mixture实际pilot脚本与输出（见第6节）。完整weak-profile理论及实际实现/保存统计分别见第7/8节，未确认G8的最优weak值超过G1。随机jitter核最终稿通过独审，详见本目录RANDOM_JITTER_AUDIT.md。本人未重复代理模拟、修改旧文件/Git或启动其他代理。

整数预算给出可保留现有原来源充分资格的非空多项式窗口0<beta<1/9，其中v约n^beta，K约n^((1−beta)/2)。它只扩大了可供混合探索的尺度预算，原单分量一阶和仍精确抵消为K sqrt(v)=c_p sqrt(n)。本轮小例没有提供新的渐近下界。

## 1. 已实际运行的分块小例

block_exact_probe.py采用4坐标、4个等概率状态、每bin含2个隐藏状态、2相位和2共同查询；所有Y正且均值精确1。它是验证一般代数的有理小例，查询Y没有被声称来自同一个cube平均来源。

用Fraction穷举16组bin条件及全部256个真实状态，得到96个条件矩检查：块条件均值与直接枚举严格相等；较松的有理Cauchy方差帽严格通过，计划中的较紧Minkowski方差帽另以浮点逐组核对通过。后者小例核对不替代下面的一般证明。全G=1,2,4与两查询的带潜变量提议RN相等、总质量精确为1，W≤2。错误地先按坐标平均相位会得到不同来源，所测条件均值最大差为97/256。

真实共同查询弱系数也非随G单调：lambda=2时G1为31/128，G2/G4为1/4；lambda=1时分别为33/64、127/256、47/128，方向相反。这只排除「块数越多越好」的无条件直觉，不支持原来源的增益定理。

实际计算了块各自选查询的乘积U_G与真实共同查询V_G之比。G2/G4的U_G=V_G概率分别为17/32、1/8；平均log(U_G/V_G)约.138284/.428520，最大比为37/20、9/4。可将此作为真实pilot的描述诊断，以区分块内相位平均的收益与共同L损失；不允许以U_G尾代替V_G尾。

## 2. 一般条件矩与提议顺序

固定G个坐标块，F_G=prod_g avg_r prod_(i in g)rho(x_i−theta_r)。每块密度积分1，故F_G积分1。cube在各坐标的积分可用Tonelli和乘积结构展开，精确得到同一L下的prod_g avg_r prod_(i in g)Y_L(x_i−theta_r)。所有块必须使用共同L；也不能将avg_r与坐标乘积交换。

保存8个最细块的counts，并将相邻细块合并为G=1,2,4,8的coarse块，足以重建各G的条件矩。不能先合成一张总直方图后恢复分块来源。给定这些fine counts，原律下不同细块的隐藏bin位置独立；coarse块由不相交细块组成，故仍条件独立。块内不同相位共享这些位置，一般相关。

每coarse块/查询/相位的条件均值A_gjr为对应m_qjr的计数乘积。块均值M_gj=avg_r A_gjr；令a_gjr=A_gjr/(R M_gj)，各相位条件相对方差帽为r_gjr，则L2三角给块相对方差帽r_gj=(sum_r a_gjr sqrt(r_gjr))²，不使用相位独立。不同coarse块的实际条件二阶矩乘积，给总相对方差帽prod_g(1+r_gj)−1。先形成此共同L总均值/方差，再对每个预定lambda求Cantelli，并在共同L中选择最大下界，合法。按均值选L也合法但可能较松；分别给块选择L不合法。

提议的潜变量顺序应为：选uniform分支，或选G和一个共同anchor j，再对每个coarse块选相位，此块内所有细块使用同一相位。相对原细块counts律的RN为

    D = 1/2 + 1/(2 H J) sum_(G,j) prod_(coarse g)
            [ R^-1 sum_r prod_(fine h in g,q)(m_qjr/Z_jr)^N_hq ],

其中H是候选G数。multinomial组合系数在RN相除时准确抵消。这里是整计数样本的混合，不能把单坐标混合再取n次，也不能只用实际抽中的一个latent标签作为边际D。D≥1/2，故W≤2。

若所有目标L在proposal中，C_Z=exp(max(0,n max_jr log Z_jr))，则proposal中匹配G/j项至少M_j/(2 H J C_Z)。在该Cantelli下界为正时M_j>lambda，故lambda W B≤2 H J C_Z；另有2lambda帽，取两者最小。每块相位可能不同仍只付Zmax^n，因为块长和为n。此帽不含R^G，但依赖proposal确实积分了所有R^G种块相位分配。baseline上界项一般只有2lambda帽，不能误套active下方向帽。其条件矩可用总counts，因为它是一个真正逐坐标乘积源。

## 3. 实际整数/Fraction参数筛查

polynomial_budget.py选p=1/2、n=2^m，m=720/1440/2880；beta取1/16、1/12、1/10、1/9、1/8、1/6、1/5、1/4、1/3，保证beta*m为整数。选择使v=9n/(116K²)≥n^beta的最大K=2^d−1，B为不小于实际v的最小二幂，M=ceil(sqrt(v))。K、v、B、M、H(B)、capacity及误差上界均用整数/Fraction；输出的log2仅作显示。

Cx使用保守有理上界：a_B≤1、sum m^-2<2、H_M≤1+ceil(log2 M)、pi<22/7、sqrt(K+1)向上取整。有限面显示上界用sqrt(v)≤M、sqrt(log n)≤ceil(sqrt(m))，没有浮点决定gate。没有生成这些巨大维数的显式本原多项式或标签；q>H只应用已有存在性引理。

| beta | m=720时d | log2(q/H(B)) | q>H | 必要比率无邻接capacity |
|---|---:|---:|---|---|
| 1/16 | 335 | 145.8301 | 是 | 是 |
| 1/12 | 328 | 82.8301 | 是 | 是 |
| 1/10 | 322 | 28.8301 | 是 | 是 |
| 1/9 | 318 | −7.1699 | 否 | 是 |
| 1/8 | 313 | −52.1699 | 否 | 是 |
| 1/6 | 298 | −187.1699 | 否 | 是 |
| 1/5 | 286 | −295.1699 | 否 | 否 |

三组m上，前三个beta均通过q>H，而beta=1/9均失败。beta≥1/5的所测常数下连R_B=(B−1)(B−2)个不同非相邻循环比率所需K≥2R_B+1也失败。不能据q>H失败宣告原类型标签不存在；1/9≤beta<1/5仍可能需要更有效的显式构造/存在性方法。beta=1/5的边界还受v的前导常数影响，不能仅凭此表作一般不可能性断言；beta>1/5时该capacity的幂次则确实冲突。

## 4. 原证明预算逐项变化

记v约n^beta、K约n^((1−beta)/2)，保留原B约v、M约sqrt(v)。下表是沿用已有估计的充分消失条件，不将估计失效当成命题不成立。

| 项目 | 量级或要求 | 原预算的充分范围 |
|---|---|---|
| finite-field存在性 | K/v^4→infty | beta<1/9 |
| 真cross乘v | v/sqrt(K)+v log M/K | beta<1/5为主限制 |
| safe删除占比 | v²/K | beta<1/5 |
| tilt三阶误差a_n v | v/K | beta<1/3 |
| 均值平方费用 | v²/n | beta<1/2 |
| 有限面旧sigma=n^-1/8 | v n^-3/8+sqrt(v)n^-1/8 sqrt(log n) | beta<1/4为主限制 |
| 谱尾乘v | v T_M=O(v^-1/2) | 任意beta>0 |
| quantization乘v | v/B²+v/K | beta<1/3 |
| M<K | sqrt(v)<K | beta<1/2 |

fixed-sigma的1/4不是本质概率障碍；改变平滑参数可重新优化。但在beta<1/9窗口它已自动满足，无需改概率引理。原文曾以「K/v^r对全部固定r」一笔带过，新的多项式处方不再满足该笼统陈述，必须逐处替换成上表实际使用的条件；已核对source/probability/main三份detailed稿中的使用处，未发现需要比K/v^4更苛刻的已用幂次。

平滑峰留在原格内部依赖B≤wv，B<2v与w≥4足够。正背景依赖H_M/K→0，在此窗口成立。原主阈值必须改成h=v/2−2约n^beta/2，不能继续用固定h30。q_L=v+O_w(1)的核预算仍成立；near pair经r=Kv|L−J|/2重标后可求和的Ce^-cr/8仍适用。middle/far pair的总费用Kv exp(−v/25)在任意beta>0更快消失，因为其log为((1+beta)/2)log n−Theta(n^beta)。固定有限局部查询、随后T/A/D等顺序极限仍须保持，不能把固定维CLT直接用到随n增大的查询数量。

例如m720、beta1/10时，所算a_n v、vCx上界、safe v²/K和旧有限面上界的log2分别约−248.33、−85.54、−177.38、−49.09，真实有限数值上都很小；这仍非全套有限n下界证书。

## 5. 第一单分量抵消与本轮可检验目标

所有27组参数都以Fraction核对K²v/n=9/116。故单分量的原局部簇数约Kv、每簇Mills费用1/sqrt(v)，一阶和仍为K sqrt(v)=3/(2sqrt(29))sqrt(n)。局部Brownian强度58/(9L²)也不随v增大；增加v没有自动改变这一处抵消。原下界证明在通过资格后仍只给sqrt(n)量级。

新strong预算Kv约n^((1+beta)/2)确为超过1/2的多项式幂，beta<1/9时可留下小于5/9的上方空间；这只说明已有受限来源上界没有排除某些混合收益，不是下界，也不是任意输入C_n上界。

本轮应优先用同一来源、同一预定lambda比较各G的真实共同尺度Cantelli下积分与baseline真实尾上积分；若只能改善条件均值、块各自选L的U或扫描最优代理分数，不登记弱收益。可以同时记录块内相位份额与块间共同L损失，但这些只是解释信号的描述量。N128时尾贡献ESS、正负抵消和逐阈值family预算必须明示；任何择优目标的独立确认需另行冻结。固定CDF差仍不比较最优weak值，故更直接的筛选是第7节的完整profile围栏。

## 6. 实际分块pilot只读审查

已只读round034_block_mixture/block_probe.py、registration.json、RESULTS.json及REPORT.md，未重跑其MC。用本机Python标准库核对脚本SHA与注册值相等，并从保存lambda=e²的4组paired_D重算均值与表值相符。源GF(2^8)/0x11d、B4基(1,6)资格、解析bin均值、transition方差帽和精确有理query并集沿用既审来源接口；K255/v8、v16分别为n6704800/J233/Q48960/M3与n13409600/J460/Q97920/M4，均含区间两端。未满足原v≥16log n、B≥v及最大周期规则，代码与报告均未称原匹配来源。

实现保存(sample,8,Q)细块counts；每细块行和断言为n/8。coarse phase logs和方差指数仅合并相邻细块，随后Minkowski形成块帽、再将(1+r_block)相乘。相位0 baseline使用所有细块对应phase0之和，没有套用分块相位平均。各G、各阈值先形成总模型的共同L条件矩，再max Cantelli；between_G和对product的差均是lower-minus-upper。

sampling先选uniform或G/j，再每coarse块独立选phase、向其全部fine块复制该相位；RN准确使用全4J候选、每块相位边际以及Z_jr的coarse块长度次幂。没有把实际抽中的latent标签当作边际分母。归一化C_Z和2·4J·C_Z下方向帽均保留；upper项一律只用2lambda。逐行独立抽样，前后64分割的discovery/holdout互不复用随机样本；预定第二点成本gate不按统计阳性决定。

small_control的直接rho卷积采样、4^G latent展开及bin内控制顺序正确。其控制数值是探路证据，不是一般方差证明或浮点区间证书。脚本释放counts后再形成coarse工作表，内存报告明确为主要数组估算而非RSS实测。没有发现来源平均顺序、block条件独立或实际提议RN的实质错误。

每点N128、alpha=.05/1000。固定7阈值每个保留4个非负下CI、4个对product signed CI及12个有序G差CI，共140；另有5个holdout下CI和4个holdout差CI，共149/点，两点298项，1000预算足够。上积分仅报点估计、未给上CI，因此现有表不能围住weak最优值。有符号贡献取[−2lambda,2lambda]，其Hoeffding费用4lambda sqrt(log(1/alpha)/(2N))方向正确，不将signed消去比解释为尾ESS。

v16在lambda=e²的G1/2/4/8对product积分差分别为1.50532、2.24245、2.78109、3.51526，对应非负尾贡献ESS约71.85、82.16、88.00、94.11；这些固定阈值比优化扫描稳定。但全部保存的保守signed下端未正。heldout v16 G1仅约1.7 ESS，G8约19.1；phase0 product约1，不能把这些扫描点估计排名当真实weak排名。v8也没有随G统一改善。报告保留了这些失败与低ESS，没有以描述SE代替确认。

## 7. 完整weak-profile围栏的数学预审

此节先审核主线提出的新预定筛选公式，实际实现核对见下一节。令Z_Gj=A_(L_j)F_G≥0，T_G=max_j Z_Gj。展开块相位为R^G种assignment的概率平均，square Jensen给Z_Gj²≤avg_assignment product_i Y_j(x_i−theta_i)²。原均匀查询律下坐标独立，且相位平移不改变一维二阶矩，故

    E Z_Gj² ≤ (E Y_j²)^n,
    E T_G² ≤ sum_j(E Y_j²)^n ≤ H2.

这对G1/G8、甚至所有固定块数都成立，不要求块phase相关性消失。H2可由bin二阶帽作确定上界：E Y_j²≤sum_q p_q m_qj0²(1+t_qj0)，再取n次方、对j求和。若用ordinary float计算该上界，仍需明确数值未认证；不能把接近下端的浮点值当确定H2。理论上H2≥1。

弱函数w_G(lambda)=lambda P(T_G>lambda)在0<lambda≤1时≤1；在lambda≥H2时由二阶Markov亦≤H2/lambda≤1。因此仅需围住[1,H2]。固定h步长.25，lambda_k=e^(.25k)，直到lambda_last≥H2。CDF单调给每段[lambda_k,lambda_(k+1)]

    w_G(lambda) ≤ lambda_(k+1) P(T_G>lambda_k).

对所有k同时有效的尾概率上CI形成整段上围栏，连同外部值1得到完整有限J模型weak上界；任意一个固定网格阈值的有效下CI给weak下界。严格尾符号在端点保持该不等式。上/下CI及两模型若联合比较，应分配整个家族alpha；即使上下不分离也必须报告，不能把固定CDF阳性改写成weak比较。

有效上贡献帽不需要lower active条件。若模型的条件union/Markov上界U_G≤sum_j M_Gj/lambda，且该模型实际提议D_G≥(sum_j M_Gj)/(2J C_Z)，则lambda W U_G≤2J C_Z，同时≤2lambda。可以对该非负有界贡献作KL/Chernoff上CI，再乘lambda_(k+1)/lambda_k覆盖阈值段；使用下CI则方向错误。

归一化有两种易混接口：均匀phase抽取各自normalized tilt，得到prod_blocks avg_r A_gr/Z_r^m，并不等于raw M_Gj整体归一化；若真抽normalized raw M_Gj，phase需按Z_r^m加权。前者仍由C_Z支配得到上述上帽，故数学可行，但sampling、RN和文字必须一致。另2J C_Z帽仅适用于被提议支配的同一个G模型；单相位baseline在G8提议下不能无条件套此帽，其latent项可能只占R^-G，需保留2lambda或另证支配。新weak筛选比较G1/G8各自专用提议时，应为各自模型上界使用各自的M_Gj。

该围栏比较有限确定J尺度的真正weak-profile，不是所有连续L的weak范数；单纯有限网格上界不足以上界连续尺度baseline，因而即使分离，也不能直接声称连续cube weak范数或C_n最优常数改善。它比固定lambda尾差更接近本轮目标，剩余差距仍须准确标明。

## 8. 完整weak-profile实际实现及保存统计

只读profile_probe.py、profile_registration.json、profile_results.json，无新增MC。本目录profile_saved_checks.py仅读保存均值/帽/矩摘要，使用Python标准库独立反演全部CI并重算H2与整曲线包络，结果在profile_saved_checks.json。脚本、注册及全部原pilot文件SHA256均匹配；所有CI的独立复算绝对差至多5.7e−25，保存log H2与独立logsum复算差为0。这验证了保存摘要的内部算术一致性，未恢复随机counts或重新验证三角函数/条件矩。

两个专用模型均为K255/v16/n13409600/J460/Q97920/R4，G1/G8各N512，分别使用冻结独立seed202610290343/344、batch64。源query并集来自既审精确有理配置，包含两端。每个模型直接保存其G块counts；没有为G8先合并成单块。给定counts的phase均值与方差按正确块顺序计算；每个阈值重新在共同L取Cantelli最大下界，同一模型的条件union upper含Markov截断。

每行先半均匀，否则选共同anchor L及每块独立phase，再按各phase自己的normalized bin tilt抽multinomial。RN的每块phase边际、长度n/G次方归一化及全J anchor边际与实际sampling一致。上下贡献均被同模型Markov支配，2J C_Z帽无需lower active条件，G1/G8分开抽样所以此处没有pilot的候选G数因子。

文字精度小项：registration和结果将proposal简称为“normalized model block-mixture tilt”；准确说法应是「每块均匀混合各自已归一化phase tilt，再对共同anchor混合」。当各raw Z_jr不完全等于1时，它不等于raw M_Gj整体除以其期望。实际代码是前者，RN也按前者，C_Z支配仍足以证明上下帽，因此这是描述精度问题，不是实现方向错误。

预计算H2=sum_j[mean_q m_qj0²(1+t_qj0)]^n=2348496379.4705815，log H2=21.577041121795297。h网格有88点，从0至21.75、步.25，末lambda约2.791932932e9已超过H2；因此结合外部帽1覆盖全部正阈值。上CI在lambda_k处先与lambda_k、H2/lambda_k取小，再乘相邻lambda比，并与整段moment帽H2/lambda_k取小，方向均正确。所得max只用于有限J模型，未将它当作连续尺度上界。

每个模型88点的lower和upper共176项，两模型352项；alpha_each=.05/352。CI是非负有界变量的单侧Bernoulli-KL反演：lower返回下根的保守左端，upper返回上根的保守右端；高阈值零样本仍保留正的upper CI。家族同时覆盖后才取最大lower及最坏upper区间，允许该数据选择而不再付一次费用。ordinary float的alpha求和显示.05000000000000001，只属舍入；严格实现可按有理1/(20·352)及有向舍入处理。本文均按名义联合95%统计口径报告，未称区间认证。

| 模型 | 完整有限J weak下CI | 完整有限J weak上CI | 下CI赢家h | 该点尾贡献ESS |
|---|---:|---:|---:|---:|
| G1 | 6.220996019515227 | 77.07827059550362 | 4.25 | 123.1250 |
| G8 | 9.187364329875424 | 58.94412068449272 | 3.5 | 260.5064 |

两带明显重叠，未证明任何方向的最优weak分离。G1的最大lower积分点估计为16.092（h8.5），G8为14.569（h4.5）；它们也不是weak真实排名。质量权重均值分别1.019996±.032214、1.024709±.018713，与1相容。各模型batch顺序不破坏独立rows，同一模型各阈值相关不妨碍union；两个seed分离也没有被错误用于免除家族数量。

残留FP缺口包括：三角函数和解析bin均值、相对方差帽向上有效性、H2逐坐标二阶数值向上有效性及n次方累积、实际multinomial概率和RN对应、C_Z/贡献帽/CI求值舍入。保存摘要和hash通过不支付这些误差，普通MC置信也不支付；目前只得到以积分量/提议/帽正确为前提的完整有限模型统计带。没有全连续尺度weak、一般C_n或新渐近幂次/对数下界结论。
