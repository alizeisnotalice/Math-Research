# upper_adapt082：两批最高绝对gap同步细分

**按完整080账本，先细分最高gap的64个parent，再细分当前physical frontier最高gap的256个parent；两批全部成功并停止。** eps=1/512、q=5/4及旧072 mean包络保持不变，未消费085理论或代码。

最终完整F上的截时积分证书为

\[
L_H=38079990466813616787770097274397379916416525/1461501637330902918203684832716283019655932542976,\qquad
U_H=13132353454328594128110747804495780001655001/182687704666362864775460604089535377456991567872<71884167/1000000000000.
\]

即约[2.6055386798e−5,7.1884166908e−5]。这里L_H/U_H描述 **H截时积分**，不是实际I的双边评估。用076完整E上端，B=U_H^(4/5)+E，严格得到

\[
I\le B_U^{5/4}\le 37726567/500000000000\quad(\text{内部上端约 }7.5453133514\times10^{-5}).
\]

不加064 tube。相较080的UH约2.16259124e−4和I上包约2.20694960e−4，本轮包络明显下降，但尚未达到I≤4e−5。**没有正式Cn更新。** g=1e−4按本轮当时冻结登记仅作conditional；所得conditional cost上端约56.41437562，严格<32门未过。后续g资格或085新包络必须另锁源、另登记，不回写082。

## 1. 冻结来源、数学对象与有限候选

只新增本独立clone的upper_adapt082。079 manifest为28e42068c0619a0769ed2130966c29d4f99860222a436d8298cd622c24cab3db；080 manifest为28eeeaa38ae47399413ab76f16630ffb43318d77dde3bb53bd2ee1dffc8ac55b。其全部文件逐hash验证未变。初始完整52360个L5 fullbox、准确weight及向外量化积分端点全部读自080，无空间抽样替代。

节点(L,i)是∏[i_j/2^L,(i_j+1)/2^L]的真实完整box，i非增，w=24/∏重复次数!。selfL/selfU存 **w乘真实积分端点**，共同分母2^160，故跨level直接整数相加，不遗漏体积或orbit因子。没有使用H_eps同eps空间齐次。

每批选择尚未展开的physical leaves，按准确self weighted gap=(selfU−selfL)降序，再(level,indices)升序冻结有限prefix。首批64正是080全层最高gap；次批256依据首批全部完成后的真实frontier。每个parent对四轴同步二分，canonical child=sort_nonincreasing(2i+bits)，去重16个bits。child唯一parent=floor(child/2)，准确核对Σchild weights=16·parent weight。

generic anchor只用于0阶DD资格：先真实box中心，再中心+eta beta^j（eta=width/4096，beta=1/8,1/16,1/32），再代码锁定的7个小负jitter，共11个有限候选。280投影互异与严格box包含准确验证；真实box没有移动。冻结mechanism误写“8个负jitter”，errata单独澄清为7，未改变原机制、代码、登记、候选或数值。代码与每job候选资格是准确登记源。

两个batch在任何H前保存具体parent/child列表、准确gap、源state/receipt hash、代码与输入hash、预定H数量。对每个新anchor，在SHA42c67156078cb0eb134910a5248821eb3c826c96981afb3d1d5b2fe3242a871b的066 cache作准确point lookup；本轮没有任何匹配，故3779个新anchor各73个端点都是新对象，合计275867个新H，未重算旧点H或旧jet。

## 2. 连续积分kernel保持不变

每node以017的全部3504锥，只生产0阶权w[k²]+=−det/(2∏node differences)，准确聚合140权并验Σw=0；H_t=D⁴/[144(2π)⁵]Σw log((144Dt)²+k²)。384-bit Arb中完整消去后取球端点。没有jet，没有逐锥绝对误差外推。

全部73正时间端点和H_1无穷尾给严格正质量，覆盖[1/512,∞)。真实box offset给16角点max uᵀQu、h∞、精确μ及η=E uᵀQu，Q=I/4−J/36。075 expLip与旧asinh比值向上端取min，给range；075均值方向score ||μ||∞/a及072系数21/8·ηR/a²给mean，再与range逐band取min。最后用真实[L,U]及mean上端的凸弦积分，所有range/mean/积分端点向外量化2^−160。旧globalHcap3仍合法。没有消费085的25/24或其它未登记系数。

独立读回全部3779个stencil用符号递归DD重建，逐权等于生产的重心公式；验证零和、二/四阶矩的符号与半因子。全部275867个band严格正。另用Chebyshev式old R、半exp平方、直接二次均匀矩和另一凸弦表达式读回；不重算H。两批各自读回通过后才派下一批。

## 3. Physical frontier与数值hierarchy严格分开

physical frontier是实际未展开terminal fullboxes。只有parent的 **全部登记children成功**，才退休parent并安装所有children；半完成children只写cache，不改有效frontier。source080一直覆盖F，完整child orbit又精确替代parent orbit，因此任何完整状态的weighted体积均为1。

数值hierarchy保留所有历史self bounds与parent/child边。独立递归规则为

\[
\operatorname{effU}(n)=\min(\operatorname{selfU}(n),\sum_c\operatorname{effU}(c)),\qquad
\operatorname{effL}(n)=\max(\operatorname{selfL}(n),\sum_c\operatorname{effL}(c)).
\]

self与完整child sum都界定同一个parent full orbit的积分，其交集有效，且必须effL≤effU。数值上端不要求由单一physical frontier实现；numeric_hierarchy保存两类不等式。即使raw childsum变差，也保留所有children及cache、继续该child，并以旧self支付上界，不丢数据或重跑。本轮观察到raw childsum更差parent数为0；这不冒称已经压力测试退化或恢复场景。

两批累计展开299个L5 parent、21个L6 parent。数值节点56139；physical leaves55819，分为L5=52061、L6=3526、L7=232。全部边的同步守恒、parent唯一性、selection prefix和antichain更新由独立迭代路径核验；weighted volume准确为1。原080 self bounds和所有新child self bounds都保留。当前physical raw积分和与有效根和相同，但未来有变差分支时两者可不同。

## 4. 两批实际结果与条件32测试

|阶段|parents|children|新H|全域L_H|全域U_H|I上包|conditional cost(g=1e−4)|
|---|---:|---:|---:|---:|---:|---:|---:|
|080输入|0|0|0|2.588366755e−5|2.162591241e−4|2.206949604e−4|165.0079700|
|082首批|64|677|49421|2.594702204e−5|9.030144069e−5|9.403398330e−5|70.30680116|
|082次批|256|3102|226446|2.605538680e−5|7.188416691e−5|7.545313351e−5|56.41437562|

采用B_U的准确有理外端，严格32门是5^5 B_U^5 < 32^4 g^5，两次均失败。另用截时L_H与 **固定E_upper** 构造该budget表达式的下诊断；当前没有排除该固定E预算路线未来<32。该诊断不是实际I或C5的下界，不能用来声称真实常数≥32。

剩余最大未展开weighted gap为726638414397412287595329828789934094169/22835963083295358096932575511191922182123945984≈3.1819915444e−8，来自L5[31,31,17,17]。最高十项如下；本轮到两批边界停止，没有再派第三批。

|node|weighted gap|
|---|---:|
|5:31,31,17,17|3.18199154442e-08|
|5:30,28,26,24|3.16914981213e-08|
|5:31,20,20,19|3.16433976875e-08|
|5:31,31,19,18|3.08069774972e-08|
|5:31,30,28,24|3.07575372777e-08|
|5:31,29,27,24|3.07203098365e-08|
|6:62,61,61,61|3.04285457412e-08|
|5:31,19,19,19|3.03928534076e-08|
|5:24,23,23,23|3.02581202385e-08|
|5:31,20,19,19|3.02080045935e-08|

## 5. 运行、读回、恢复与资源

新程序无wall alarm/timeout。首生产主PID72229，12worker72244–72255，实测wall11.482493秒；次生产主PID72620，12worker72651–72662，40.867541秒。逐512锥/8端点原子保存同point精确prefix；同源码/fingerprint恢复只复用完成的cell，不静默重放。两个程序完成后均close/join全部12worker并退出，没有无限controller。

首读回PID72430，6.804796秒；次读回PID72882，28.763485秒，全部递归DD/alternate包络/层级通过，两个12worker pools均join。可携读回另实测约0.6655秒，检查3779条量化凸弦、完整selection与同步orbit、physical frontier和数值层级以及10个端点样本SHA。全部production/verification worker PID现在均不存在，partial checkpoint为0。实际未发生续算或算术故障，不声称已测试中断恢复。

源码/输入identity：global registration SHA c03abd61687227db1925a6fc190a7c61b8adfba5d3950c683840cbf45662001e，manager SHA e9171690f008b8fd6ab4ed6ecaca01b5f4e94588a93a00dfab6279f12a66ba77，kernel SHA 428a25eff59858f0a26ced44728eb3850889cfe3a3092525194b10028e4e9be3。首批reg9be73b44c01b20b572077c32ec8e698df987559d0f26d250bdd22bc13e0856a9；次批regaf7d749fb71dd6e13aecc2be73f9ddd9b6b3d19cabb49430fdf244a3f2975497。receipt与readback hash全部归档。

本机Python/python-flint0.9.0提供数学计算；没有MATLAB、提交、暂存或推送。解析部分为已锁理论与本轮同步partition/intersection证明，数值为准确有理/Arb包围。087在独立clone审计首批，作者读回不替代其独审结论；本报告不预报087结果。

## 6. 最小可携包和公开重建

`python -B verify_portable.py`默认只读，依赖本目录及相邻冻结079/080包，不需要旧659MiB cache或本机raw endpoint目录。作者本次--record结果在portable_verification.json。batch_frontier是物理对象，numeric_hierarchy是数值对象，child cells/registration保存全部新局部积分与完整展开边，可从080准确重建，不复制原52360条root数据。

完整3779份stencil/73端点在本机 `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_adapt082/evidence/` 保留，mathematical_hashes保存移除runtime后的确定性SHA，endpoint_samples保存10个完整样本。

公开单cell重建：`python -B reconstruct.py --key LEVEL:I,I,I,I --output NEW_DIRECTORY`，只接受新输出目录并比较数学SHA。入口保留原登记point、eps和73端点，绝不对H_eps使用空间齐次。本轮没有为了测试入口重跑已成功H。

只有这两批已获执行；manager再次运行会拒绝第三批。085新mean若以后采用，须以新阶段/源锁/登记复用本轮band masses，不改082旧self证据。必要文件为报告、冻结机制/errata/登记、生产/读回/重建源码、两批输入/证书/receipt/frontier、数值层级与少量端点/数学hash/summary/manifest；大raw/stencil留在外部scratch。
