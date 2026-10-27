# 第038轮：非边配对支付、完整审计与三份LP反证

本轮没有新的Cn上界。当前已归档严格算子界仍R036；另R042新阶上界已在独立clone获审，按独立轮次另封存。本轮只比较原有限见证族的可行支付证书，不把小于1的个案负载称作完整算子常数。

## 机制与有限族的一般输入意义

对同一输入f、水平λ及真实见证V_i=aQi，固定h_i=|V_i|，任意q_i≥0且Σh_iq_i=U=|并Qi|都给λU≤∫|f|Σq_i1Vi≤L(q)||f||1。开轴平行盒Helly使L(q)=max_clique Σq_i，所以重新分配支付可严格改善同一有限实例的上证书。要变成Cn新界仍须控制所有原有限族的sup；仅79实例不提供这个uniform步骤。一般L1及全部正共同半径仍是原紧水平集有限覆盖/内正则性桥接，没有换算子。

R037集中只合并cliques，可能因全部正q在maxh而全停。本轮对每条非边{i,j}保留P=hiqi+hjqj，外部系数不变。完整cliques按含i、含j、二者皆无分三类，外部最高价Ai,Aj,Arest；真实目标是max(Arest,Ai+x,Aj+(P−hi x)/hj)，0≤x≤P/hi。两条斜线交点clip到区间可得精确最小m；全部minimizer区间由每条线≤m解得。其上能量E=Σh q²严格凸，最优x为P/(hi+hj)再clip到该区间。因旧点可行，新pair的L不增；若L平则只接受E严格降。每步完整扫描所有非边（zero-total pair跳过），按(Lnew,ΔE,i,j)排序，最多128次接受。只有完整末次扫描无改善才称pair-coordinate stop，封顶不称收敛。

这些公式对任意正预算也成立；本批固定原witness-volume预算与原f，不能把独立margins自由调整。全程用Fraction与完整正体积clique图，浮点仅用于三个后续LP的提案；所有接收/比较有准确证书。

## 注册身份、实际批次与准确比较

精确79个输入来自R037失败标签的cached q_before、U及原history身份。没有重新跑原atoms、selection、flow或baseline LP；只恢复锁定旧generator的几何/h/complete witness graph，并准确核对Σhq=U、原L和hash。旧R037注册差异仍属于实际输入身份，不因此变成新的独立样本。

实际程序/15动态进程/79任务/128步上限在启动前由phasea_registration.json锁source及input hash。soft理论SPEC在启动后才写成，并只读比对为同一pair与plateau规则；不能把后成SPEC的hash称作事前注册。三个LP候选选择脚本及其三输入在79结果完成后另行冻结，再执行三次，属于有范围的候选认证，不是预注册频率推断。

79/79完成；准确比较62个L下降、17个相等、0个上升。50个full nonedge scan坐标停止，29个在128步上限停止；共3827次接收。A/S0/j6 partial一步越过旧集中零步障碍，其原/新完整分数在compact bundle；展示约3.94191593e−5→3.05180438e−5只作读数。没有从这些小数重造证书，也没有证明该例新支付为全局最优。

旧LP标量interval只作历史诊断：23/56的半gap分类不是79个新LP严证。仅三例重新得到完整有理primal/dual与每条约束，独立审计全部PASS：C/S2/j5 full为非边坐标停，C/S1/j6 full及A/S3/j15 partial是128-step cap。若Mlo≤M≤Mhi，则Klo=U/Mhi≤Lopt≤U/Mlo；三例均准确满足L0−Khi≥10^-4 L0且2Lfinal>L0+Khi，故该有限算法输出未达到半gap目标。严格分数、完整向量与界方向见compute/lp_cert_phasea_*.json及../../audit/round038/lp/audit_lp.json。不能从非边坐标停推断相邻pair或一般多点无法改善，更不能从封顶推断局部最优。

## 独立认证与本机资源

独立Terra审计从有理几何重建严格相交图、完整极大团，验全部3827步非邻、payment、真实maxclique负载及strict lex下降；50个终态额外独立枚举129033条非边，均无改善。79件输入最终hash不变，无图截断或数学失败。审计没有验证每个中间chosen pair在全部候选中全局最优，代码read-only比对另承担执行规则一致性。Python整数位数上限的3次未完成尝试修后只续未完成件，误匹配的registration/summary不是案例，未计数。三个LP独审只检查新原/对偶向量，没有重解旧LP。

生产wall约328.33秒，累计worker CPU约1691.02秒，15个worker PID，每worker峰值RSS约92MB。启动采样15worker各87.7–95.4%CPU；主另实际观察约1425.4%合计，随队列耗尽降为4个忙worker。primary_resource_snapshot.json保留后续一份真实889.8%采样及原stdout，不伪称它是最高/全程占用。随后独审利用尾部闲核最多10进程并行验已完成件。没有填内存、空转或扩大注册任务以制造负载。所有程序按任务/步数/图节点有限上限结束，无总墙时间限制；没有n7数学工作或MATLAB会话。

## 最小复核、来源及结果限制

原79个完整结果约71MB留原本机独立clone不改；其实际raw SHA与最终q canonical SHA均记录。必要数据压为phasea_compact.json.gz约200KB，保存所有几何、h、U、q_before、完整cliques、接收pair序列及停止状态；反复的巨大每步分数由标准库replay_compact.py重建。因此它保留复核所需数据，不是原71MB的逐字节压缩替代品。

在本compute目录执行python3 replay_compact.py，只需标准库且不依赖原raw文件/SciPy，即可重放79件3827步、检查真实L/支付/lex及终q hash，并检查三份LP全部原对偶约束和H严格不等式。此standalone replay由producer编写，不能当作第二作者；独立原件审核另在audit/round038/phasea与lp。可选--terminal会新做50终态scan，本次归档未重复执行，已有独审全scan证据。全文/数据hash由MANIFEST及compute/compact_manifest、archive_checksums锁定。

公开/内部Cn比较均未变化；既未宣称公开发表/全球最强，也未用有限模型cost替代一般算子界。R039另分析允许Σq下降而E暂升的plateau与三cube反例，R040/042转新选择/共同mass及更强截断。下界和七维专项保持停止，goal持续active。论文并行修订，未在本轮声称已替换主PDF。所有成果每轮独立commit，网络未确认时保留待推送。
