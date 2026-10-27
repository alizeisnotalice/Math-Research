# 173：172 bands + 153 split 的有限consumer（仅准备）

固定171输入不改：455个EPS64同点donor、35872 pairs、1120 receivers（witness64，其余32）；71个pair曾以同source旧155q失败。15个mass worker各独占30/31个source，join后再15个receiver-owned corner worker；任一时刻最多15 workers。每source仅一次49-band计算，共22295新加权；最多573952新闭角点。无新H/log/gradient/mean/exp/split，所有同目标Delta/gamma从已审155准确复制。旧4096 massonly批不运行。

## 真实字段和单位/模型门

原g球来自151 coefficients `jet_id` 对应的原minimal jet，原H64区间为 `values[0][0:2]`，必须校对model、epsilon=1/64、五坐标、raw_record_sha256、minimal_projection_sha256。master只在真实授权后流读151 jets一次，提取455个原记录到各source-owned运行输入；准备过程只对jet压缩原件作字节SHA，未解码或访问其中的H区间字段。

Donor统一字段 `endpoints[0..72]`，每项 `{t,H:{lower,upper},...}`，73个timekeys须逐字等于锁定registrations的完整列表；其中index24=1/64，index32=1/32，最后1，末尾∞必须另加。point与原g同5坐标准确有理相等，固定第0坐标1，不用置换/齐次缩放、均值/ratio替代H值。

411 donor是range092 compact：锁packet SHA/key/point/完整73time，核`original_raw_sha256`及原global080 raw bytes SHA，并要求原raw端点数组和point与compact完全一致。5个095和39个103是原完整evidence：锁packet SHA/key/point、原run/kernel/registration合成fingerprint。模型绑定同时锁producer registration+manifest和共同bodytable SHA8dd2...，准确限定n5/alpha1/body16/unnormalised inverse Fourier(2pi)^-5；不同历史拼写由明确的冻结stage映射绑定，不凭字段缺失推定模型。所有source byte SHA实际核；174/主仍须审接受这些原始端点的资格来源。

对donor同time H64=[dL,dU]与原g64=[gL,gU]先要求 `max(gL,dL)<=min(gU,dU)`；不交立即保存诊断停止，不clip/不混核。两份原球不改。band分母的合法下端取L64=max(gL,dL)，U仍取原gU。H32仅取donor的真实lower，核0<=L32<gU。

## 新准确公式与16角点

从原73端点中取49个time>=1/64，各相邻[a,b]及完整∞尾：m_U=HU(a)-HL(b)，∞端HL=0；有限w=(2/3)*(a^-6-b^-6)/(a^-4-b^-4)，尾w=2/(3a²)，K_U=sum49(m_U*w)。Qband=K_U/(16*L64)，分母必须lower。新q=min(oldqualified155q,Qband,512/3,256*(7/10-8*L32/(15*gU)))，输出独立rule172及全部49项/来源/上下端方向。保存0<newq<=old155q。

角点phase单source完整16闭角点，原Delta115/gamma不变，允许盒外anchor/负Delta；margin110=margin115+同已存155 Ulog，逐source恒等式Delta115+Ulog=Delta110先核。71个旧pair只在newq严格下降才可重测，否则记skip_oldpair_q_unchanged；首个PASS115后其余pair记skip_after_firstPASS。新q未下降但pair此前未测时，可作新跨锚点对象。

best保留原158及所有实际已测（包括未PASS115）的最大整盒margin110；不做不同source分角点拼接。另110081叶完整继承155+158，遍历所有111201准确leaf索引/原node、52360roots prefix-free+Kraft，全域min/新floor候选与target115 pending分开。正式g/C5/P8合成均由主另审。

## 执行、prefix与授权接口

本包没有运行band/q/角点；只有metadata拷贝、原系数字段匹配、SHA、AST。真实调用 `run173.py --dispatch-reviewed --reviewed-manifest-sha <preparedSHA> --authorization <primaryJSON> --authorization-sha <SHA>`。primaryJSON须绑定reviewed_manifest_sha256、sources455、fixed_pairs35872、workers15，并含三个`{path,sha256,accepted:true}`门：qualified168、qualified172_174、qualified_donor_endpoint_provenance。此为待主确认的机器接口；无此授权不能开任何新band/corner gate。

两phase各先登记15真实PID再开gate；nice+5、无wallcap，各阶段固定scope结束关闭join，故障只留成功prefix。mass每source持久hashchain/checkpoint；corner每128pair checkpoint且每事件fsync。显式--resume-reviewed才可用同输入/代码/源/数值指纹读取成功prefix；不盲重启，不重求H、不重做已记band/corner。原g balls运行投影、455 source-q registry、全部49terms、pair事件、best及退出/资源收据保存外scratch，最后另作最小包。预计总RSS约1.5GiB以内，实际CPU/RSS由每worker/master记录；未设总CPU/墙钟限制。
