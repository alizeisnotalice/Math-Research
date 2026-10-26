# upper_mixed098：194残余叶的H1/128准备方案（尚未求值）

089已准确排除固定H1/64>105微目标；完整G未被反证，旧g100资格保留。本阶段只完成准备，所有H128/physical-gradient/kernel/log求值数为0；等待主与099审阅后另授权首194新H，不自动派发。

## 已冻结来源

- 086完整部分g105状态：SHA8647a4543ef6a869e747ec06d93f5137cf11feb7b3a64a002995d25a681c4f9a。
- 089 counter终局overlay：SHAba9ba3ac4f8e7b6788ed75005643dae86ed4e2f139c3a11bf4f12a986b491e74。
- 089 terminal_manifest：SHAc8b639b5ff46b36a602c1f38f2fc291d76a7b87a9b76218a6c808b9f875012c7。
- 075 report：SHA3308d1e5c660cd11d73d47a013acf794bd1f596afdb5a113573a5bfd7481f225，曲率对每个正epsilon成立。
- 090独立核值/结构/新g105资格report：SHAfb83d5017e162bf593d199036939f88a4b78bccb0bdd8f5f9667a133636731ff，manifest8a61660e74a83c8a25e7a58e9a72385daf3e620d1a8727fa93d96cf900c1aaab。四个H64障碍已独立确认。

所有旧source与终局只读，未修改046/070/077/086/089/090。新输入仅由这些冻结对象转换，不重包原H64基线。

## 完整mixed-G初态

`/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_mixed098/initial_mixed_state.json`，SHA `19294a9bd493524fdc23ec2abe2bb3e70577f04b762b08d4fc57e78af5377480`，764,441,089 bytes。

- 100804通过叶保留原certificate字典及node，新增显式leaf_epsilon=1/64、component=H_epsilon、lower_bound_object=G说明。
- 旧190 queue及4 retained counter inflight完整映到194新queue；新inflight为空，4个旧counter原球仍保留。
- 48638历史split及52360根保持，accepted100804+pending194=100998=roots52360+splits48638。所有root闭Kraft、prefix antichain、祖先/closed-child均核过。
- accepted_H64_mapper.json给全部通过叶来源及旧leaf/certificate SHA；pending_transfer_mapper.json给190源queue索引与4原inflight ordinal的新queue索引，包含物理闭盒。没有遗漏counter区域。

cache唯一键为(canonical epsilon,canonical point)，实现为cache_by_epsilon的分离字典。H64桶184796条，含全部旧180040、新正常4752及4个有效但不达105微的H64观测；H128桶为空。新lookup严格只查请求epsilon桶，无flat-cache fallback，不能把同point的H64值/gradient充当H128。

## 待生产的首194请求

首194闭盒全部改用H1/128，仅本叶epsilon变化；共同完整G、target21/200000、body方向/权重/Fourier归一化均不变。现有EPS64通过叶不求新H、不重算gradient或Delta。

所有194固定候选已作整数generic280预检，H128 cache hit=0。四点虽存在H64观测，pair key不同，仍是194个新的H128 anchor。execution_tasks_194.json绑定每个node、精确坐标、epsilon128 pair key、候选/投影检查、专用point checkpoint路径；四个旧H64球仅供单调性一致性guard，不替代新H128计算。

首批预算最多194新anchor、12本地worker，无wall cap。启动前另登记真实master/worker PID和source/code/table/request hashes，有限194输入完成、严格H128目标障碍、真实数值/程序fault或STOP即停止新派。已派任务保留完整point/precision/cone前缀，池正常close/join；不设无限controller。

## 独立epsilon128源码与可审核diff

`epsilon_kernel.py`独立复制046的kernel段，锁定为H1/128。`kernel.diff`准确显示仅如下变化：默认及合法epsilon约束1/64→1/128；时间参数明确改为scaled=144D/128，并assert它等于144D*epsilon。原body vertices280、simplices3504、五分量DD递推、D^4/(144*(2pi)^5)归一化与物理derivative slopes均不变。046原文件SHA仍df3fb4165aea29e4aa3d1fb7ecacdb9cc5ad4a531bc03e3cc43a0820671fe028。

kernel SHA `128d98c94e55632336edafe024aa3a27c1486898670a52be89cb39cc22ed6cc0`；cache_model.py显式EPS pair key。jet_job.py仅定义待调用worker，任何import/准备均不求H；记录epsilon128及kernel SHA，256/384/512只为数值质量/目标跨球retry，绝不因闭盒margin失败升级精度。

certificate.py仅定义待调用新H128判据：075 Q_full=I4/(16epsilon²)，本叶epsilon1/128给准确罚项系数1024；Delta来自新H128_L/target的向下log，gradient包络只来自新H128四个physical gradients/H128，同盒16角点Fraction最小值严格正才认证。每新叶显式记录epsilon128及G≥H128推论。未执行任何该证书/新log调用。

所有Python源码AST解析已通过，prepare.py实际运行PID86050、13.317743秒，OS已查无，未初始化pool或执行任何kernel.evaluate。

## 理论与放行边界

共同模型中G≥H_epsilon逐点，故旧H64通过叶与新H128通过叶都给同一个G的下界；不能把H64球换标签成H128，也不宣称混合H拥有同epsilon的齐次性。全域尺度推广只能对同一完整G使用齐次与cube symmetry，须经099/主终审。新G105目前仍pending194，未称Cn/C5/I突破。

本准备文件、source SHA、完整mapper、typed-cache和kernel diff供主审。没有生产H128授权，也没有自动启动代码；主审核后再建立有限首批的实际运行登记与pool调用。
