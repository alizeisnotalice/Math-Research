# upper_g105089：固定H1/64的g105方案严格障碍与完整终局

固定目标H1/64>21/200000在第9有限批被四个准确H上球严格排除。本链已停止派发，全部worker正常close/join退出；这不反证完整G，也不撤销已合格g100结果。

## 输入与生产范围

唯一源086部分g105状态SHA `8647a4543ef6a869e747ec06d93f5137cf11feb7b3a64a002995d25a681c4f9a`，96104个epsilon1/64通过叶、1856待证原闭盒、cache180040。先直接消费1856已经准确失败的g105 fullQ证书，零H预split成3712 closed children，不对同盒同anchor重包。预split映射、完整root Kraft及旧对象恒定先核后派。

模型始终同一alpha1/16条未单位化body方向、Fourier因子(2pi)^-5；target始终21/200000、epsilon始终1/64。旧cache精确命中直接复用，所有新点严格280非零互异整数投影guard、全部3504锥和H+四physical gradients。256/384/512仅作数值retry；本次全部4756个实际新H记录在256bit完成，每条3504锥完整且numeric_quality=true。

每程序至多4096新anchor、有限登记当时准确queue，12worker、无wall deadline/alarm，每批完成close/join。没有无限输入controller。所有source070/077/086/075/081及历史stage字节锁保持。

## 有限批结果

| 批 | 登记/派发 | 新H | cache hit | 新pass | 新split | 待证状态 |
|---|---:|---:|---:|---:|---:|---|
|1|3712/3712|2588|1124|3426|286|572 queue|
|2|572/572|566|6|299|273|546 queue|
|3|546/546|546|0|388|158|316 queue|
|4|316/316|316|0|188|128|256 queue|
|5|256/256|256|0|168|88|176 queue|
|6|176/176|176|0|112|64|128 queue|
|7|128/128|128|0|44|84|168 queue|
|8|168/168|168|0|75|93|186 queue|
|9|186/12|12|0|0|8|190 queue+4 counter inflight|

累计4756实际新H观测，1130 cache复用；其中4752新值在目标上方正常入positive cache，另4准确反例值完整保留在raw及counter inflight。正cache逻辑总184792，计入4个有效原H64反例球后全H64观测184796。新增通过叶4700，总passed100804；含零H预split1856，总split48638。完整pending194=190 queue+4保留counter节点，未遗漏被反例点所在区域。

## 四个固定H64目标障碍

最强点为

`(1,2095105/2097152,16760833/16777216,134086657/134217728,2097153/1073741824)`。

它位于normalized cube face且满足全部generic guard，准确H上端约0.00010496009654751766577672612549，小于0.000105。判定来自完整有理上端与21/200000之差，显示小数未作为证书。其余三个点同前四坐标，最后坐标分别6291457、10485761、14680065除以1073741824，也均H_U<目标。全部四球来自完整3504锥、256bit、numeric_quality=true。

`fixed_H64_obstruction.json`及外部同名batch0009文件保存完整五分量球、point checkpoint/来源SHA和严格差值。只否定固定截时1/64的H目标；G≥H1/64没有反向上界，因此不能称G或C5反例。g100原闭覆盖及其075/029接入不变。

## 保存、回放和作者核验

第1..8批共17630日志事件，所有source/parent hash、真实PID、point/cone checkpoint、完整cache、16角点fullQ资格及闭root/Kraft由作者逐批零H复核后才接下一批；原字节manifest均保持。第9批仅派首12个点，在观察到严格H_U<目标后停止新派，等待已派任务完成，12worker均正常close/join exit0。最终总17667事件，新16角点证书5882条，共94112角点；4个反例不伪装成正资格。

作者verify_stop.py另用零H过程核第9批37事件并准确回放到同SHA，核全部12原球的物理盒/280投影/质量/3504锥来源，四个严格有理比较，以及8个正常点的768bit Delta和独立可分端点判据。精确root数为52360，每root Kraft=1；passed100804+pending194=52360+split48638=100998。作者stop PID83864记录2.754秒。所有9个master、108个worker及stop核验PID均OS查无。

冻结终局overlay `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_g105089/batch0009/state.json` SHA `ba9ba3ac4f8e7b6788ed75005643dae86ed4e2f139c3a11bf4f12a986b491e74`；journal SHA `f048b91154bb667701d877ead87767b15222b8207ded654a1794da0353cf1654`；该状态190 queue+4 retained counter inflight是数学待证前沿，不代表仍有活worker。source086+完整累计overlay足以重建全部旧cache及新值/树；没有逐批复制731MB旧cache。

## 终止范围与新098准备

本固定H64方案在严格障碍条件终止，不因研究总时限停止。主已选新098仅准备：保留100804个H64通过叶作为G逐叶下界，194完整转新queue，残余改用H1/128。H128必须独立epsilon cache桶，旧H64同point不能覆盖或充当新H128；新kernel独立源码将明确144D/128和Q罚项1024。当前没有任何H128求值，等待主与099理论终审。

090正在独立复核本链新值（含四障碍），作者结果不替代该外部资格。完整G推广只依赖共同模型的G齐次与G≥H_eps逐叶，不能宣称混合截时H同eps齐次。整体I仍缺，本轮不更新Cn/C5。未提交、推送、使用MATLAB或修改他人项目。
