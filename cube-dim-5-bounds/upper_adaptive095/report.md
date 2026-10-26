# 095：最高绝对 gap 的首 1024 父盒同步细分

已完成且作者独立读回通过；仅执行登记的这一批，不自动追加。未改变正式 C5=32。合法 g=1/10000 已由主证明链接受，g105 未消费。

## 范围与数学接口

初始源为已由主和094接受的092，manifest `05cd878504800f4a66673d9a408666c6480b92a1db5828de01f6ba5ec0cd801c`，其56139个历史节点、55819物理叶盒完全保留。选取当前物理叶盒准确加权 selfU−selfL 最大的1024个，依 level/indices 打破平局；名单和全部唯一子盒在新H前冻结于 preflight_registration.json。每个父盒四轴同时二分，child=sort_nonincreasing(2i+bits)，唯一parent=floor(child/2)，子orbit权之和=16父权。结果13406个唯一子盒（L6=11476、L7=1914、L8=16）。不将单轴二分套用旧orbit权。

核固定 n5 alpha1 body16、eps=1/512、q=5/4、384bit、73个端点；3504真锥的0阶DD全部有符号聚合后取球，使用017的真实280ray资格。每点 exact280投影互异且锚点严格位于真实满盒内部；未使用H的eps齐次，不计算旧jet。新空间子盒保持真实满盒，offset/体积准确。

range仅使用091新chain ratio的192bit向外端；97项catalog中保留092已有89项、仅新增8项。mean使用085面内25/48·Emax|u|²余项和真实非对称 μ；凸弦使用真实range+mean。global Hcap=3合法全域截断。新L下端、U/mean/chord上端在2^-160上向外量化。完整父盒所有子证书成功前仍支付092旧self界；成功后所有历史层级继续取 effU=min(selfU,sumchild effU)、effL=max(selfL,sumchild effL)。物理覆盖与数值层级分别保存；不因childsum较差丢弃缓存。本批1024个childsum都严格改善对应父self上端，且都实际降低root上端。

## 查重与真实执行

完整旧080+082 H73注册中0命中。H64注册采用089接受的180040 base+4744累计新cache，32个 exact(point,time=1/64)命中直接复用；pending9全部12 raw point keys与095新点精确交为空。1/64与1/128键区分。g105目标失败不作为H值失败，也不作为全局g105合格。

主PID 87715；12实际worker [87726, 87727, 87728, 87729, 87730, 87731, 87732, 87733, 87734, 87735, 87736, 87737]；wall 320.582432s。新增H 978606，复用H64 32，新增纯0阶锥DD 46974624，0jet；全部worker正常close/join。没有wall alarm/timeout，工作量有限。锥prefix每512保存，每个完成端点保存，完成点按同fingerprint复用；真实诊断停止派发并保存前缀，本批无故障。

累计点内计时：stencil 603.926s、endpoint 3081.485s、envelope 58.040s；这些是多worker累计时间。单worker最高RSS 43286528 bytes。完整raw/stencil/endpoint仅在外scratch的 evidence/，不把全cache传worker。

## 全域结果与32门

最终数值节点69545，物理frontier68201，每次完整父替代的加权体积准确为1。共同量化分母2^160：

- LH = 38413887677479752724438204644362255442485903/1461501637330902918203684832716283019655932542976。
- UH = 75063710980782909270112604401565474168288679/1461501637330902918203684832716283019655932542976。

LH仅为截断H的层级积分下界，不能当作实际完整I的下界。固定076完整短时E，使用整个F直接积分，不额外添加064tube。B=UH^(4/5)+E，I≤B^(5/4)。准确有理上端与所有幂门保存在 run_receipt.json/summary.json。

- 原092 I上端约 6.01026698256809e-05；新上端约 5.470251657598456e-05。
- 可引用简单有理 cap：I ≤ 54702517/1000000000000，严格优于原cap6010267/10^11。
- 对内部有理上端下降约 8.98488081%。
- g100 cost上端约 40.8996707458，exact门5^5 B_U^5 <32^4 g^5为false；正式C5仍32。
- 最大剩余准确weighted self gap位于 `5:29,24,23,20`，gap=592662956803590025831955495223660946209/182687704666362864775460604089535377456991567872。top64详见 summary.json；这只决定后续候选排序，没有执行下一批。

## 核验与证据边界

作者独立读回PID 89178、161.472293s：所有13406点使用与生产拉格朗日权不同的递归DD路径核全3504锥，校验Σw=Σwk²=0及四阶det矩；全978638个保存band正性、32个缓存球精确同一、真实offset/h/μ、条件最大坐标积分计算M2、97项重标度224bit整数根、所有224bit四次根凸弦、1024 hash-chain完整父前缀的层级和volume1、最终全域B/I/cost幂门全部PASS。0新H/log重算，浮点只用于展示；原Arb log区间生成与017/075/085/091解析来源仍是明确共享输入。

verify_portable.py由080/082/088/092量化证书重建初始数值图，再独立整数检验本批全部轨道、凸弦、最终层级/体积/预算；PID 89342，2.241839s，PASS，完全不读外scratch。它不单独证明保存的H球，完整raw作者检查另见 verification.json。

本次repo仅保存源码、有限登记/小输入、全部量化子证书、物理/数值层级摘要、数学hash、少量端点样本及必要报告。大raw位于 `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_adaptive095/evidence/`；final state位于同scratch的state.json.gz，SHA `1a91e0d8ddbe69f7ae3a203547b9ddf61ee233fab716b7fd0f2c1c1f67b07c09`。parent_events.jsonl.gz冻结1024事件链。所有旧source/登记/manifest字节未修改；未提交、未用MATLAB。

## 复现

在该clone根使用已有venv Python运行 `-B cube-dim-5-bounds/upper_adaptive095/verify_portable.py` 可零H核账本。保存外scratch时运行 verify.py可零H核全raw。显式单点从零重建入口为 reconstruct.py --key KEY --output 新目录；它仅在明确手动调用时重新计算指定数学对象并比较确定性数学SHA，本轮没有调用或重跑成功点。run.py无参生产已完成；不要再次启动，后续批必须另行登记并由主选择。
