# upper_cover070：有限12进程新H批闭合077剩余前沿

实际完成五个有限新H批，077剩余684叶全部闭合。最终roots52,360、passed97,960、split45,600、queue0、inflight0；没有严格H反例或数值/程序fault。该结果是待主及独审接入全域g的闭合覆盖候选；整体I仍未由本轮支付，不声称C5突破或改动正式常数。

## 唯一输入、冻结对象与登记

完整唯一基态 `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_framecover077/final/state.json`，SHA `c55918ec39ef62629d53d20e1aaed7a5faf454d52212c5897bd16fd15056ea9b`；原passed96,959、split45,283、queue684、cache178,722全部只读复用，未重算任何旧通过点。075理论报告/manifest、078 final接入审计、原046准确API、原3504锥表、077剩余困难顺序均在计算前后核对完全相同SHA。各精确绝对路径与锁列在registration/manifest。

模型为alpha1、16条未单位化body Cauchy方向，逆Fourier因子 $(2\pi)^{-5}$；目标保持 $H_{1/64}>1/10000$。075全frame曲率 $Q_{full}=I_4/16$ 给出尾位移闭盒证书

$$m(v)=\Delta_L+\sum_j\min(g_{j,L}v_j,g_{j,U}v_j)-256\sum_jv_j^2.$$

有向gradient仿射最小值全局凹，再减正二次仍凹，因此准确检查16物理角点足够。新Delta来自缓存H下端与目标之比的向下Arb log端点；gradient来自H及四导数球的有向比值包络。每个角点margin均用Fraction准确计算。

首次684点完全沿冻结077剩余清单的困难顺序。后继split子节点按父缓存H_lower、父fullQ失败margin、root/path升序派发，优先解决低H/低margin区域。闭盒margin失败仅按最长边、最低坐标优先二分成两个闭child，不能靠提高精度修复几何。

## 实际核与有限程序

使用已存在环境 `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python`，Python3.15、python-flint0.9.0。没有安装工具或使用MATLAB。每个有限批先写完整数学/代码/源/候选前缀/资源/停止条件/真实master PID登记；spawn12本地worker初始化后，真实worker PID清单先落盘，再派第一点。所有池完成后close并join，全部worker exit0。

每批最多4096新anchor，本次实际批长684、420、174、32、8。用户最新取消程序wall上限，因此新程序没有alarm或总wall硬截止。每个程序只有一个有限输入前缀，后批由代理另起进程，不存在无限输入controller。

核严格使用锁定046 FirstJetKernel：每个新anchor准确非零且互异的280整数投影，全部3504锥、H及四个物理尾导数。选点仍为盒中心，然后按固定三组generic微扰；不对零/重复节点运行普通DD。仅在非有限、宽度不足或H球跨目标时使用256/384/512数值精度梯级；本次所有1318点在256bit一次完整通过数值质量，没有精度retry、cache命中或partial返回。实际累计4,618,272个锥求和。

worker在点前写可恢复point checkpoint；每个合作分段最多约1秒后按128锥步长返回精确dyadic partial并保存，再继续同点同精度，不遗漏锥，也不构成程序wall截止。本次每点均单段3504锥完成。严格H_U<目标才保存固定epsilon目标反例并停止新派发；真实诊断/fault同点保存并停止。未触发这些停止条件，绝不将其解释为G反证。

## 五批结果

| 批 | master PID | 新anchor | 新passed | 新split | 剩余queue |
|---|---:|---:|---:|---:|---:|
| 1 | 69217 | 684 | 474 | 210 | 420 |
| 2 | 69398 | 420 | 333 | 87 | 174 |
| 3 | 69650 | 174 | 158 | 16 | 32 |
| 4 | 69812 | 32 | 28 | 4 | 8 |
| 5 | 69932 | 8 | 8 | 0 | 0 |
| 合计 | — | 1318 | 1001 | 317 | 0 |

每批summary在最终完整物化之前记录wall，分别16.124、18.948、16.890、14.329、10.360秒；第5程序随后额外执行单次完整物化，因此该行不冒充最后程序总wall。当前不设wall硬截止。每批都是有限输入，末批闭合后未再派核。

全部1318新缓存H_lower严格大于1/10000；最小值有理下舍入 `53115711/500000000000`。1001个新通过叶的最小log margin有理下舍入 `52803633/1000000000000`，严格正。新增节点深度2/3/4/5/6分别664/440/174/32/8。最终cache180,040，旧178,722项全部保留。

## 日志、恢复与完整核验

每批保存小overlay，而非重复复制674MB未变cache。overlay含累积新cache/pass/split、当前queue/hardness/inflight与计数，唯一旧源仍077。单写入者日志hash串联dispatch、received、merge、finish；结果收到即fsync，再合并。五批共3,959事件，所有日志从各自initial_overlay零H回放，产物均与本批state逐字节同SHA。恢复不派核、不盲目重启故障点。

作者另运行verify.py，PID70182，9.026793秒，0 H/gradient/kernel调用。它用交错二进制独立重建物理盒，重现固定generic选点、280投影guard、point checkpoint/任务hash/返回球/缓存的一致性，核全部1318新记录的精确质量宽度与gradient比值包络、512bit log下端比较、21,088个fullQ角点以及通过/失败状态。所有raw attempt记录完整3504锥。

每批和最终前沿均核52,360个root的准确Kraft和1、antichain、closed child分裂、祖先闭合及节点不重复。最终97,960叶等于52,360+45,600，queue/inflight为空。verify.py还仅从不可变077源和最后小overlay重建最终canonical字节SHA，与完整物化文件相同，从而核旧cache、passed/splits及全部旧元数据前缀保留。

OS另实查五master、六十worker、作者verify共66个PID全部不存在，所有生产池均close/join/exit0。具体PID、每批文件hash、各点checkpoint hash清单见batch_index及外部point_manifest。

## 冻结候选与复现

完整物化 `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover070/final/state.json`，680,238,201 bytes，SHA `8727af764ea088b1c30522f0204d6630ad6e0a4d804c7698c5ed2d646cedf343`。顶层保留源元数据，仅cache/passed/splits追加新对象、queue置空，新增task070/status070/fingerprint070/conversion070。fingerprint绑定run registration及最后overlay SHA。

只复核结果、不重新求H：

```bash
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B /Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper/cube-dim-5-bounds/upper_cover070/verify.py
```

只回放某一批journal、不派核：

```bash
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B /Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper/cube-dim-5-bounds/upper_cover070/run.py --work-dir /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover070 --recover-only /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover070/batch0005
```

重新求值必须使用新空scratch并逐次运行run.py的一个有限批，所有原源锁需仍匹配。当前已闭合scratch不会自动启动下一批。

主与外部独审已获冻结完整state和全部链路径。全域G的homogeneity、029转移资格和整体I支付仍由主/独审终核；本轮不把闭合H覆盖直接宣称为C5突破。所有旧锁保持不变，仅写新070及外部scratch070；未提交或推送。
