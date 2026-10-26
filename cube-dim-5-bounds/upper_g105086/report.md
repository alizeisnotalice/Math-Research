# upper_g105086：固定H1/64的新g105零核资格

新目标21/200000=0.000105：全部180040缓存及97960有效anchor准确扫描后，没有已知H_U<目标的障碍；新target重包96104叶通过、1856仍未证。无新H、physical-gradient或核求值，不把局部通过推为全域g105。070/077及全部旧链字节锁定。

## 输入与准确球诊断

唯一完整输入为070 final/state.json，SHA `8727af764ea088b1c30522f0204d6630ad6e0a4d804c7698c5ed2d646cedf343`，原cache180040、有效叶97960、历史split45600、roots52360。所有球均为同一alpha1模型、epsilon=1/64的H及四个physical gradients；没有消费积分侧epsilon1/512的球。

扫描先登记后执行，PID72679，11.061秒。完整Fraction比较每个H球，cache和active均无H_U<21/200000，且全部有效HL>目标。缓存和有效叶的最小HL/最小HU同属source叶97953、node `[52328,"111011"]`，anchor为

`1,522241/524288,4177921/4194304,33292289/33554432,2097153/268435456`。

闭盒是[127/128,1]×[127/128,1]×[63/64,1]×[0,1/64]。完整原球和准确比较存diagnostics/mapper，HL下舍入为 `53115711/500000000000`，约106.231422微；证书从未消费显示小数。这只是缓存anchor的最小球，不能当作连续H的全域实际最小值。

075曲率与081原值独审来源列在登记/manifest。081可继承生产值真实性及旧g100接入，但g105通过性由本轮另资格，不从g100直接升级。

## 新target的独立证书

075 full frame给Q_full=I4/16。每个有效源叶直接读取原HL/HU及physical-gradient球，以512bit重新包装gradient/H，并重新取log(HL/g105)向下端点。旧g100 Delta及旧margin没有重放。

准确判据是 Phi(v)=Delta105,L + Σ min(gj,L*vj,gj,U*vj) −256 Σ vj²。该函数全局凹，16个物理角点严格正才认证整个闭盒。全部计算使用Fraction；新certificate独立标明target/epsilon/fullQ16以及source叶索引/type/hash。原070叶类型和原证书保留在不可变源及一一mapper中，没有原地改写。

风险排序为HL有理升序，再root/path。48个有限块每块≤2048叶、每64行fsync前缀字节offset和SHA检查点，完整块串联父checkpoint、登记和证书文件hash。无wall alarm、单进程0worker。首2048风险叶短反馈239通过/1809未证，另一个有固定47块上限的进程完成剩余。两程序记录2.973秒/34.652秒，峰值约2.80GB。

## 全结果及准确残余

- 测试97960叶，1567360角点；96104通过、1856未证，覆盖904个roots。
- 原baseQ叶77453、signed81叶16341、fullQ叶4166；本次未证1856全部来自原full_frame_16，前两类全部通过。
- 未证叶深度0/1/2/3/4/5/6分别423/474/443/327/153/28/8。
- 每个残余anchor HL仍严格>105微；失败是闭盒下界不足，不是H反例。

准确集合 `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_g105086/uncertified_boxes.jsonl`，SHA `0cfee52567e81ec39af81176efd7a792c593dfa5000e242aa2bf3477db8daec7`。每行含source叶索引/旧type/hash、准确闭盒、HL/U以及独立新Delta/gradient/fullQ failure margin。盒并集体积13155/16777216（忽略零测边界重叠）；这不是H反例集合或积分误差证书。最小已通过margin有理下舍入 `1759433/250000000000`，完整分数保留。

## 完整新target状态与一一mapper

独立部分状态 `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_g105086/final/state.json`，731,890,653 bytes，SHA `8647a4543ef6a869e747ec06d93f5137cf11feb7b3a64a002995d25a681c4f9a`。顶层target明确21/200000、epsilon1/64；96104 passed使用本轮新证书，1856 queue保持原闭盒，cache180040/roots52360/history splits45600原样保留，inflight/resume为空。历史split注明源g100，仅作为树来源，不充当g105叶资格。

accepted及pending证书各保存原source070 leaf index/type/hash。`state_leaf_mapper.jsonl`逐一将全部97960源叶映到新passed/queue索引，保留原物理盒、anchor、H球及旧type/hash；SHA `af6d5d4f7d5b617c8baee8e839796ddb015d0b8c24d6003d68cdcc89404f310f`。pending_order089.json给确切困难顺序，SHA `68f4ef52ed6dd1360d8487266e9a744a82d4298000677092abd67b026bf0c8c8`。

根表被准确核为32阶全部有序整数四元组。任意1≥y1≥y2≥y3≥y4≥0，取ij=min(31,floor(32yj))即落在相应有序root闭盒，包含0/1边界。二叉闭分裂、prefix antichain及每root准确Kraft和1保证passed+pending完整覆盖，未证区域没有被丢弃。

## 作者资格、来源和后续

verify.py PID73370，19.439秒，复核全部97960新target Delta（768bit log下端交叉比较）、physical-gradient/H精确四比值包络、独立可分端点最小值、最坏角点/通过性、mapper/orderedroots/闭Kraft、48块链及完整残余/status身份。无H/physical-gradient/kernel调用，未重跑g100。原source070 SHA前后相同。

物化独立从source070+新证书构造，并重算canonical文件SHA相同，所有复用cache/roots/history对象保持原内容。最终物化PID73878记录16.890秒。所有扫描、资格、验证及物化PID已OS查无。

089仅准备方案和来源，未派新H。1856当前叶自身已有generic anchor缓存，标准固定选点应先直接复用，不能重算旧H；是否先零H扩展几何，再对未缓存子点派H，由主另确认。新g105仍pending，旧g100资格不受影响。整体I仍优先，不混epsilon，不推Cn/C5。

最小入口为scan.py、requal.py、verify.py、materialize.py；完整新证书块、mapper、准确残余及完整状态留在仓库外。所有路径/hash见manifest。未提交、推送、使用MATLAB或修改其他项目。
