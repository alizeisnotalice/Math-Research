# 094 wholephase 有限审计协议

本阶段只消费最终冻结的 092 manifest 及其全部 14 块、56139 节点、每节点 73 个旧 H 端点；原 52360 roots、320 internal parents、55819 physical leaves 不改变。数学计算开始前另写 whole_scope.json，锁本协议、代码、catalog 结果、091/092/088/085/076/075/087 来源和准确 initial/final state SHA。旧 catalog 登记不改。固定 14 块（前13块各4096，最后2891），8 个本机 Python worker，逐块持久 checkpoint；无 wallcap。成功块只在 code/scope/block SHA 完全相同时恢复复用，不重复已成功全量基线。

1. 逐项核所有冻结 manifest 和输入 hash，逐节点将 092 old fields 与已接受 088 精确记录、记录 SHA、初始 state 逐项对应。旧 M2 的全量独立资格由 088 冻结证书消费；094 不重跑其成功 moment 基线。对每个真正 anchor/box 独立重建 h∞ 和 ||μ||∞，检查 u1=0。
2. 读取并核全部 original raw 字节 SHA 与 endpoint-only packet SHA，准确比较 point/endpoints/indices/weight/level，确认原端点没有换 anchor、模型或时间。这是输入身份审计，不是 H 数值的独立重生产。原 H 区间与 Arb 内核仍属共享接受输入。
3. 由相邻 H 区间作正 band 差（最后以0支付∞），mL>0。每个 x=h∞/a 消费已经 094 整数审计通过的 catalog RU。计算 mL/RU、mU RU 及 mU min(RU,1+||μ||∞/a+(25/48)M2 RU/a²)。每项只用整数/Fraction，以分母2^256向下/向上取整并相加，确定准确的2^-160 floor/ceil bin；若某 bin 区间有歧义，仅该节点回退准确 Fraction 总和，不新评 H 或比值。允许生产 Arb 端点比真实准确 bin 多付一格，记录是否逐一相等；危险方向一律失败。
4. 检查 L=max(oldL,newL)、U=min(oldU,3,newU)、M=min(oldM,U,newM) 的安全方向及 0<L≤M≤U。用 li^5 与 ui^5 的纯整数四次根，向外证 z^(5/4) 凸弦；检查 selfL=max(oldselfL,weight·newlower)、selfU=min(oldselfU,weight·newchord)。不把小于新 chord 的旧 selfU 当成错误。
5. 逐字核 initial/final 所有边、roots、physical frontier 不变，精确核 permutation weights、parent 16倍权重和 physical volume=1。空间闭覆盖沿用 087 已证同一 frozen 082 geometry；094 不重复 087 所有边界案例。
6. 用递归 DFS（生产用 level 逆序）独立重建 14 个完整前缀，每个前缀明确检查所有未更新节点仍支付旧 selfL/U，再核 effective L=max(selfL,sum children L)、U=min(selfU,sum children U)，根总和与 receipt 完全一致，UH 单调不增。末状态逐节点与冻结 final state 完全比较。
7. 只用整数有理幂核 B>E、(B-E)^5≥UH^4、Iupper^4≥B^5、costupper^4≥(5B/g)^5。g=1/10000 和 E076 是冻结同模型接受输入，不重审 g；E 是完整短时间 L^(5/4) 范数，不能直接加进 I，也不再加064 tube。准确测试 Iupper<70301737/10^12 及 5^5 B^5<32^4 g^5。未过32门则不能称有新 C5/Cn。

全部新增数学计算为本机标准库 Python；无 MATLAB、Arb、log、exp、asinh、新 H/stencil/jet/锥。有限对象和恢复策略先注册；失败保留 exact rank/key 与堆栈并立即报告。只写 audit_range094，不改生产目录，不提交。
