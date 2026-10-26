# 145：next234 零 H 跨锚点有限独立审计

唯一授权阶段 `next234` 已完成独立审计：108 个 first-PASS 证书的全部 1728 个整数角点与生产结果严格匹配；3904 条日志、2346 个选择来源及完整 52360 根 forest 均通过。该阶段最终为 111075 个 passed、126 个 pending、58841 个 split。此结论是本阶段的组件资格，`globalg110_qualified=false`、S credit 为零；没有由该 partial frontier 推导新的全域结论或 C5 数值。

## 实际执行与准确范围

主审授权的准确计划为 `audit_plan_ready.json`，SHA256：

`81d01037fe1412cb38cac481c644aeb80c6b430e2e2a6f4641072cad985c1e27`。

外授权同字节副本 `primary_authorization.json` 的 SHA256 为：

`856764bf56f5dbd2b06f471d367a323bcd1c4cf6d14ec00e079d52fd7dd4d433`。

唯一实际审计 PID 为 **64149**，退出码 **0**，完整 subprocess wall time 为 **55.18057820800459 秒**，stage body 为 **52.50235324999085 秒**。新排他输出为 `run001/`，注册和退出回执分别为 `execution001_registration.json`、`execution001_receipt.json`；引擎自己的运行前注册位于 `run001/execution_registration.json` 和 `run001/next234/run_registration.json`。控制台为空且无异常。没有重启、追加阶段或复跑成功角点。

本次执行为单 worker、无 wall cap，只执行 next234 的 108 证 × 16 角点。新 H、log、gradient、kernel 调用、失败角点重算和新 split 均为零；137、141 的已过角点未重算。准备回执和原未就绪配置均保留原字节，明确与此次实际结果分开。

## 通过的数值与协议门

108 个通过证中，96 个使用 EPS64 的 q=256，12 个使用 EPS128 的已冻结 q 规则；108 个 anchor 都位于目标 child 外。仍是每证一个 anchor、一个 epsilon、同一 raw H128 记录的原 Delta 和整组 directed gradient，通过全部 16 个角点取得正 margin，不存在跨 anchor 或跨 epsilon 拼接。

整数核 `exact137.py` 与 137/141 成功版本逐字节相同，与 130 的算术实现独立：先将整证全部项清至统一正分母，再用整数分子核全部 16 个角点、严格正性、最小 margin 与字典次序最坏角点。新 H128 q 直接按正分母的整数等式核验 `1024 - 768*(21/200000)/H128_U`。实际 108 条整数结果位于 `run001/next234/winning_integer_checks.jsonl`，SHA256：

`95474f8de1182e305feaa02c1ec8b0eb444c0426d3450b79a21025eff62f20c6`。

15 个生产 shard 的完整 journal 共 3904 条，分为 2418 条 evaluated 与 1486 条 skip。独审核验每条 hash chain、seq、offset、checkpoint/fingerprint、全 pair index 一次消费、first-PASS 后的 skip 语义、每 child 唯一首次通过证以及 residual 的保存失败证和 gap order；没有重新计算失败角点。

## 来源身份与数值资格

实际来源库为 117284 个唯一对象：旧 131 库 111331、138 新库 4274 和 143 新库 1679。全部池身份、完整计数、source index 和 pair 池标签均检查；从中实际绑定的 2346 个选择来源，包含 138 origin 的 476 个和 143 origin 的 567 个，其余为旧来源。旧来源即使未被 137/141 选择，也直接回绑原 105、107、123 certificate 与原 raw record，不以旧审计的 selected 标签代替原身份。

每个来源均核 source/core/object hash、原 certificate hash、raw record hash、point、H interval、原 Delta/gradient 和 penalty。使用更小 H128 q 的来源，均回读原 098 完整 g105 的 EPS64 closed leaf，核 target `21/200000`、严格正 pass margin、原叶身份、祖先路径和 point 闭包含；同时把新 ancestry proof 绑定到原生产 state 的准确 SHA。原 schema 的 inner certificate 若没有冗余 target，仍要求外层 target/EPS64/严格 pass 和完整主审 g105 来源全部成立。

143 的数值资格直接锁定真实最终 manifest

`7cc2db15c61472fda0898499e3ecb5900f032f386e30ada64aabb54f80620db0`，

其成功 final check、structure、source lock 和原生产 state 分别为

`015981b7eae3d510e20dfc39b1b581f885d39d6dc96d6235e0fa3ea3c67fc655`、

`0dd4df9d929dac977012ae51a599e67188c94784d13bdb7adb0d18fde0cd5a39`、

`c61c957ece2640be306ac0f73f3fbb6fac427e009ecb59236964ae3e4d0c415d`、

`ad5b35f4aad937086a164047209867b63e5cc67e547b0d6bbc4eb62dcec79882`。

它已独审全部 1679 点、8395 分量、26864 角点及 5038 事件；145 不重跑这些值。138、123、130、137、141 等继承资格、原源文件与完整生产输入都由 ready plan 的 **6565 个实际输入 SHA** 明列，完整清单以该计划为准。引擎实际执行前后逐项核验，均未改变。

## 全树和边界

起始已冻结 frontier 的 1679 个 pending 被 143 的 1679 个 parent 恰好一次消费，1562 个新 passed 和 117 个新 split 产生 next234 的全部 234 个 child。每个 child 都核 direct parent、原 parent source、精确二等分盒、正确分割坐标和闭边界。

加入本次 108 个通过 child 后，全部 52360 根的 heap 父子关系和每个 split 的两个 child 均闭合；所有 leaves 无重复且为 antichain，每根 Kraft 和精确等于一。实际最终计数为 111075 passed + 126 pending = 111201 leaves，58841 split，满足 `leaves = 52360 + splits`。126 个 residual 与日志和冻结 gap-order 文件逐项相符。

## 冻结与复现

`stage_list_ready.json` 只含明确命名的 next234；没有读取或发现未知未来阶段。`audit145.py` 允许以后对另一个主审后的显式有限清单使用同一核心，但本次资格只覆盖准确 ready plan。未来阶段不属于本次授权或结果。

全部代码、配置、授权、准备材料、实际注册/退出、完整整数结果和结果摘要列于 `manifest.json` 的精简白名单。`final_check.json` 记录最终源锁、代码锁、旧冻结件及退出核验。旧 141 manifest 始终为

`4d312e5db169cd34d3d21afc209fe9ecdbab3f6f3e407a2b78b587cda689fdbe`。

本轮没有 Astra、MATLAB、下界工作、提交或推送。主代理保留最终接纳和完整证明责任。
