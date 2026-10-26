# 128 授权零 H 有限批实际结果

主已完整审阅冻结脚本、来源资格及新证明后，授权唯一 **11330 个 child / 181280 个闭角点** 批。授权原件：

`/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds/round014/primary_zeroH128_authorization.json`

SHA `eff9c488e02b4759a64035b1dd15533fe4d7a578c591b37929577ab42b290a9a`。本文件是授权后执行记录；`report.md`、`source_qualification.json` 的“尚未实际测试”是保留的授权前冻结状态。

## 结果与资源

实际 evaluator PID **38496**，exit **0**，单进程，无 pool，约 **5.237 秒**。完整 finite 对象全部完成：

- **5903 通过，5427 待证**；检查 181280 个闭角点。
- 0 新 H、物理梯度、核、log、exp、split；没有新空间节点或后续派发。
- 最小严格通过余量约 **0.000016068982340623138**，child `[51199,"1001"]`。准确 Fraction 保存在 `results.json`。
- 待证余量范围约 **−0.21965780407200233 至 −0.000007473619023601817**，分别在 child `[51516,"1"]` 与 `[51941,"1101"]`。

所有来源、冻结脚本、准备包、qualification、授权原件和完整 52 chunk g105 packet 在执行前后逐文件 SHA 校验一致。wrapper 的 before/after 锁列表相同，成功输出已保存；没有故障重跑或覆盖。

输入仍是同一 H128 父锚点保存的 Delta/四个物理 log-gradient；只用该锚点的 H128 upper 和旧 H64 g105 整叶下界改为
\[
q_U(c)=1024-768\frac{21/200000}{U_{128}(c)}.
\]
每个成功 child 均由同一锚点的全部 16 闭角点严格通过；没有按角点切换锚点。

## 待证诊断的准确范围

5427 个待证 child 分布于 **1369 个 root**，path 深度 1–6；准确物理体积和已保存在结果 JSON。其本次实际 q 系数约为 **347.2167–554.4776**。5665 对 sibling 中，2756 对两者通过、391 对恰一者通过、2518 对均未通过。

这只说明本次真实新 q 的单父锚点整盒证书尚未分离；不能据此断言 H128 或 G 小于目标。最坏 margin 的节点与逐一来源已冻结到外置 `pending_gap_order.jsonl`，没有重新运行其 H、旧 fullQ 或追加角点测试。

若全部既有与本次共享资格最终通过，候选前沿可有 **104482 通过、5427 待证、57549 split、109909 活叶**。11330 个 child 的新结论本身不改变树，全部 5903 成功叶与 5427 待证叶仍覆盖原 11330 child 前沿。完整树及来源独立审计由 130 负责；123 raw/Delta 的独立数值与主接纳仍待完成。

**全域 g110、C5、I 与全局 S 均未更新。** 没有新 H 或细分授权，不继续尝试剩余盒。

## 冻结原件与复现入口

外置唯一运行：

`/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_curvature128/evaluation001/`

|文件|SHA-256|
|---|---|
|registration.json|`a701dfa9615f7992ec204aad33ceb71bfd5d3fbf29d0a2e40b82f2aea04127e2`|
|candidate_results.jsonl|`765200e4186f5f519691d1516fec85bf1aaafa10c54f4ff87b33d7aa0ae36d67`|
|diagnostics.json|`e2449c1efae0905868492a0006e491987524f85bf28ebb5db4a596e7e7787bdc`|
|pending_gap_order.jsonl|`d24de9b7a86155fdbd17e04bdfa9d0bada64973f11d10a5511931347645dfbfd`|

完整来源前后锁、actual PID 与退出码位于外置相邻 `dispatch001/`。精简归档保存源码、原来源资格、授权/运行登记与结果摘要；不复制约 761 MB seed、123 raw 或全量角点派生数据。

实际执行命令如下（记录用途；已完成，不应覆写 `evaluation001`）：

```sh
python3 -B /Users/zhengzhihao/Desktop/T/Math-Research-cube5-curvature128-sol/cube-dim-5-bounds/upper_curvature128/evaluate_zeroH.py \
  --reviewed-preparation-sha 53621408cb415fede12bba827cac2b57b87a9ba678ec9abef2e583021f289c67 \
  --output /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_curvature128/evaluation001
```

本批的独立复核可直接读取冻结候选与结果，按完整 source locks 重建每个 child 的 q 与 16 个闭角点。若另行重放 evaluator，必须使用全新输出目录并保留本次 immutable 运行；无需导入 H 数值核。

未提交、未推送，未修改任何他人文件；未使用 MATLAB，未新/续派 Astra 或其他代理。
