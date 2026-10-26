# 167：完整 I 输入的便携化方案，等待实际导出授权

本项已完成来源、存储和固定导出代码准备，未大批导出、未运行数学 checker，未重求 H/R/mean/factor/旧 power，也未合成新的 I/C5。156 仍是待 161 与主接纳的候选；正式 I150/g192 状态不由本项改变。主代理保留最终导出授权、源区间资格及完整 I 证明责任。

## 来源与可删除的继承链

121 的全量独审实际消费共享的旧 H 输入并重审 112 新 mean/chord/tree；其 `whole_protocol.md`、`check_whole_v2.py` 和结果明确不重审旧 range/geometry 数学。既有 `verify_portable.py` 主要闭合整数账本与凸弦，并不因名字含 portable 就补齐所有 H→range 来源。

095/103 的 kernel 明确从**同一个 anchor 的完整 73 个 H_t 端点**、091 新 R 上端、真实 offset 几何和接收面 cap3 构造 L/U/M。没有额外核角点或 jet；旧 B² 的 16 角点是几何二次式计算，不是额外 H 求值。

092 起初还与 080/082→088 的旧 range 取 max/min，并与旧 self 取 min。因此仅凭它的报告或 121 的 PASS 不能直接删除这些旧来源。本项实际读取所有 56139 个 092 结果字段，检查选择分支：

- 全部 56139 个新 rangeL 严格大于继承 rangeL。
- 56138 个新 rangeU 严格小于继承 rangeU。
- 唯一 rangeU、mean、selfU 与旧值相等的节点均为 `5:31,31,31,31`。

故对其余 56138 个 092 节点，原代码的严格 max/min 分支证明 L/U 来自新 091 链。准确解释为：原 384 位有向累加下端先向下 floor 到 `2^-160`，上端先向上 ceil 再与 3 取 min；不能声称它们必然等于理想有理和的最近 floor/ceil。便携 checker 应直接验证保存 L 不大于准确合法下和、保存 U 不小于准确合法上和或由 cap3 支付，不必重现 Arb 中间舍入。

特殊节点的旧 range/self 资格可由真实 cut 排除。新 cut metadata 程序只比较保存 mean/self/eff 字段并做局部已存 child-eff 整数和，没有重算 mean、power 或全局 UH。实际 PID 86244：全部 3392 个内部节点均选 childsum（相等时选 children），末端恰为 94580 个物理叶；所选 self 全属于 97971 个 112 严格新弦节点。特殊根有 5 个子节点，其新 childsum/new eff 为

`26143870058526986308648060983955316061135`，严格小于旧 selfU `5503018354736479757810956912611615042231038`。

因此这份 cut 不使用特殊根的 legacy self，也不使用其旧 range；其余节点全部 092/095/103 range 可退回同 anchor 091 新链。原 080/082 exp/asinh R、旧角点 B² 和旧 self 数值不再是本 packet 的 active 算术来源，无须复制其大账本。

156 的 97907 个新 mean 与复用 150 的 64 个 mean 均在保存字段上严格小于 m112，所以当前 mean 确为新 band 和的有向上端，继承 mean 的 min 分支不 active。保留 m112 作为原弦的参考横坐标与 trace 字段即可；未来读者可直接证明保存 self112U 不小于同 L/U/m112 的向上弦，再证明当前 newself 不小于同 L/U/newmean 的弦。无需再证明旧 Mold 是真实均值上包。112 全 97971 个 strict self 也排除了更早 self 的 min 选择；唯一 nonactive 节点已由 cut 排除。

## cap3 的同模型资格

132 已接纳的准确 portable cap 是

`M132=488568596392024500738533860719/250000000000000000000000000000 < 3`，

它约束同模型、同 ε=1/512、整个闭正面 `F={(1,y):y∈[0,1]^4}` 的 H512，139 已明确消费该资格。所有接收 dyadic 盒均在 F 内，故 Ui 的 cap3 合法。这个 cap **只约束接收盒的 H512**；并不施加到盒外 anchor 或 anchor-to-box 线段。线段比值与 Taylor 仍由 085/091 全空间定理支付。包中保存 132 原准确组件、证明及原数值输出，未重优化 M 或重算其 log/root。

## 精确计数与体积估计

只读清点 PID 83983 实际 stat 原 97972 个 packet，总 gzip 字节为 **1,540,460,673**：092 56139 个/465506822 B，095 13406 个/343560282 B，103 28427 个/731393569 B。共 7151956 个 H 区间、完整 73 时间端点；不是 73 个全局 H 数值。

固定导出排除唯一特殊根，保留 **97971 个 packet、7151883 个区间、14303766 个有理端点字段**。来源计数 092=56138、095=13406、103=28427。按固定 id 模 15 分片，6 片各 6532 条、9 片各 6531 条。完整原树仍保留 97972 节点、52360 根、3392 内部节点、94580 叶；不通过省略内部 H 输入改变已授权对象。

全局时间键只存一次。每个 H 区间编码为 `(lower_integer,width_integer,exponent)`，准确表示 `[L/2^E,(L+W)/2^E]`；宽度非负，不用浮点、压缩精度或网格舍入。每条记录还存固定 id、规范有理投影 SHA、原 gzip SHA 与原 raw ancestry SHA。092 的 ancestry 取原 `original_raw_sha256`；095/103 保留原 packet gzip 身份。源字段映射及规范投影的完整规则见 `certificate.json`。

10 个事先固定的代表 packet 的字节布局估计约 4.0–4.23KB 压缩/节点，正式记录再含来源/投影头。故估计 15 个二进制 payload 合计约 0.42–0.44GB，完整必要 packet 约 **0.45–0.48GB**；这是有限样本估计，尚非导出实测。每片硬门为 90000000 B，严格小于 100MB，超门即停，不自动改片或扩输入。

当前固定 node index 为 **21182444 B**，含准确几何/offset 身份、L/U、m112/self112、当前 mean/self、原 point、源 packet SHA/路径与分片 id。绝对旧路径只供本机这次锁定导出读取，便携 reader 不打开这些路径。端点载荷正式导出后不再依赖外置 H cache。

## v2 并行导出与 transport 检查

原 `export.py`、`codec.py`、`export_plan.json` 原字节保留。原 plan SHA 为 `790353f4f15f9a546bfb1282465e111b561342283e932ac8ee15e47fddc83c4b`；它是串行版本，不实际运行。

新 `export_v2.py` 固定 15 个显式子进程，每进程只负责一个 sid。每至多 256 条记录写一个临时 gzip member，关闭/fsync/rename 后保存 chunk SHA 与写前投影；完成 chunks 可持久保留。每个实际新片完成后，该 worker 必须实际解码整片，逐条把读回规范投影与原写前投影比较，同时检查 id、record 数、完整 EOF、片前后字节 SHA，以及每个原 packet 写前/片完成后的 SHA。仅是精确 transport，不运行 `verify_packet.py` 的数学入口。

master 只登记真 master/15 worker PID、调度这 15 个固定进程、复制必要 metadata 和汇总。任一 worker 故障即停止剩余进程并有界 terminate/kill/join，保留完成 chunk/片，不自动重启。程序无 wall cap。main 的全部 metadata/source SHA 在前后核，H 原 packet 由各 worker 逐个核写前及片后 SHA。正式输出采用明确文件白名单；logs、故障临时文件和重复 staging chunk 不进证明归档。

只做过一个新布局测试：PID 92810，9 个已登记代表 packet（排除不用的 legacy 根）、657 区间，在内存中 encode/decode，全部规范化有理端点准确相等；正式片导出数为 0。该测试不重新证明旧 H 区间真值。

## 便携复核、复用与尚存边界

`verify_packet.py` 仅准备，未运行。decode-only 检查精确 transport；读者另行调用数学模式时可用原 111 权/112 J recipe 和保存 R 目录核 factor，按准确 band 差核 range/mean，按原 112 的整数四次根上包重建末端弦及树 min。这里不重跑那些旧 power。完整 I 的短时增量用原 139 正积分证书，原 `139/check.py` 可仅用保存根的四次幂检查 128 个面板/257 个根外端，不重求 H/R。新的 UH 与原合格 J139 的相加留给后续明确授权或读者，不在准备阶段求一个新 I。

库复用已按实际字节核：已有公共原件可复用 **10393677 B**，其中 139 原 packet/报告/checker 为 **209720 B**。`repository_reuse.json` 列逐文件 SHA 和从 `cube-dim-5-bounds/upper_integral167/packet/` 出发的仓库 relative 引用。公共归档可省去这些重复 metadata，使用 `stage_packet.py` 从同仓库原库进行 SHA 校验后的 hardlink/copy 到读者新目录；不做数学。独立 standalone 导出则只包含运行所必需的约 31.86MB metadata，绝不复制 1.54GB 原 raw/cache、DD 权缓存或无关日志。

必须明确证书边界：此方案把所有既有 H 球的**数值载荷**带入自包含 packet，并让 range/mean/chord/tree 算术可复核，解决了“只有 scratch 路径及自报 hash”的缺口。但原 H 球本身仍是已接纳的共享数值依赖，投影 SHA 证明身份，不证明它包含真实 H。若要求完全不信任任何旧 H 区间，还需要独立 H 核复算或完整 log/DD 数值见证；目前没有把它们伪装成已具备。本轮和计划导出均不做这种新求值。源解析定理也须阅读原证明，不能仅导入 PASS 标签。

当前待主审的唯一实际导出入口：

```text
python -B export_v2.py run --plan-sha 2893e7dd3df59bfe9fc9ae91e7610abfdeef6a0a69022e04d8ea63ede3e6c43a --approval integral167-parallel15-transport-export-approved --output NEW_EXTERNAL_DIRECTORY
```

`plan_v2.json` SHA 为 `2893e7dd3df59bfe9fc9ae91e7610abfdeef6a0a69022e04d8ea63ede3e6c43a`，锁 480 source 与 8 code。先交主审，未运行上述入口。未改 160/147/132 等冻结原件，未使用 MATLAB/Astra/下界，不提交或推送。
