# upper_cover048：冻结046余队列的8进程有限批

## 实际结果与范围

固定 \(\alpha=1\)、\(\epsilon=1/64\)、\(g_0=10^{-4}\)，直接 import 已冻结的046
`FirstJetKernel` / `choose_anchor` / `certificate`，没有重新求043资格点或046的128个旧锚点。
本轮唯一生产批最多4096新锚点或120秒，实际先达到4096新锚点上限：

| 项目 | 实际结果 |
|---|---:|
| 新锚点 starts / 完整完成 | 4096 / 4096 |
| 新通过闭盒 / 新 split | 4096 / 0 |
| 加上046后累计通过叶 | 4224 |
| 完整剩余 queue | 48136 |
| 批尾 inflight / partial resume jobs | 0 / 0 |
| 新 kernel 次数 / 精度 | 4096 / 全256-bit |
| 资格交叉重跑 / 旧点重算 / cache 命中 | 0 / 0 / 0 |
| 新锚点采用首个 jitter | 4096 |
| H 目标反例 / 数值或 worker 失败 | 无 / 无 |
| 全批墙钟 | 40.929237208秒 |
| 完整 state compaction / journal事件 | 66 / 12289 |

新通过叶的最差准确有理 log 余量为
\[
 \frac{975981153904265447934476086315595512766167433344965245718851472835292823233741}
 {3705346855594118253554271520278013051304639509300498049262642688253220148477952}
 \approx0.2633980547410253>0.
\]
最差叶编号3860，根索引 \((15,15,15,0)\)。新中心 \(H\) 包围下端的最小值约为
\(0.000358034798332182\)。所用门槛检验不是小数比较，而是16角点准确有理 log 余量严格正。

实际吞吐为约100.08个新完整锚点/秒（包含调度、校验、I/O和池初始化/关闭的全批墙钟）。
各 worker 报告的完整 kernel 墙钟之和约279.3826秒；它是各调用墙钟之和，不是 CPU 计时。
主进程报告的峰值 RSS 为1098547200 bytes，单个 worker 报告峰值最大42434560 bytes。
没有 MATLAB，也没有操作其它任务的进程。

选定节点只是046原剩余 queue 的前4096个：根编号128--4223，首索引 \((6,1,1,0)\)，
末索引 \((16,11,10,6)\)，下一待处理索引为 \((16,11,10,7)\)。
所有新根的最大尾索引仍不超过16，没有到达018中心附近的困难区域。
本轮只证明 **partial certified cover + 完整剩余 queue**；没有全域 \(g\) 或 \(C_5\) 改进声明，
不能用本批通过率/吞吐估算最终细分量、困难区域精度或总完成时间。

## 1. 冻结父状态及新唯一活句柄

输入只读父状态为
`/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover046/batch001_state.json`，
SHA-256
`b23f4ee65bc7ce1d8bdee7f0050b4e492b453241a6594cce8e79f96e7bcfbb2f`。
046的 `state.json` 与该批快照均保持原字节，未调用046 `--resume`。
原128个通过叶、cache和资格记录直接复用；本轮结束还逐项检查128个旧 cache 记录完全未变。

登记前准确取原剩余 queue 的前4096节点，其 canonical JSON hash 为
`389af47a843955abc366b8b36eadb5d7134b1eabeef68710cab5e45fa2691e46`。
父 fingerprint
`5ac95d23a36f4bfea5fabd7a2d834e9045e774a78d2cc318d7124e46ad82ca9e`
按046源/登记/上游/043来源重建并核对，不仅比较一个存储字符串。
046四个源文件、017/018三输入、043四原件也逐项核 SHA 后才启动池。

新048 fingerprint 为
`9e73e52f6f4acda7f378c542096cd921cdcf96961a48489bb99d2fa0c8cd7ef3`，
同时锁定新源/登记、046四源、父状态字节、上游、043及选定队列。
新唯一活覆盖句柄为
`/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover048/state.json`。
它显式记录 `046 frozen parent -> only active048 handle`、父 hash及父 fingerprint。

本轮五个新文件为登记、生产驱动、结果摘要、只读验证摘要、本报告；没有修改任何旧锁或他人目录，
未提交。临时队列、journal、完整 cache 及池产物全部在仓库外持久 scratch。

## 2. 并行生产与唯一 writer

主进程 PID48042；8个 spawn worker 的 PID 为48051、48052、48053、48054、48055、48056、48057、48058。
每个 worker 初始化时独立核被冻结 API 和017表 hash，读表一次，然后只评其独占节点的 jet。
pool ready后先报告以上PID及 fingerprint，再派发独立节点。

主进程最多保持8个在算 job，同节点不派给两个 worker；锚点坐标 cache 和在算坐标也检查重复。
保留046规则：中点不合格时固定三个微小有理 jitter 中选首个280投影非零且互异者，
256-bit主精度，仅数值不充分或中心目标未分离才384/512-bit重试。
实际4096中点均退化，首个预登记 jitter 均合格，未求其它候选的密度、未升位。
本轮没有改变 \(\epsilon\)、基、门槛或节点顺序。

完成顺序允许乱序。主进程立即把收到的完整结果写入 journal，再严格按预登记 ordinal 合并
pass/split/diagnostic；worker不能改 queue、cover树、源文件或总账本。
整盒通过仍使用046/043的下 log 抛物面：\(Q=I_4/4-J_4/36\)，准确梯度有向端点，
相对于真实 jitter 锚点计算原闭盒的16个角点，最小准确有理余量严格正才通过。
整个盒内出现零/重复投影不影响正混合定理；普通 DD 只在已验证合格的中心执行。

若整盒几何证书失败，最长边/最低坐标 tie的二分子盒追加到queue，留下一批，本批不再派发子盒。
若收到 \(H_U<g_0\)、数值/几何诊断或 worker error，则停止新派发，已经派发者收尾且保留真实结果。
固定 \(\epsilon\) 的 H 目标失败不能冒充完整 Green 或 \(C_5\) 失败。
本批没有上述失败或split；相关恢复分支不能称为已实际压力测试。

批尾明确执行 `pool.close()`、`pool.join()`。保留8个进程对象检查：8个 exitcode均为0，
8个 `is_alive` 均为false；这些证据随结果和最终 journal 保存。没有遗留该池 worker。

## 3. 低重复I/O的账本及完整覆盖分割

使用单主进程 hash-chain append journal。每个派发记录有节点、ordinal、实际有理锚点、资格测试；
每个收到结果记录有完整5元球的准确 dyadic端点、精度和尝试记录。
每个合并记录包括决定及完整角点证书。完成的 `received` 小记录作 fsync，
其前的 dispatch/merge也随此前缀持久化；完整状态每64合并作一次原子compaction和fsync。
实际初始1次、64个周期compaction、最终1次，共66次，没有逐叶重写18MB状态。

compaction时把所有 root沿split重建叶集合，检查严格分割
\[
 \mathrm{current\ leaves}=\mathrm{passed}\ \dot\cup\ \mathrm{queue}\ \dot\cup\ \mathrm{inflight}.
\]
已经乱序收到、尚待有序合并者仍算 inflight，不能从分割中失踪。
局部派发检查 queue头准确等于登记节点，合并检查ordinal及父叶，split保留两个闭子盒。
本批实际没有split，最终准确有
\[
 52360=4224+48136+0.
\]
额外只读检查确认：通过叶准确为根编号0--4223，剩余 queue准确为4224--52359，无重复、无遗漏。
完整根覆盖闭Weyl域及置换传递的证明沿用046；本轮并未把剩余48136盒当作通过。

worker使用冻结 API 的128锥期限检查。若到120秒，保存已完成的完整锥前缀及对应位数/尝试；
部分和只作 `partial` 状态，不能作为完整H使用。停止派发后先等已有job返回，
未解节点按原ordinal重新放到queue前端，前缀放入 `resume_jobs`，不丢盒。
本批先达到4096上限，无partial/inflight遗留；没有实际触发120秒partial路径。

## 4. 已实际完成的无密度恢复与证书复核

初始048完整状态和全 journal均保留；`--recover-only` 只执行hash-chain验证、事件重放和叶分割检查，
不创建池、不执行密度计算。实际从初始状态重放全部12289事件，得到恢复状态字节 SHA
与生产最终状态 **完全一致**：
`edd97dc22763426c36427e152f263c5090a2114bfa2a4695d0a9b347f1992ed7`。
这验证的是本批实际完整事件流的可恢复性，不是对未触发的中途进程崩溃/partial路径作压力测试。

另外对4096个新叶的全部65536角点做了只读准确分数复核，耗时约1.5915秒，无新密度求值：

* 每个采用锚点属于原闭盒；128旧cache及旧资格记录逐项未变；
* 384-bit标量log确认存储 \(\Delta_L\le\log(H_L/g_0)\) 的方向；
* 准确端点除法确认存储log梯度区间包含五元球给出的各端点商范围；
* 用存储有向梯度和准确Q独立重算每个角点余量，最小值与叶证书准确相等并严格正。

复核摘要见 `validation.json`。本轮资格重跑数仍为0，新完整kernel数仍为4096。

## 5. 持久文件、hash与复现

下列完整数据在仓库外目录
`/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover048/`：

| 对象 | 字节 / SHA-256 |
|---|---|
| `state.json`及同字节`batch001_state.json` | 18658424 / `edd97dc22763426c36427e152f263c5090a2114bfa2a4695d0a9b347f1992ed7` |
| `journal.jsonl` | 20928635 / `046925cfad5879ac7ffc99fc7ee596091630b780ab5d6fb9c2ca73848e913a0b` |
| 完整cover树 canonical JSON | `03411d1d93df485eefcbbd7f726c84f0075f179525ecb5db1edf80f31fe6e5bd` |
| 剩余完整queue canonical JSON | `ad891102cdc431938b75b1dbfe951010d0a2f38af2f5723d51d7eb6b29247865` |

journal最终链hash为
`255156dd4bd42014e39e511eb4fa7fc5420c1a5a21ddbf49ea0e309f19a27132`。
`initial_state.json`的准确hash记录在结果JSON，`recovered_state.json`与最终state同字节。
最终state已经包括完整树/queue/cache、父迁移记录、048源指纹、池关闭证据及批次计数，
不是只有一个摘要计数。它是主后续批准新有限批的唯一活句柄；本驱动不自动开启下一批。

生产命令：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B \
  cube-dim-5-bounds/upper_cover048/parallel_cover.py \
  --source046 /Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper/cube-dim-5-bounds/upper_cover046 \
  --upstream-root /Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-upper/cube-dim-5-bounds \
  --parent-state /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover046/batch001_state.json \
  --work-dir /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover048
```

实际无密度重放使用同一命令再加
`--recover-only --recovery-output /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover048/recovered_state.json`。
若要重做生产，应提供全新的仓库外work/output目录；已有048状态会拒绝再次派发。
本批停止后没有第二个4096批，没有046续算，也没有自行提交。
