# upper_continue123：两有限批6237新H128的冻结结果

已按同一主审版本完成4096+2141两个原前沿批，8worker/nice5/no wallcap；没有派任何新生child。6237个真实新H128及24948个物理梯度分量，新增572通过、5665几何split，最终**98579通过、57549split、11330待证**。fault/obstruction/inflight均0；所有新H128下端>110µ，最小约119.151880165570µ。仍为partial g110，S记0，没有I或C5更新。

|批|master|worker PID|新H|pass/split|wall秒|
|---|---|---|---:|---:|---:|
|1|30663|30704–30711|4096|307/3789|86.130|
|2|33802|33890–33897|2141|265/1876|60.805|

两pool均正常close/join、exit0，OS已无。原98007通过叶及189086 typedcache不改；新cache按EPS128/point独立追加6237项。全部新点256bit、3504cones、280generic、1segment。生产仍用已审target110、物理同EPS梯度和**Q1024**；没有混入后来的新q128理论。未提交、MATLAB或新/续派Astra。

作者每批门实际成功PID33007/34968，分别11.912/11.030秒，零H复核6237个768bit Delta、99792角点、18713 journal events、cache/source及52360根闭Kraft。第二批只在第一批作者门全过后启动。第一次作者脚本在写最终receipt时因argparse变量被端点变量遮蔽而退出；保真摘要另存。仅修复未冻结的作者helper并重跑零H门，**不是生产故障，没有H重启或重复生产**。

## 最小归档与原始SHA

`packet/`保存两份原始point manifests、两份冻结输入及最终累计增量state的确定性gzip；每份压缩SHA、解压原字节SHA和字节数在packet manifest中。没有复制761MB seed，也不依赖活pool。第一批state通过最终state的精确子集+原header恢复，`restore_packet.py`已实际字节验证两个原始SHA：

- batch1 state `75aaf1ef83f62b2e72c2b524a18af9749ab1ff50f5c276a0610be2ba6834e4a9`。
- batch2 state `6c1d4c6e075053c427d797798b5ccfc20e62e46a46151511e01a330b38213301`。
- batch1 point manifest `40342fbe1d0bf596dc3a31aac741237ad926b17539c8f96a5988f1e3e331b656`。
- batch2 point manifest `faf4de812b216875bebbddc802843ca8f28a9523b0f5e93a108a563374ff04e2`。

归档的实际registration/pool/summary/author receipts与source refs绑定096/098既有packet和105/107/118来源哈希。096 packet为160bit向外量化的既有g105叶证据，不伪称它能恢复原761MB rawcache；本包自含本6237个**未量化原始新值/梯度/target110证书**及两个精确增量state，可单独复原这些新结果。全局基线资格继续依赖已冻结来源与主/独审。

```sh
python3 cube-dim-5-bounds/upper_continue123/restore_packet.py
```

该命令仅解压、record逐项比较和原始state字节SHA校验，0H/gradient/log/corner重算。`archive_whitelist.json`为本次精简白名单；外置完整journals/point checkpoints仍保持原锁供audit123协议审计，未机械复制所有重复cache或旧baseline。131只作为下一轮准备入口，不在本轮声称真实有限认证已完成。
