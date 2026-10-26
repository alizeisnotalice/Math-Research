# 167 执行附录：精确 transport 已完成

冻结的准备报告、v1/v2 源码和计划原字节保留。本附录记录后来取得的唯一实际授权与运行，不把历史准备文字改写为已执行。

主授权文件 SHA `0ba083b1c0f5d782aad9a5be0307e2593cb9cec034482c2b59834f563dcd667e`，唯一 plan v2 SHA `2893e7dd3df59bfe9fc9ae91e7610abfdeef6a0a69022e04d8ea63ede3e6c43a`。新 external 输出为 `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_integral167/export001/`，没有重复启动或自动重启。

实际 master PID **99171**，15 个显式 worker PID **99172–99186**。master 保存的实际 wall 为 **17.183933332999004 秒**；15 个 worker 均经 master 的 wait/poll 取得真实 **exit 0** 并 join。完成后逐 PID 的 OS 查询确认 16 个进程均已不存在，所有 fault 文件数为 0。

运行完整转换了 **97971 个 packet、7151883 个 H 区间**。15 片 payload 总计 **405181695 B**，最小片 **27000016 B**，最大片 **27027971 B**，均远低于 90MB 硬门。每片完成后，其 worker 已实际整片 decode，逐条重建规范有理 projection 并与写前保存 projection 对比，核 id/record 数/EOF/片 SHA，且核全部原 H packet 写前和片后的 SHA。15 片全部通过。

最终另作字节 readback：全部 94 个输出白名单文件 SHA、480 个 metadata/source 锁及 8 个代码锁均一致；没有再次解码全 payload，也没有运行数学 checker。实际 H/R/factor/mean/range/chord/旧 power 求值数均为 **0**。此工作没有重新证明旧 H 球真值或源解析定理。

资源记录有明确缺项：首个 `ps` 采样请求到达时 17 秒作业已全部退出，故**没有活 CPU/RSS 样本**；原程序也未保存历史 CPU/RSS。短 launcher 未保留 master 的 OS wait status，所以该字段为 null，不能把完成 receipt 冒充实际捕获的 master exit code。15 个子进程的 exit 0 是 master 实际取得并写入的结果。没有补造资源值，也不为补采样重跑成功导出。

源资格在本次运行前已由主更新：`round018/primary_integral156_acceptance.json` 实际主 PID **91014** 接纳了 I156；该文件中的 UH 与这份 packet 的既定 candidate UH 完全相同。原准备状态保留历史；**新的数学资格由 161/主证明承担，不由这次 export 证明或改变**。本附录未重新合成 I 或 C5。

## 必要归档清单

- `archive_standalone.json`：95 个必要文件（含原 completion），共 **447783353 B**，可在一个目录运行 reader。
- `archive_repository.json`：只 copy **57 个文件、437389676 B**；其余 **38 个文件、10393677 B** 使用已实际核同 SHA 的仓库 relative 库引用，其中 139 原件仅 **209720 B**。`stage_packet.py` 可在读者新目录从这些 relative 原件作 SHA 验证后的 hardlink/copy，形成 standalone reader 目录。
- 原 21182444 B node index 只应在 `packet/node_index.jsonl.gz` 归档一份；若同时导入准备源码，需要根部输入路径时可采用相对 symlink `node_index.jsonl.gz -> packet/node_index.jsonl.gz`，不再复制第二份。

不归档 `worker*/chunks/**`、worker 私有重复 projection map、空日志或临时文件。390 个 staging chunk 仅在 external scratch 保留作故障/生产审计，不复制到 Git。没有把完整 1.54GB 原 cache 搬入证书，也未提交或推送。

原 output completion SHA：`555e4cfb1b66d0f6270988e3fed6d9548c1b84c1630917d64101d3ac58728770`。

执行 receipt SHA：`cc437aa30f3b845ecdef51fb90aa85820a18593befc60416dd7ce9c55d4eac6d`。

仓库 relative 复用清单 SHA：`f2860e6e66b5a7c14635295dc6a7fcf5405d9ada70afff16a21972e5a77675ad`。

standalone 清单 SHA：`bf2941be78cbfbedfe4f56bd8b68f70cc8cccc29aec6e0c47a6a6bf2c3ed89c4`。

数据、receipt 与清单可供主归档；本代理不提交、不推送，旧 160/147 等冻结原件未改。
