# F061 独立梯度参考准备：九配置通过，等待完整批 release

已接受新 physical time-pole PF 路线的准备证据。九项固定配置一次执行通过，
完整准备程序墙钟 7.253755 秒；1 次编译、0 失败、0 重试。尚未执行任何完整
两 case 参考，C016 的真实六梯度独立全群缺口目前仍未闭合。

## 已保存 pilot 的独立审查

作者 pilot manifest 精确为
`e41ac3b40e1b2010b1d0a5eaa8fd1e7cd463e06f8304e553b5f2c09d8dbb2ec3`。
其 32 个产物及 input manifest 的 26 个 source/runtime 哈希均吻合。二次原
native 的 argv、receipt、jet 都是相同原 body32 unit-weight 模型、相同
B_{1/250,1/2}、各自准确中心，六个 physical gradient 精确选择后五个 face
字段，没有窗口或数据混用。

新准备配置 1 逐项将实际 native 的完整 mantissa/exponent 转为 Fraction，
与保存 jet 的 F 和全部五个 face intervals 精确比较。四盒共 128 个角点的
位移、gradient dot、共享分母 clamp、Q32、exponent、向外端点次序、阈值标志、
各盒最小端点全部重建吻合，没有调用作者 kernel，也没有重新计算 exponential。
A057 的 exponential 算法已有 F057 独审，本次核对其实际输入/输出连接。

两个中心的 h=1/500 盒均 UNRESOLVED，h=1/1000 盒均条件 PASS；原始调用计数
是 2 native、0 fallback、4 boxes、128 corners。zero/repeated S3 terms 在
native 中保留；不以 Euler 单独证明六分量。原 release 有遗留字段
`template_only:true`，但同时有明确 main authorization 和 `production_authorized:true`；
源码使用后者及完整哈希门，且主线程确曾放行，故不是本次数学阻断项。

## 独立算法与新标量证据

`physical_pf_proof.md` 从每个 cone 的时间积分直接物理求导：

    Re integral_a^T i product_j(t-i b_j)^-1 sum_j d_jk/(t-i b_j) dt.

准确相等的 b 先聚类；每簇保存 multiplicity 及各 physical 方向的标签导数和。
原函数的 simple pole 系数纯虚，实部仅含两个独立 atan 之差；高阶 pole 用
复数倒数幂的端点差，b=0 无特殊删项。C++ 用局部 residue-series 卷积，完全
不调用 C016 DD/Hermite 或 normalized-jet 递推。global signs 显式遍历 64 个、
cone 权重 alpha 无倍增，最后除以 (2pi)^6。

配置 2 用符号计算核对 pole orders 1 到 7 的原函数导数。配置 3 到 7 为五组
全新 nodes：互异、两对重复含零、全零、六重 2/11、3/2/1 多簇。每组用 signed
Hilbert matrix 的六个独立方向；独立 SymPy apart 系数及 exact higher
primitive/Arb1024 参考，全部被 C++ Arb768 的 30 个 scalar field 区间包含。
这些 synthetic fields 不是实际 body32 gradient 证据。

## 真实 partial benchmark 与计数

每个真实 case 只取准确 index=7559*k（k=0,...,127），tuple 顺序为 lexicographic
permutation、64 sign masks、21 cones。两次共 **256 个真实 tuple**，每次共享
遍历同时计算六个 physical fields，共 **1536 个 tuple-field contributions**。
不能把这些工作报为“0 数学/kernel 工作”，也不能称为 12 个完整字段。

| case | partial tuples | repeated | zero | native kernel seconds | 剩余工作的线性估计 |
| --- | ---: | ---: | ---: | ---: | ---: |
| A038-low3 | 128 | 0 | 0 | 0.00894796 | 67.6376 s |
| S3 | 128 | 52 | 35 | 0.00532050 | 40.2177 s |

估计不是严格 runtime bound。S3 sample 实际经过重复和零 pole 路径；两者允许
重叠。每个 partial 保存全部六个经 (2pi)^6 归一后的区间，未来必须复用而非
重算这些 tuple。原样运行记录和 owned PID 均在 preparation receipt。

唯一编译成功，有一条 `-Wempty-body` 提示：line 44 的空 for 循环本来就是
按准确 key 查找已有簇的搜索，随后的 if 区分新簇/旧簇。这不是错误，未修改
已编译源文件或追加编译。运行 scalar/benchmark 的 stderr 均空、exit 0。

## 固定完整批，尚未启动

主线程已选定每 case 一个 remaining runner，全部程序硬限 180 秒。两个固定
准确中心与窗口见 `full_reference_plan.json`。每次新计算 **967552 tuples**，
跳过保存的 128 个；以准确有理下/上端点分别相加。于是每个完整字段恰含
967680 tuples，每 case 6 字段，两 case 总共 12 个完整物理梯度字段。

新 remaining 共享遍历合计 1935104 tuples / 11610624 tuple-field contributions；
加先前 1536，完整计划总数 11612160。scalar fixture 工作另计。不存在重复
tuple、旧 048 reference 重跑、旧 full-normalization 基线重算或作者 native 调用。

验收要求所有六个完整参考区间分别被对应 C016 Arb512 gradient 区间包含；
同时检查完整 repeated/zero counts 与保存 native 一致。S3 的每个 zero 或
repeated tuple 均参与聚合。未完成、超时或任一比较失败即停止，只 kill+wait
自己的子进程、保留输出与收据，不自动 retry 或改用分块，不宣称完整 gradient。

`run_reference_once.py` 仅完成 AST 语法读取，未执行。它要求主线程最终 release
绑定 preparation seal；固定顺序 low-primary 后 S3-primary，case 目录 exclusive
防止重复，且 S3 只允许在 low 全六分量比较成功后启动。完整参考仍为待放行，
没有 global g/coverage/C6 新结论。

## 源锁

- PF source: `500cb027288dbe1b68ec8530c0234cb9cd3b5c922011904c404eefd2800b7e92`
- Compiled binary: `9ffd29c7a2b31e4f1f5c1985bf9dd572f7fc18aaff8f3809d4eaf1bc14674f5c`
- Full runner: `f6a109ac421fee1e54f5c5d575562ec429e6a491c6b242744e6cec8221a815b2`
- low partial: `1a5a286af10069d20e9cfc9a87cf62e5ece3dabae2b10b86e23c8601728378a6`
- S3 partial: `ff140f22d7c81d8fff63440cd7676f4bdbaec150631d4fa2168d4e1d72c390b1`

完整精简白名单和哈希在 preparation_manifest.json。仅写本目录，其他 source、
F057/F048 只读，无 MATLAB、agents、stage、commit 或 push。没有活跃计算进程。
