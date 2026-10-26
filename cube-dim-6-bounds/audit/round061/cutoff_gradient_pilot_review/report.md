# F061：两个真实中心的完整六梯度独立交叉核验通过

结论：**PASS_TWO_CASES_ALL_TWELVE_PHYSICAL_GRADIENT_FIELDS**。
对固定 A038-low3 与 S3=(1,1,1,1,3/5,0)，原 body32 unit-weight 模型、相同
band B_{1/250,1/2} 的全部六个 physical gradient，独立 Arb768 参考区间均被
对应 C016 Arb512 区间完整包含。这两个实际 pilot case 的全梯度独立交叉核验
缺口已闭合。该结论不扩张为任意点/窗口的逐分量计算认证，也没有 global g、
全域 coverage 或新的 C6 上界。

## 独立算法与实际完成范围

参考直接对每锥时间有理核的物理坐标求导，并用 confluent physical-pole
partial fractions 积分。完整公式与标签重数处理见 `physical_pf_proof.md`。
不调用 C016 DD/Hermite 实现，不以 Euler 验证代替逐分量参考，也不重算旧
F048 reference 或旧 full-normalization 基线。

simple pole 的系数纯虚，实部原函数为两个独立 atan 的差；高阶 pole 是
复数倒数幂端点差，准确零节点也直接保留。每个精确重根簇保存标签导数之和，
保留重数。全群显式遍历 720 permutations、64 sign masks 和 21 cones，
使用 alpha=|det R|/product q_j，最后统一除以 (2pi)^6；不复用原 native 的
32-sign 权2实现。

准备批的五组新 scalar、每组六方向，共30字段，已由独立 SymPy apart/
exact higher primitive/Arb1024 对照通过。完整计算源和二进制在放行前已冻结，
之后未改源码、重新编译、调精度、换输入或增加 case。

## 每个 tuple 恰计算一次

准备阶段每 case 保存 128 个预登记 tuple 的六分量 normalized partial sum。
正式放行后每 case 只计算补集 967552 个 tuple，再将两份区间的准确 rational
下/上端点分别相加。

| case | 已保存 tuple | 新 remaining tuple | 每完整字段 tuple | 完整 physical 字段 | 结果 |
| --- | ---: | ---: | ---: | ---: | --- |
| A038-low3 | 128 | 967552 | 967680 | 6 | 全部被 matching native 包含 |
| S3 | 128 | 967552 | 967680 | 6 | 全部被 matching native 包含 |

两 case 合计 **1935360 个共享 tuple evaluations**。每个共享遍历同时计算六个
physical fields，共 **11612160 个 tuple-field contributions**，相当于12个
完整字段。不能把共享遍历误报成12次独立 full traversal。准确计数拆分为：

- 已保存 benchmark：256 tuple、1536 tuple-field contributions。
- 新 remaining：1935104 tuple、11610624 tuple-field contributions。
- 真实参考 native 进程总数4：两个 preparation benchmark，加两个 remaining；
  另有一个 five-fixture scalar 进程，与真实参考工作分别计数。

没有重复计算 benchmark tuple。raw remaining 输出仍写 `full_reference:false`，
这是因为其自身仅含补集；完整性和六分量通过由 `comparison.json` 中准确聚合的
967680-tuple reference 建立，不能把单个补集输出当作完整 reference。

S3 的独立全群合计 **405600 repeated-base tuples**、**191808 zero-base tuples**，
分别与原 native 的完整计数吻合；其中 benchmark 为52/35，remaining为
405548/191773。这两类可以重叠，全部参与积分及聚合，没有删项。low-primary
两类均为0。

每个 `full_reference/<case>/comparison.json` 保存六个完整参考区间和六个对应
native 区间的全部准确端点，以及逐分量包含判定。未从 display 小数或短 radius
重建端点。

## 原四盒同窗审计保持有效

准备配置1已将两份真实 native value/gradient 完整 dyadic 端点精确映射到 jet，
并复核四盒共128角点的位移、梯度下 dot、共享 F 分母 clamp、Q32、exponent、
端点次序及 exact target 判断。F、五个 face derivatives、center、model 和
window 同源；整盒接口是已审 A057 Q=I5/32。

| case | h=1/500 | h=1/1000 |
| --- | --- | --- |
| A038-low3 | UNRESOLVED，保存下端点显示约1.5278877211203848e-5 | PASS，下端点显示约1.5677763798311330e-5 |
| S3 | UNRESOLVED，保存下端点显示约1.5236195399366246e-5 | PASS，下端点显示约1.5632187853588424e-5 |

显示小数仅供阅读；判定使用保存的准确有理下端点与1/65000比较。本轮没有重跑
这些 exponential 或作者 kernel，亦未搜索更小半宽。完整梯度参考补上了这两个
matching jet 的独立梯度实算证据；它不会使两盒自动成为全域覆盖，也不把较宽
盒的 UNRESOLVED 改写为真实反例。

## 实际运行、release 与封存

主实际 release SHA256：
`797b894c7ed3ca14b25bc07c567b5af244ec13538ae1e04a84633803b86e73e8`。
它绑定准备 seal
`0d2ab3026f3fa735bd2656e804c76fe8fba60369496a1bee78eaf8a88649634a`。

| case | 原 tool session | owned native PID | exit | native seconds | runner wall seconds |
| --- | ---: | ---: | ---: | ---: | ---: |
| low-primary | 6619 | 66082 | 0 | 24.9248 | 24.971672 |
| S3-primary | 77263 | 66184 | 0 | 20.1137 | 20.154127 |

两个原 session 均经同一 handle 观察到 exit0，runner 已 wait 自己的 native。
两程序均在180秒硬限内，无超时、retry、chunk fallback 或残留工作。S3 仅在
low-primary 全六字段通过后启动。`full_tool_execution_receipts.json` 和各 case
的 receipt 保留真实 command、session、PID、哈希和结果，不用事后重跑补收据。

PF source SHA256：
`500cb027288dbe1b68ec8530c0234cb9cd3b5c922011904c404eefd2800b7e92`。
runner SHA256：
`f6a109ac421fee1e54f5c5d575562ec429e6a491c6b242744e6cec8221a815b2`。
所有准备源码、runtime、作者 pilot 及输入哈希在封存时仍匹配。原作者 pilot
的真实调用计数依旧为2 native、0 fallback；本审计新增作者 native 调用为0。

一次编译发生在准备期，成功并仅有已审的 intentional empty-loop warning；
正式完整批没有编译。MATLAB、old baseline/reference rerun、agents、stage、
commit、push 全为0。主线程保留理论与全局证明整合责任；本有限批完成并停止，
无进一步 kernel 或新窗口授权。
