# 第247轮：实际动态分母与共同测试容量

一般 n^(1/2+o(1)) 一致弱型上界尚未闭合。本轮独审原动态测试矩：在原完整正交通dP上，R=T e_t/(qD)满足∫R dP≤2Λ(16/3)X，R≥H的交通可按该常数/H吸收。它同价包含R226静态测试门，实际新信息是保留动态归一化分母；D_Z在原terminal-good槽等于旧D，不重复登记。补域R<H仍没有一般根阶来源预算。

[研究包](research.zip)保留详细LaTeX文本、作者与独审、辅助路线缺口、8例三阶段本机数值和完整复现驱动。不编译PDF。有限clock动态门有新增交通，不认证FIRST/真实cube/FULL；40/70位对比不是区间证明。

## 复现

解压后，用带mpmath的Python运行：

```sh
python3 reproduce.py --verify-only
python3 reproduce.py --output-dir /absolute/new/path/outside-package
```

驱动核manifest后在包外重算完整8例JSON，逐字节比较，再核manifest；结果写receipt.json。当前机器依赖为Python3.12.14、mpmath1.3.0，具体命令见包内numeric_notes.md。24个精确D槽、75个完整后验中点检查、24组门对照通过；普通多精度最大归一化差4.000896423e-41。

仅新增本README及research.zip，更新项目总览。完整本地检查点与复现回执保存在仓库外；无原始PDF或重复历史归档。下一步仍需用全局极值几何支付原联合补域，不能用此分支代替一般闭合。

研究包：35,245字节，11个成员；SHA256 `fa0c89150cdcc23ab13a3d7c0f4953281137bacb24322e7577303d8180bbb0c2`。主代理完整包外复跑通过，整个JSON逐字节一致，前后10项manifest哈希一致。
