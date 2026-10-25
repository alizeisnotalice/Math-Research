# 012 — C 的最小可携复核入口

`012_portable_C_verify.py` 是为 C 归档准备的单文件复核器。它只直接导入
同目录第 002 轮 audit-owned 的 \(\mathbb Q(u,v)\) scalar/root class；C 的
1024 长度向量、42 盒筛选与 slab-bitmask union 均在本文件中实现。它不导入
生产的 `field`、box 或 union 函数，也不依赖第 010/011 的调用链。

输入按项目根相对路径读取：

* `lower_round003/a2_equals_a3_results.json`
* `lower_round003/lock_c_a2_equals_a3.json`
* `lower_round003/lock_b_1_32.json`
* `lower/rectangular_certificate.json`

默认项目根是本脚本上一层的 `cube-dim-5-bounds/`。本次实际复跑时，为访问
下界 continuation checkout，显式传入了该 checkout 的项目根：

```text
venv/bin/python cube-dim-5-bounds/audit/012_portable_C_verify.py \
  --project-root /.../Math-Research-cube5-lower-continue/cube-dim-5-bounds
```

实际用时 3.11 秒并通过：C 系数匹配、42 盒、审计区间包含在 C lock 区间内，
且严格超过 B、001 及 \(19943/10000\)。输出写为新文件
`012_portable_C_verify_results.json`，没有改写 011 的运行证据。

独立性限制未变：它共享第 002 轮审计标量/root 实现，故不是第三个独立代数
内核；其几何和联合体积代码独立于生产器。
