---
name: math-g04-rubio-de-francia-square-function-scope-check
description: "用于树上 Rubio de Francia 平方函数界的数学研究：核对输入与假设，组织方法和证据，说明中心立方体极大算子的迁移条件；不凭名称补造内部接口。"
---

> **当前状态（2026-10-07）**：部分可用。经典频率区间平方函数的条件与 p=2 正交投影有限实例已核；加权/稀疏端点保留来源各自的空间和参数限制。没有从频率区间结论迁移到任意空间树或中心立方体最大函数。
> 外部文献证据随完整证据包提供；包根须同时含 `workspace_revised/`、`audit_current/` 与 `evidence/`。`EVIDENCE_ROOT` 和 `MATH64_SOURCE_PACKAGE_ROOT` 均指包根。单独安装的 Skill 不含共享 PDF/claims/cases；按[本地访问说明](references/portable-audit-access.md)从 registry root `math64-20261007` 解析包根以定位绝对 resolver，或显式提供包根。用 `--skill-dir` 与来源 ID 核验来源 SHA、用 `--artifact-path` 定位 claims/cases 并核验当前 SHA。registry、root 或 resolver 缺失时报告配置缺口，不猜路径/内容。resolver 的 PASS 只验证路径和当前字节，不代表全文阅读、命题证明或数学验收。


# G04 · 树上 Rubio de Francia 平方函数界


## 输入与产出

频率区间或树索引族、投影 P_I、输入 f、指数 p，以及不交/嵌套等精确假设。

产出应包含适用性判断、使用的假设、数学步骤、证书或界的方向、常数依赖，以及尚未解决的缺口。

## 执行步骤

1. 固定Fourier投影P_I及S(f)=(Σ_I|P_If|²)^{1/2}，列全频率区间族的两两不交或重叠重数。
2. 经典实线Rubio de Francia适用强型2≤p<∞的两两不交区间；p=2由Parseval常数1，p<2不能对一般任意不交区间族沿用同一uniform bound。
3. P-03df4e0c85ffb707的粗糙连续Rubio weighted p=2为weak端点，Walsh群或偶径向递减A₁权有各自strong结论；稀疏族依赖f与频率族。精确条件与[来源定位](references/frozen-source-locators.md)另列，不把加权端点或weak结果当经典全范围。
4. 树版本先定义空间树或频率树，证正交、有限重叠或Carleson packing；名称不提供这些前提。
5. 有限频率分割Parseval核常数，重复区间使平方能量按重数增加；任意重叠次数不能当统一常数。
6. 到空间中心立方体的迁移仍需频率分块与全尺度上确界连接；不交频率投影本身不控制该几何极大值。

## 证据与失败处理

Rubio de Francia 名称不能替代任意树平方函数的证明；区间不交条件、指数和端点都不可省略。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [适用案例](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。
