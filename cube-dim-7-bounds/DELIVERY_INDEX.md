# 一般高维覆盖：当前结果与逐轮交付索引

项目目录：`cube-dim-7-bounds/`。目标active；下界暂停、七维专项停止、原资料包/旧暂停项目只读。

当前最强已接受eventual完整成本是[R064](upper/round064/REPORT.md)：Cn≤Phi64<T−n/5，未知起点。已显式认证[R063](upper/round063/REPORT.md) log n≥2^33的n/10和[R056](upper/round056/REPORT.md) log n≥2^20的n/40。所有比较使用完整F和各自valid域；首项系数仍1、内部未发表。

当前[24页PDF](paper/main.pdf)和[TeX](paper/main.tex)是R060冻结稿，全文数学FINAL ACCEPT、主双编译和逐页图审，证明范围R056/R059。R061/R063/R064是后续报告，不冒称已写入该PDF；后续论文将统一修订。paper/BUILD.json与MANIFEST锁源/PDFhash、原运行身份和旧20/30页Gitblob恢复。

| 轮次 | 独立提交主题 | 实际提交 | 完整报告 | 必要变更文件数 |
|---|---|---|---|---|
| R041 | 一般高维论文定稿与终审 | `38786ecb77ff244047a13c314dc5e2e08433f421` | [报告](upper/round041/REPORT.md) | 25 |
| R042 | 二阶矩截断与显式高维上界 | `5abe7f03346219ed123a7252fed8b0cc8b42ca2b` | [报告](upper/round042/REPORT.md) | 27 |
| R043 | 准确矩与截断模板范围审计 | `3bd7f83b9572bc155f9aa02ed9979e9e29cc802d` | [报告](upper/round043/REPORT.md) | 12 |
| R044 | 联合支付与小块重分配证书 | `324af2638a6cb18c84b0a630f9faf6ee673fb03e` | [报告](upper/round044/REPORT.md) | 37 |
| R045 | 逆列加权与负线性高维修正 | `f7ce8805bfd5c2a5140612082c07aa358cfb4b0e` | [报告](upper/round045/REPORT.md) | 28 |
| R046 | 逆列行积分有限数值诊断 | `718972058247d636550d1683ce0eba8d64e266e5` | [报告](upper/round046/REPORT.md) | 20 |
| R047 | 完整固定面交叉矩与本机满核诊断 | `4dc93fb9da32b47d07ec2f9fc5fea7fd4f73ad12` | [报告](upper/round047/REPORT.md) | 23 |
| R048 | 逆列主定理论文与全文终审 | `5122f5eba011430a857692896772ddffd02fba5d` | [报告](upper/round048/REPORT.md) | 22 |
| R049 | 真实算子显式门槛与全尺度闭合 | `38801fbd17a9019e37aea8359f3d422ecdb3a8f6` | [报告](upper/round049/REPORT.md) | 12 |
| R051 | 灵活逆列预算族与准确实例激活 | `3bfd3249d30dfc7b91ca69b21434575480e1377f` | [报告](upper/round051/REPORT.md) | 11 |
| R052 | 准确预算搜索与独立保留证书 | `e81f8d0f8a2643dbab5410a6ea2c95c4d572916e` | [报告](upper/round052/REPORT.md) | 25 |
| R053 | 准确参数完整算子门槛与根式比较 | `08a4dbf005b589522f943cf6d129ac1361dbef44` | [报告](upper/round053/REPORT.md) | 13 |
| R054 | 最新准确上界二十页论文与全文终审 | `4474b8be0cc47c904b23a5327108e37b5e6aea79` | [报告](upper/round054/REPORT.md) | 24 |
| R055 | 简单准确构造与更早真实门槛 | `324a616a9cbf2323bfefcc905055f5ddf0350fca` | [报告](upper/round055/REPORT.md) | 13 |
| R056 | 计数校准与逐面方差严格上界改进 | `2a9b922514efbcb00cba9c4c7cfba39e0e1b31ae` | [报告](upper/round056/REPORT.md) | 41 |
| R057 | 捕获边界导数诊断与精度失败审计 | `464637579f8d71b4200a2876a2d6cd9f5d6d1941` | [报告](upper/round057/REPORT.md) | 25 |
| R058 | 宽来源带模型有限搜索与准确条件证书 | `174473ceb8c1f7287ccc422df841085945b1cc1b` | [报告](upper/round058/REPORT.md) | 24 |
| R059 | 宽来源带完整严格上界与一成线性修正 | `0937192e6b5cef98679a144176fe4be1fd23a3dd` | [报告](upper/round059/REPORT.md) | 30 |
| R060 | 精化覆盖与宽来源带二十四页论文终审 | `2f393a17f443f10411fa9ea435fd88fd5350c14a` | [报告](upper/round060/REPORT.md) | 24 |
| R061 | 同一真实历史不交分层的严格线性改进 | `c33eff40c5a91305db69232db81c36cbddb21112` | [报告](upper/round061/REPORT.md) | 20 |
| R062 | 分层直接和前沿本机计算与准确层证书 | `86473c652e19d4b7abdd83f938786f28e71a3e18` | [报告](upper/round062/REPORT.md) | 33 |
| R063 | 宽来源带一成修正的完整显式门槛 | `a4da5e5d5744e31a8536923043e3e00c8a4b8830` | [报告](upper/round063/REPORT.md) | 17 |
| R065 | 乘法相对核有限诊断与真实错误符号审计 | `87904a9de0961cba253397fe9021fccbc57769df` | [报告](upper/round065/REPORT.md) | 25 |
| R064 | 乘法匹配与径向校准 | 本轮commit见Git log及最终readback | [报告](upper/round064/REPORT.md) | 以本轮精确白名单为准 |

以上历史提交均已确认远端包含；本轮R064只有真实同步receipt成功后才报告已推送，避免把自身hash写入自身。每轮原件/来源SHA、最小code、原actualreceipt、独审及清理范围在REPORT/MANIFEST；过去轮次另见rounds与准确账本。

R062本机v2实际136.02秒、累计CPU1846.22秒、平均13.57核，采样峰14.588核；R065为19.7459秒、累计CPU247.9498秒、采样峰14.488核。R067新批次数据由自身正式审核/报告单独记，不以运行意图或worker数量冒称持续满核。数值候选/条件证书/计算失败/完整Cn严格认证分开；错误tau符号是真实假失败而非数学反例。

归档采用具体文件白名单，保留必要原UTF8收据压缩单份，不上传全中间表、缓存、截图或重复旧论文。未知来源同名2未跟踪文件保持原样。R066尚待full审，不把更好targetprofile模型自动升级为Cn。
