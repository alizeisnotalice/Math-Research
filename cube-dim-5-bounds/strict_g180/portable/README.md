# 180 精简归档与仓库内复核入口

本目录保存实际180审计的精选原件、原80源锁与原71输出hash身份，以及仅依赖仓库已有151/155数据的数学复核入口。原数值运行是PID31006及15worker，完整结果与独立性限制见 `evidence/report.md`。`evidence/manifest.json` 原样保留当时71件输出身份；其中15份journal不在本小归档内，其长度、整字节SHA与当时实核的末链见 `omitted_local_evidence.json`。原324MB输入包及268MBjournal留在原本机目录，均不重复入仓。

归档动作不重新计算EXP、q、除法、角点或forest。`archive180.py` 只导入参数字段模块，逐一比较原15输入包与从repo151/155恢复的107530源、3812支撑、111201叶，随后字节复制及回读。实际结果由 `equivalence_result.json`、15份equivalence-shard和 `portable_manifest.json` 记录。若字段不等，程序停下保留失败记录，不用数值推算补齐。

| 原输入中的资格字段 | 仓库恢复来源与纯字段处理 |
|---|---|
| canonical完整集合、同源point/epsilon/Delta110/四个gamma/oldq/规则 | 151 coefficients的完整core，原样复制 |
| H及四个物理导数区间、model、typed epsilon、point、bits、raw SHA | 151 jets完整记录，按jet_id关联，原样复制 |
| 参与源的old_H_interval与typed_cache_key | 同jet的values[0]；core epsilon及未舍入point字符串拼接 |
| 新q及其原事件hash绑定 | 155 packet/new_q完整行，原样复制；本次不算q |
| H64质量L/U/epsilon/point及真实support node | L的原资格常数字符串、同jet U、同core point、原proof alias对应151 support；不作除法或包含计算 |
| mass完整来源hash与布尔资格字段 | 151 support原源hash；证书规范JSON SHA；恢复原149 alias后作JSON SHA；无数值运算 |
| 原proof alias选择 | 保持151 provenance原行序，沿用155的首alias/首个curvature优先字段选择规则 |
| 151中省略的149 alias margin/target | 按同source_alias_index取151 leaf的claimed margin，及同core target字符串 |
| 闭叶node/真实box上下端/epsilon/canonical与旧余量 | 151 leaves原字段；与原155 source leaf逐字段完全相等 |
| 新margin110/margin115、worst bits、passed115及原事件hash | 155 packet/leaf_margins完整行，原样复制；本次不算角点 |
| 3812个质量支撑完整证书、原叶hash及其jet | 151 supports+jets原样关联、复制 |
| Delta115的意义 | 其输入视图为同source Delta110与同common_log_U；原相减身份继承实际31006。本次不执行减法，未来显式数学入口才重新核叶的shifted margin |

完整mass参数与来源hash均可由仓库字段恢复，无需新增数值补片。80个原source锁仍完整存在 `evidence/plan.json`，包括当时9个raw state的身份；portable入口不要求这些外置state，也不冒充重放了raw回连。109446个raw记录逐字节回连是实际31006的继承证据。原rigorous H/物理导数enclosure仍是共同依赖，尤其H64_inherited052的原值不会因小归档而变成新的独立计算。

未来如需重新执行数学复核，在仓库根目录运行：

```sh
python3 cube-dim-5-bounds/strict_g180/portable/portable180.py --verify-mathematics
```

它仅依赖Python标准库与仓库已有 `upper_global151/`、`upper_retarget155/packet/`，使用15进程；从精选原源码复用实际测试过的numeric_core/full_corners/box/roots/common_log及其有向整数/Fraction算法。新适配入口的数学执行本轮为0，不能将AST检查和字段等价说成新入口已完成数值测试。省略该显式选项只打印说明，不做数学。

未来入口验证的H128统一floor仅在F，随后G≥H128及完整G齐次性才扩到Q。P8/down192继承178的实际独立字段核验；本入口不重算I、P8、根或下界。原C5上端28936079/1000000及原研究未完成条件均不改变。
