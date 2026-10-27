# 五维中心方体极大算子的准确上下界

当前为第022轮按用户PDF逐段复核的写作修订稿，12页。论文 `cube5_upper_bound_zh.pdf` 保留原文件名，源为 `main.tex`。按照用户六页标准组织为引言、预备知识、主要证明、应用；原运行与复审过程转留工作记录。研究保持暂停，没有新界或数值搜索。第020轮原稿可由Git提交 e5b1c4187038104ef17c4ff4e5ee73d14c68a20b 恢复。

采用本机 TeX Live 2026 XeLaTeX/latexmk，ctex/Fandol。最终源码及PDF指纹、编译和视觉检查见 `final_qa.json`，本轮修改范围及一致性检查见 `../round021/` 和 `../round022/`。原 `compile020*` 文件仅对应旧稿，不作为新稿运行记录。

编译：

```sh
latexmk -norc -xelatex -interaction=nonstopmode -halt-on-error -outdir=build main.tex
```

正式结论仍为 C_cert ≤ C5 ≤ U175 < 14460339/500000，C_cert > 19943/10000。准确下界25系数来自 `../lower_round003/a2_equals_a3_results.json`，上界完整分数、根区间及准确比较在 `../round020/primary_upper175_acceptance.json`。新版将论文中的U175简记为U_*；没有改证书。

## 数学补充材料

- `../primary/finite_source_transfer.md`：全空间有限源、同一接收集、一般输入及全部尺度证明。
- `../final_bounds020/lower/witnesses.json`：42见证盒的长度、原子约束及端点；同目录为1122不交盒重建和准确25系数。
- `../round020/finite_collar_result.json`：有限光滑见证的严格余量。
- `../strict_g180/portable/`：全域闭覆盖验证及仓库数据参数视图。
- `../second_audit178/`：完整积分验证器、时间尾及准确结果。
- `../upper_secant175/ARCHIVE_EXECUTION.md` 和 `../audit_secant175/`：新弦的原始根证书、独立整数复核和未来显式复核命令。
- `../upper_integral167/packet/`：完整原H端点输入。
- `../upper_combined_archive179/`：组合积分数据的无损恢复、来源映射与hash。
- `../round020/paper_constants.json`：短整数常数；其程序仅检查最后标量不等式。
- `../results/bounds_ledger.json`：完整准确账本及公开、原内部、接续比较。

上述数值验证均为此前真实运行。本轮只运行编辑一致性检查、最后标量检查和编译，不重复它们。原H/物理导数区间仍为共享输入，不同验证路径亦共享Python整数算术，未声称第二套H引擎或形式化验证。173仍是未经独立数值审计的候选，不进入论文。

第020轮的审读、身份字段勘误、执行回执和排版等价记录保留在相应归档。第021轮移出正文的流程信息可由原提交恢复，不删除历史证据。
