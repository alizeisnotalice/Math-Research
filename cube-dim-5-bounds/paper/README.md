# 五维中心方体极大算子的准确上下界（暂停后的审计归档）

当前论文 `cube5_upper_bound_zh.pdf` 为12页完整上下界稿；保留历史PDF文件名以兼容原链接。源为 `main.tex`，实际采用本机TeX Live2026的XeLaTeX/latexmk，ctex/Fandol。源码与PDF的准确SHA及视觉检查在 `final_qa.json`，实际编译证据在 `compile020_final_receipt.json`。

从本目录编译：

```sh
latexmk -norc -xelatex -interaction=nonstopmode -halt-on-error -outdir=build main.tex
```

结论：C_cert ≤ C5 ≤ U175 < 14460339/500000=28.920678，C_cert>19943/10000=1.9943。C_cert约1.9943164615684859、U175约28.920677505176677仅展示。准确25系数见论文及 `../lower_round003/a2_equals_a3_results.json`；完整上界四次幂、根区间和差值见 `../round020/primary_upper175_acceptance.json`。不从展示小数重造证书。

暂停后现有上下界审计及代码/实际结果在 `../final_bounds020/`。新下界审计采用1122个互不相交盒，重建原42见证全部25coeff/原区间；明确finite collar与严格平滑余量，不提高下界。新上界审计核39项完整准确比较，并审175实际独立根/支付/全树及继承分析链。最终稿另经上下界分工审读，新短常数经独立整数解析通过。排版后的源能准确逆向还原审读锁定版本，见 `../round020/paper_format_equivalence.json`。

最小下界几何复核可把 `audit/012_portable_C_verify.py` 与 `audit/002_independent_lower_audit.py` 复制到仓库外新目录，再用前者 `--project-root /ABS/cube-dim-5-bounds`，其输出落在新目录，避免覆盖历史结果。它与本次全新disjoint-union审计不同，本次未重跑该旧入口。175原根证书及未来显式复核命令在 `../upper_secant175/ARCHIVE_EXECUTION.md`；本次参数等价归档没有执行未来数学入口。最终标量gate的冻结ledger快照恢复见 `../round020/replay_final_bounds_gate.py`（本次未再执行）。

完整g/H、时间端点与原有向区间资格仍有共享依赖；不同审计程序亦共享Python整数库，不声称重新生成第二套H引擎或形式化验证。173仅候选，独立数值审计从未开始，未入论文。研究目标保持用户要求的暂停状态；原双侧突破及后三轮条件未满足。
