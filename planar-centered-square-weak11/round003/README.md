# 第3轮复现与证据范围

本轮区间沿用[第2轮精确端点](../round002/strict-interval.json)，没有新一般数值界。[本轮报告](../rounds/003.md)说明反例、有限第一变分与mixed见证必要性。

- [统一删除反例](uniform/proof.md) / [独审](uniform/review/audit.md)：全部正尺度、概率有理原子族及直接L¹函数。
- [固定有限族第一变分](theory/proof.md) / [独审](theory/audit.md)：不声称真实halo导数或统一原子数余项。
- [机制算核说明](mechanism/README.md)：两端精确匹配R1/R2证书；顺序归属不唯一。

从项目根运行：

```sh
python3 round003/uniform/check_family.py
python3 round003/uniform/review/probability_family_checks.py
wolframscript -file round003/uniform/review/independent_checks.wl
python3 round003/mechanism/decompose.py --output-dir /tmp/r003-reproduction
OPENBLAS_NUM_THREADS=1 python3 round003/mechanism/independent_verify.py --output-dir /tmp/r003-reproduction
python3 round003/mechanism/check_mechanism.py --output-dir /tmp/r003-reproduction
wolframscript -file round003/mechanism/scalar.wl /tmp/r003-reproduction
wolframscript -file round003/theory/two_atom_first_order.wl
```

Python需支持fractions和现代pathlib；独立坐标压缩面积核还需NumPy。最后一个直接质量核另检查纯新点与mixed分类，数学字段已在移位后实跑对比，见`mechanism/relocated-review-receipt.json`。若只核标量不重新生成记录，可将scalar.wl的目录参数指向本包mechanism/，其中已保存summary.json。

反例Python脚本按原始设计在同目录写sample或JSON，可重建已归档的有理证书。机制脚本把2.85MB全记录写入指定生成目录；该可再生文件、缓存、完整原研究scratch不纳入Git。`mechanism/SHA256SUMS`是保留的历史身份清单，包含原/移位全记录和带project/前缀的输入位置，不能作为当前目录的逐文件check命令；公开文件的可核SHA见项目verification-manifest.json及本轮source-bindings.json。

冻结正文及独审保留原文与原始来源路径，包含冻结当时“待审核”等历史状态；当前裁决以配对的最终audit为准。公开路径映射、便携化改动和省略记录的理由见[source-bindings.json](source-bindings.json)。证明与关键审文均字节保留；派生脚本只调整路径、输出位置及摘要，有限数学循环与精确比较不变。没有将哈希一致当作证明，也未声称全程CPU满载。
