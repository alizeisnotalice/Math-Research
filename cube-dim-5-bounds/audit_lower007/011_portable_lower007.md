# 011 — 可携 lower007 归档复核入口

`011_portable_lower007.py` 只接受一个项目根目录，并仅读取：

- `lower/certify_rectangular.py`；
- `lower_round003/a2_equals_a3_results.json`；
- `lower_round007/common_tail_event_intervals.py`（用于锁定生产来源hash）；
- 八个 `lower_round007/exactC_event_interval_*_results.json`。

它不读取生产clone的绝对路径、不需要48个独立sample JSON，也不导入生产事件、插值或Bernstein实现。八个final JSON已内嵌所需48个准确并集证书。

实际从只含上述11个输入文件的临时独立staging树运行。`011_portable_lower007_receipt.json`保留本次实际结果、完整输入SHA-256、175个原始端点、369个事件和八段通过记录。临时路径不是可携接口的一部分。

重放：

```sh
python -B audit_lower007/011_portable_lower007.py . \
  --output audit_lower007/011_portable_lower007_receipt.json
```

它仍共享项目中归档的 \(\mathbb Q(u,v)\) 算术和隔离器，且复用final JSON内的准确并集证书；它不重跑48次并集基线。
