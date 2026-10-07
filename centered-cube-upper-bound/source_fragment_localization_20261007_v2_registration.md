# 开区间端点检查修复：单独预登记 v2

2026-10-07；登记先于 v2 执行。原 guard 与 failure.json 保留，原 failed run 不改成 PASS。

修复仅将远离族 R1 的端点检查 left−d/2>1/2 改为 ≥1/2。原 E_R=(9/16,15/16) 是开 interval，其每个内部点仍有严格距离；闭端点至多与 source 面接触，L¹ 捕获质量仍为0。因此不用改变原输入、尺度、阈值或预定几何结论。

三轮仍 J=12,24,48；其余全部断言原样。新 guard_v2.py 单独执行一次，保存 v2_results/v2_receipt，并核本登记 SHA256。此登记继承原 registration 的范围约束：只有一维原有限 winner 的空间组件证书，非一般 order 或原 soft gate 证据。
