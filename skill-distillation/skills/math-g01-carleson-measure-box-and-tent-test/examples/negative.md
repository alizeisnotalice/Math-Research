# 根测试不足反例

根测度为 1、左右子节点测度各 1/2、孙节点测度各 1/4。只在 LL 赋系数 α_LL=1/2，则根比值为 1/2，而 LL 子树比值为 2；真实全子树 packing 常数为 2。只检查根节点不能推出 Carleson 条件。

该有限例已实际计算，见 `audit_current/g/case_results_20261007.json`（G01-negative-root-only-test）。
