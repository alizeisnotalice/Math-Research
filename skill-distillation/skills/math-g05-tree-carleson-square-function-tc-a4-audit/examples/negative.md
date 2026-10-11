# 根 packing 不能代替逐子树检查

深度 1 二进树根测度 1、左子节点测度 1/2，仅在左节点赋 α_L=1。根比值为 1，左子树比值为 2，故真实全子树 packing 常数是 2。只核根节点不足以调用 Carleson 嵌入。

本例已实际运行，见 `audit_current/g/case_results_20261007.json`（G05-negative-root-only-packing-gate）。TC-A4 的私有规范仍未定义。
