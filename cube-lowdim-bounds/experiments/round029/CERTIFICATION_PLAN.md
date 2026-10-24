# 第29轮唯一候选全域核验预算

LP及准确成本已完成后、核验前冻结。候选beta=alpha=47/40，去掉仅元数据beta/structure，全部算术字段保持逐项相同，schema_binding.json绑定hash。

原版cert_run.py一次：Q8、maxdepth7、Kdepth6、procs6、only MKA，600秒墙钟；独立hl_check.py一次：claim3.59091、procs6，600秒墙钟。两个专属进程组可并行，超时停止，不加深、不重跑。成本已用独立整数根/Fraction严格包围，两全域实现共享Arb算术及已审公式，不能称完全独立。任何失败或超时不升级界。完整大JSON本地保留，仅摘要及hash入库。
