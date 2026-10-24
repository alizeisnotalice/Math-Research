# 第35轮固定点诊断与依赖更正

本轮只针对第34轮唯一候选，未作新LP或全域认证。历史原版first积分API的严格性说明已更正，详见dependency_correction.md；成功运行记录不予覆盖。

复现需Python3.12及python-flint0.9.0。在新的临时根目录建round01、round15、round16：第20轮scaled_certificate.json置round01；第34轮kernel_inputs.tar.xz解压到round15，并复制同轮first_report.json。当前diagnose_once.py/run_once.py放round16；第23轮diagnose.py复制为round16/dependency_round04_diagnose.py；第20轮upstream-v4/checker置round16/public_src，保留independent/hl_core.py子路径。各输入均由源码SHA绑定，不需要私人附件或旧大矩阵。

在全新输出目录执行 `python -B run_once.py`。硬预算60秒、固定192bit、两个点乘两实现共4个完整标量；不要复制已归档result/points/calls/run/log到运行目录，程序拒绝覆盖。该命令只是复现说明，本地唯一运行证据另行归档。

第一实现的外区积分采用本轮本地修正回调，所有复log传递analytic标志，未修改已归档上游。第二实现是既有闭式超几何。两者共享Arb及解析公式，不是完全独立底层算术。符号根据重读保存区间的严格端点判定；一次点评价不构成全域正性证书。

本轮没有补验旧first全域。下一轮先固定B1制作最小补丁及有预算的迁移核验，恢复其第二条证据链，再继续优化。
