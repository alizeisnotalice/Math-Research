# 第40轮最小复现

仅固定node661，2次192bit Hab积分，复用旧中心T1和AS。实际0.146398833秒；16新下界全正，但不是完整父格或新的弱型界。

临时根放round01/scaled_certificate.json（第20轮），round17/fixed_checker/hlcert.py（第36轮补丁重建），round20/leaf_cover.json及models.local/node_00661.json（第39轮saved_models.tar.xz），本目录程序放round21。锁定输入SHA由源码核对。使用Python3.12/python-flint0.9.0；空输出目录执行`python -B run_once.py`，60秒硬限，已有结果拒绝重跑。

只审保存算术时复制本目录保存输出及前述来源，然后在尚无saved_verification.json的副本运行`python -B verify_saved.py`。其90秒预算仅标准库，核两区间绝对上端、16行、来源与4调用/8事件，不重算积分。实际0.008246416秒通过。旧模型、原始源码和私人附件不重复上传。

analytic_audit和proofs/exterior_signed_taylor.md说明单侧公式；不能将该扣项用于负核系数或当作双侧误差。历史第39轮cover保留未决身份不回写，未来新cover须显式引用本次补充来源。
