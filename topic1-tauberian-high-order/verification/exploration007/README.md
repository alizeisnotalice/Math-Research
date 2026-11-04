# 主孔—卫星模型的可移植输入与结果

`satellite.py` 相对于自身所在目录创建 `satellite/` 输出，使用 Python、NumPy、SciPy；默认每二维分层256点、每三维分层32点，15进程，运行无时限。`kappa_calibration.py` 使用同一目录的 `satellite.py` 与 Python 标准库 Decimal，结果写入脚本目录下 `satellite/`。

`satellite_summary.json` 是实跑135组合摘要：126组合满足前提并完成、9组合卫星见证前提未满足。每个有效组合包含输入、分层体积与计数、Gamma估计及一个实际采样主点见证（若有）；卫星中心的缺陷和严格余量单列。采样仅限文档中的主孔半径 [δ/u,1]；全部为浮点数值结果，QMC和浮点误差无认证。

`kappa_calibration_result.json` 保存κ=2特殊校准、80位Decimal预算根及样本结果。原完整分支记录保存在上级 `satellite/branch_*.json`；本目录不含大日志。数值结果不能代替全部捕获的解析证明。

精确预算复现见 ../wolfram/satellite-budget.wl。固定输入仅有浮点见证核验，未认证QMC或舍入误差。暂停前未启动后续批次。
